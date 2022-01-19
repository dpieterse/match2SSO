#!/usr/bin/env python
# coding: utf-8

# # match2SSO
# * Running instructions are given at the start of function run_match2SSO
# * Originally written as a Jupyter Notebook using Python 3.8.10
# * Compatible with BlackBOX / ZOGY version 1.0.0 and up
# * Compatible with MeerLICHT / BlackGEM observations
# 
# Output (SSO) catalogue columns and header keywords are listed here:
# https://www.overleaf.com/read/zrhqwcbkfqns
# 
# <i>match2SSO</i> makes grateful use of the <i>lunar</i> and <i>jpl_eph</i> 
# repositories that were written by Bill Gray under Project Pluto. The core of 
# <i>match2SSO</i> is <i>astcheck</i>: a C++ script in the <i>lunar</i> 
# repository that matches detections to known solar system objects. More 
# information on <i>astcheck</i> can be found at: 
# https://www.projectpluto.com/astcheck.htm
# 
# <b>Dependencies on scripts and files</b>
# * <i>lunar</i> package (https://github.com/Bill-Gray/lunar)
# * <i>jpl_eph</i> package (https://github.com/Bill-Gray/jpl_eph)
# * JPL DE ephemeris file (ftp://ssd.jpl.nasa.gov/pub/eph/planets/Linux/)
# * MPC's Observatory codes list 
#   (https://www.minorplanetcenter.net/iau/lists/ObsCodes.html)
# 
# In addition, match2SSO uses MPCORB.DAT (MPC's asteroid database) and 
# COMET.ELEMENTS (JPL's comet database), but these are downloaded when running 
# the script and hence do not need to be pre-downloaded.
# 
# if keep_tmp is False, all temporary files and folders are removed except for 
# the most recent, unintegrated version of the asteroid & comet databases in
# the databaseFolder.
# For the night mode, we'll need to implement the cleaning function in 
# BlackBOX that is commented out at the end of the night mode section in 
# run_match2SSO. (As the night mode will be run in parallel on different
# catalogues and the individual night mode processes don't communicate.)

# ## Python packages and settings

# In[ ]:


import os
from pathlib import Path

import set_match2SSO as settingsFile #load match2SSO settings file

import numpy as np
import pandas as pd
import glob, re
from string import ascii_lowercase, ascii_uppercase

import argparse #for parsing command line arguments when starting match2SSO
import requests #for downloading databases
import subprocess #for running command line commands / C++ scripts
import shutil #to remove temporary folders

import time
from datetime import datetime, timedelta
import pytz
from pytz import timezone

from astropy.time import Time
from astropy.io import fits
from astropy.table import Table, Column
from astropy import units as u
from astropy.coordinates import SkyCoord

from multiprocessing import Pool, Manager, Lock, Queue, Array
import sys, resource, psutil, platform #for memory usage logging
import logging

#Set up log
logFormat = ('%(asctime)s.%(msecs)03d [%(levelname)s, %(process)s] %(message)s '
          '[%(funcName)s, line %(lineno)d]')
dateFormat = '%Y-%m-%dT%H:%M:%S'
logging.basicConfig(level='INFO', format=logFormat, datefmt=dateFormat)
logFormatter = logging.Formatter(logFormat, dateFormat)
logging.Formatter.converter = time.gmtime #convert time in logger to UTC
log = logging.getLogger()


# In[ ]:


#Set version
__version__ = "1.0.0"
keywords_version = '1.0.0'

#Relevant transient catalogue column names
numberColumn = 'NUMBER'
RA_column = 'RA_PSF_D'   #ra in deg
DEC_column = 'DEC_PSF_D' #dec in deg
magnitudeColumn = 'MAG_ZOGY'
dummyColumn = 'TDUMCAT'


# ## Main functions to run match2SSO

# In[ ]:


def run_match2SSO(tel, mode, cat2process, date2process, list2process,
                  logName, keep_tmp, redownloadDatabases, includeComets,
                  overwriteFiles, timeFunctions):
    """
    Run match2SSO on the input catalogue(s)/date. match2SSO can be run in 
    different mode / date2process / cat2process / list2process combinations.
    Allowed combinations are: (if not mentioned, the variable is None)
    # * Day mode
    # * Day mode + date2process
    # * Night mode + cat2process
    # * Historic mode + cat2process
    # * Historic mode + date2process
    # * Historic mode + list2process
    
    Day mode: 
    Create known objects database & CHK files for the upcoming night, or 
    - in case date2process is specified - for the specified night.
    0) Creates a run directory in preparation of the nightly processing
    1) Downloads asteroid and comet databases
    2) Integrates the asteroid database to midnight of the observation night
    3) Combines the comet and integrated asteroid databases into a SOF-formatted
       known objects database. 
    4) Runs astcheck on a fake detection in order to create the CHK files that
       astcheck will need for faster processing when running on observations. 
    5) Creates symbolic links to the used databases and the observatory codes
       list in the run directory.
    (Products of steps 1-2 are saved to the databaseFolder, those of steps 3-5
    to the run directory.)
    
    Night mode: 
    Run match2SSO on a single transient catalogue. The day mode should have
    been run once before the night mode. This allows the night mode to run in
    parallel on multiple transient catalogues of the same night, as steps that 
    cannot be parallelised (making known objects database and CHK files) have
    already been executed in the day mode. The night mode:
    1) Converts the transient catalogue into an MPC-formatted text file.
    2) Runs astcheck on that file, to find matches between the transient 
       detections and known solar system objects. 
    3) Makes an SSO catalogue containing the matches.
    4) Makes an MPC submission file of the matches.
    
    Historic mode: 
    The historic mode does the entire processing - executing steps 0-3 & 5 of 
    the day mode, followed by all steps of the night mode. It can be run on a 
    single transient catalogue, an entire night of observations or a list of 
    observations (possibly spanning multiple nights). For the first catalogue 
    that is processed of each observation night, a new known objects catalogue
    is created that is integrated to the observation midnight. The asteroid and
    comet databases used for this are only downloaded once per historic mode run
    (and only if redownloadDatabases is True). For the remaining files (or if 
    redownloadDatabases is False), the most recently downloaded versions are used. 
    
    Parameters:
    -----------
    tel: string
            Abbreviated telescope name. Can be either ML1, BG2, BG3 or BG4.
    mode: string
            Mode in which match2SSO is run. This can be 'day', 'night' or
            'historic'.
    cat2process: string
            Path to and name of the transient catalogue that is to be processed.
    date2process: string
            Formatted as yyyymmdd, yyyy-mm-dd, yyyy/mm/dd or yyyy.mm.dd. When 
            used with day mode: date for which the known objects database needs
            to be prepared. When used with historic mode: date for which all 
            light transient catalogues need to be processed.
    list2process: string
            Path to and name of the text file that contains the paths to and 
            names of transient catalogues (one per line) that need to be 
            processed.
    logName: string
            Path to and name of the log file in which comments about the run are
            stored.
    keep_tmp: string
            String that can be converted to a boolean with str2bool(). This
            boolean indicates whether the temporary files made during the 
            processing should be kept or removed at the end of the processing.
    redownloadDatabases: string
            String that can be converted to a boolean with str2bool(). This
            boolean indicates whether the asteroid and comet databases will need
            to be redownloaded when making the known objects database. 
            Alternatively, the most recently downloaded version of the databases
            are used.
    includeComets: string
            String that can be converted to a boolean with str2bool(). This
            boolean indicates whether comets should be included in the known
            objects database. There have been issues with matching to comets
            with large orbital uncertainties in astcheck, so as long as this has
            not been solved, comet matching should be avoided. (time of writing:
            23 Dec 2021)
    overwriteFiles: string
            String that can be converted to a boolean with str2bool(). This
            boolean indicates whether files are allowed to be overwritten.
    timeFunctions: string
            String that can be converted to a boolean with str2bool(). This
            boolean indicates whether functions need to be (wall-)timed.
    """
    t_glob = time.time()
    
    #Format input parameters
    mode = mode.lower()
    
    keep_tmp = str2bool(keep_tmp)
    redownloadDatabases = str2bool(redownloadDatabases)
    includeComets = str2bool(includeComets)
    overwriteFiles = str2bool(overwriteFiles)
    timeFunctions = str2bool(timeFunctions)
    localTimeZone = timezone(get_par(settingsFile.timeZoneTelescope, tel))
    if date2process != None:
        date2process = date2process.replace(".","").replace("/","").replace("-",
                                                                            "")
    #Perform checks on input parameter combinations and setting file parameters
    if not checkInputParameters(mode, cat2process, date2process, list2process):
        return
    if not checkSettingsFile(tel):
        return
    
    #Logging
    setUpLogFile(logName)
    log.info("Mode: {}".format(mode))
    
    mem_use(label='at start of run_match2SSO')
    
    #Get local noon corresponding to the night start in case date2process or
    #cat2process are specified. nightStart is a datetime object (incl. timezone
    #information).
    if cat2process!=None:
        nightStart = getNightStartFromCatalogueName(cat2process, tel)
        
    elif date2process!=None:
        nightStart = localTimeZone.localize(datetime.strptime(date2process +
                                                    " 120000", "%Y%m%d %H%M%S"))
    
    if mode=="day":
        
        log.info("Running the day mode.")
        
        #If no observation night is specified, use the upcoming local night
        if date2process==None:
            nightStart = (datetime.now(localTimeZone)).strftime("%Y%m%d 120000")
            #Add timezone info to datetime object
            nightStart = localTimeZone.localize(datetime.strptime(nightStart,
                                                               "%Y%m%d %H%M%S"))
        
        #Create a run directory corresponding to the observation night
        runDirectory = ("{}{}/"
                        .format(databaseFolder, nightStart.strftime("%Y%m%d")))
        log.info("Run directory: {}".format(runDirectory))
        if not os.path.isdir(runDirectory):
            os.makedirs(runDirectory)
        
        #Create symbolic link to the observatory codes list
        if not os.path.exists("{}ObsCodes.html".format(runDirectory)):
            log.info("Creating symbolic link to ObsCodes.html")
            os.symlink("{}ObsCodes.html".format(softwareFolder),
                       "{}ObsCodes.html".format(runDirectory))
        
        #Download and integrate known object databases
        midnight = nightStart + timedelta(days=0.5)
        createKnownObjectsDatabase(midnight, runDirectory, redownloadDatabases,
                                   includeComets, keep_tmp, timeFunctions)
        
        def create_CHKfiles(noonType):
            
            """
            Function that creates the CHK files that astcheck uses (and produces
            if they don't exist yet) when matching. These files describe the
            positions of all asteroids at the start and end of the night (both
            at noon) in UTC. By running astcheck on a fake detection, we can
            produce these CHK files in advance (which allows parallelisation of
            match2SSO runs in the night mode). Of course we will subsequently
            remove the fake detection and fake matches.
            
            As the 24-hour local observing night (between local noons) can
            overlap with two UTC nights (between UTC noons) if there is a large
            time difference between the timezone of the telescope and UTC, we'll
            have to try producing CHK files for a time close to the start of the
            local night (we take 1 min after) as well as a time close to the end
            of it (1 min before). This function should hence be run in for both
            noon types ("noonstart" and "noonend"). This will produce a total of
            2 or 3 CHK files that astcheck will subsequently use to match any 
            observation taken during the local observing night to the known 
            solar system objects.
            """
            
            #Check noonType parameter and set observation time for fake 
            #detection
            noonType = noonType.lower()
            if noonType not in ["nightstart", "nightend"]:
                log.error("Unknown noon type!")
                return
            if noonType == "nightstart":
                obstime = nightStart + timedelta(minutes=1)
            elif noonType == "nightend":
                obstime = nightStart + timedelta(days=1) - timedelta(minutes=1)
            
            #Convert observation time of fake detection to UTC
            obstime = obstime.astimezone(pytz.utc)
            
            #Create MPC-formatted file with fake detection
            MPCfileName_fake = ("{}fakedetection_{}_MPCformat.txt"
                                .format(runDirectory, noonType))
            log.info("Creating fake detection: {}".format(MPCfileName_fake))
            MPCfile_fake = open(MPCfileName_fake, "w")
            fakeDetection = ("     0000001  C{} {:0>2} {:08.5f} "
                            .format(obstime.year, obstime.month, obstime.day)
                           + "00 00 00.00 +00 00 00.0          0.00 G      L66")
            MPCfile_fake.write(fakeDetection)
            MPCfile_fake.close()
            
            #Run astcheck on fake observation to create CHK files
            log.info("Running astcheck on fake detection")
            astcheckOutputFileName_fake = MPCfileName_fake.replace("_MPCformat.txt", 
                                                             "_astcheckMatches.txt")
            runAstcheck(MPCfileName_fake, runDirectory, astcheckOutputFileName_fake,
                        timeFunctions, overwriteFiles, matchingRadius=0)

            #Remove MPC-formatted file and astcheck output related to the fake 
            #detection
            os.remove(MPCfileName_fake)
            log.info("Removed {}".format(MPCfileName_fake))
            os.remove(astcheckOutputFileName_fake)
            log.info("Removed {}".format(astcheckOutputFileName_fake))
            
            return
        
        #Create CHK files that astcheck needs in advance, to allow 
        #parallelisation
        create_CHKfiles("nightstart")
        create_CHKfiles("nightend")
        
        #Check for known object database products
        if not checkForDatabaseProducts(runDirectory):
            log_timing_memory(t_glob, label='run_match2SSO')
            logging.shutdown()
            return
        
        log.info("Day mode finished.")
    
    
    elif mode == "night":
        
        log.info("Running the night mode on transient catalogue: \n{}" 
                 .format(cat2process))
        
        runDirectory = "{}{}/".format(databaseFolder,
                                      nightStart.strftime("%Y%m%d"))
        log.info("Run directory: {}".format(runDirectory))
            
        #Check for known object database products. Stop processing if it doesn't
        #exist
        if not checkForDatabaseProducts(runDirectory):
            logging.shutdown()
            return
        
        #Check for CHK files
        utcNightStart = getNightStartFromCatalogueName(cat2process, tel, "utc")
        utcNightEnd = utcNightStart + timedelta(days=1)
        utcNightStart = utcNightStart.strftime("%Y%m%d")
        utcNightEnd = utcNightEnd.strftime("%Y%m%d")
        if not os.path.exists("{}{}.chk".format(runDirectory, utcNightStart)):
            log.critical("Missing {}.chk!".format(utcNightStart))
            logging.shutdown()
            return
        if not os.path.exists("{}{}.chk".format(runDirectory, utcNightEnd)):
            log.critical("Missing {}.chk!".format(utcNightEnd))
            logging.shutdown()
            return
        del utcNightStart
        del utcNightEnd
        
        #Check for observatory codes list. Stop processing if it doesn't exist.
        if not os.path.exists("{}ObsCodes.html".format(runDirectory)):
            log.critical("{}ObsCodes.html doesn't exist.".format(runDirectory))
            logging.shutdown()
            return
        
        _ = matchSingleCatalogue(cat2process, runDirectory, nightStart, 
                                 makeKOD=False, redownloadDatabases=False,
                                 includeComets=includeComets, keep_tmp=keep_tmp,
                                 timeFunctions=timeFunctions,
                                 overwriteFiles=overwriteFiles)
        
        #Beware that the run directory created for the processing of the 
        #catalogue is not removed. This is the case because a single parallel
        #process does not know about the rest. A cleaning function should be
        #run at the end of the nightly processing if one wants to remove the
        #runDirectory. Also delete the integrated asteroid database that was
        #made during the day mode at this time. See code below.
        #if not keep_tmp:
        #    if os.path.exists("{}MPCORB.DAT".format(runDirectory)):
        #        asteroidDBname = os.readlink("{}MPCORB.DAT".format(runDirectory))
        #        if "epoch" in asteroidDBname:
        #            os.remove(asteroidDBname)
        #            log.info("Removed {}".format(asteroidDBname))
        #    removeTemporaryFolder(runDirectory)
    
    elif mode=="historic":
        
        def matchCataloguesSingleNight(cataloguesSingleNight, startNight,
                                       redownloadDB):
            """
            Process input catalogues corresponding to observations taken on the
            same night.
            """
            if timeFunctions:
                t0 = time.time()
            
            log.info("{} catalogues to process for the night around {}."
                    .format(len(cataloguesSingleNight), 
                            startNight.strftime("%Y-%m-%d %H:%M:%S")))
            
            #Create run directory
            runDirectory = "{}{}/".format(databaseFolder,
                                          startNight.strftime("%Y%m%d"))
            log.info("Run directory: {}".format(runDirectory))
            if not os.path.exists(runDirectory):
                os.makedirs(runDirectory)
            
            #Run matching per catalogue
            makeKOD = True
            for index, catalogueName in enumerate(cataloguesSingleNight):
                log.info("Processing {}".format(catalogueName))
                
                madeKOD = matchSingleCatalogue(catalogueName, runDirectory,
                                               startNight, makeKOD,redownloadDB,
                                               includeComets, keep_tmp, 
                                               timeFunctions, overwriteFiles)
                if madeKOD:
                    makeKOD = False #Only make known objects database once
            
            #Remove the run directory after processing the last catalogue of the
            #night
            if not keep_tmp:
                #Remove integrated database made for this night
                if os.path.exists("{}MPCORB.DAT".format(runDirectory)):
                    asteroidDBname = os.readlink("{}MPCORB.DAT"
                                                 .format(runDirectory))
                    if "epoch" in asteroidDBname:
                        os.remove(asteroidDBname)
                        log.info("Removed {}".format(asteroidDBname))
                
                #Remove temporary folder made for this night
                removeTemporaryFolder(runDirectory)
            
            if timeFunctions:
                log_timing_memory(t0, label='matchCataloguesSingleNight')
            
            return
        
        
        if cat2process!=None:
            log.info("Running historic mode on transient catalogue: \n{}"
                     .format(cat2process))
            matchCataloguesSingleNight([cat2process], nightStart, 
                                       redownloadDatabases)
        
        elif date2process!=None:
            log.info("Running historic mode on night {}".format(date2process))
            catalogues2process = getTransientFileNames(nightStart, nightStart +
                                                       timedelta(days=1), tel,
                                                       timeFunctions)
            if len(catalogues2process)==0:
                log.critical("No light transient catalogues exist for night {}"
                             .format(date2process))
                log_timing_memory(t_glob, label='run_match2SSO')
                logging.shutdown()
                return
            
            matchCataloguesSingleNight(catalogues2process, nightStart,
                                       redownloadDatabases)
        
        elif list2process!=None:
            log.info("Running historic mode on catalogue list: \n{}"
                     .format(list2process))
            with open(list2process, "r") as catalogueList:
                listedCatalogues = [name.strip() for name in catalogueList 
                                    if name[0]!='#']
            
            #Order by observation date (noon that equals the start of the 
            #observation day)
            noons = []
            catalogues2process = []
            for catalogueName in listedCatalogues:
                noon = getNightStartFromCatalogueName(catalogueName, tel)
                noons.append(noon.strftime("%Y%m%d %H%M%S"))
                catalogues2process.append(catalogueName)
            
            #Process files per night
            log.info("Catalogue list spans {} nights"
                     .format(len(np.unique(noons))))
            firstNight = True
            for noon in np.unique(noons):
                log.info("Processing night that starts at {}" .format(noon))
                nightIndex = np.where(np.array(noons)==noon)[0]
                catalogues2process_1night = np.array(catalogues2process
                                                    )[nightIndex]
                
                #Use the same version of the asteroid and comet databases for
                #all data to be processed. If redownloadDatabases, this version
                #corresponds to the time of processing the first image from the
                #list. If redownloadDatabases is Fale, the latest existing 
                #downloaded version of the databases is used.
                noon = localTimeZone.localize(datetime.strptime(noon, 
                                                               "%Y%m%d %H%M%S"))
                if firstNight:
                    matchCataloguesSingleNight(catalogues2process_1night, 
                                               noon, redownloadDatabases)
                else:
                    matchCataloguesSingleNight(catalogues2process_1night, 
                                               noon, redownloadDB=False)
                firstNight = False
    
    
    log.info("Finished running match2SSO.")
    log_timing_memory(t_glob, label='run_match2SSO')
    logging.shutdown()
    
    return


# In[ ]:


def matchSingleCatalogue(catalogueName, runDirectory, nightStart, makeKOD, 
                         redownloadDatabases, includeComets, keep_tmp, 
                         timeFunctions, overwriteFiles):
    """
    Run matching routine on a single transient catalogue. Optionally, a new 
    known objects database is created where the reference epoch corresponds to 
    midnight on the observation night. The detections in the transient catalogue
    are then matched to the positions of the solar system bodies in the known
    objects catalogue. Matches are saved to an SSO catalogue.
    Function returns a boolean indicating whether a known objects catalogue was
    (successfully) made.
    
    Parameters:
    -----------
    catalogueName: string
            Name of the transient catalogue of which the detections are to be
            matched to known solar system objects.
    runDirectory: string
            Directory corresponding to observation night where all temporary
            products will be stored during the running of match2SSO.
    nightStart: datetime object, including time zone
            Noon corresponding to the start of the local night during which the
            observation corresponding to the transient catalogue was made.
    makeKOD: boolean
            Boolean indicating whether a new known objects database needs to
            be made. In night mode, this should be False. In historic mode, 
            this is only True for the first catalogue that is processed per
            observation night, as the database will need to be integrated to
            that observation night.
    redownloadDatabases: boolean
            Only used when makeKOD is True. This boolean indicates whether the
            asteroid and comet databases will need to be redownloaded before
            making the known objects database. Alternatively, the most recent,
            previously downloaded version of the databases are used.
    includeComets: boolean
            Boolean indicating whether comets are included in the known objects
            catalogue.
    keep_tmp: boolean
            Boolean indicating whether the temporary files corresponding to
            the transient catalogue that are made in this function should be 
            kept or removed at the end of the processing.
    timeFunctions: boolean
            Boolean indicating whether functions need to be (wall-)timed.
    overwriteFiles: boolean
            Boolean indicating whether files are allowed to be overwritten.
    """
    mem_use(label='at start of matchSingleCatalogue')
    if timeFunctions:
        t0 = time.time()
    
    madeKOD = False
    
    #Check if input catalogue exists and is not flagged red
    isExisting, isDummy, catalogueName = checkInputCatalogue(catalogueName)
    if not isExisting:
        return madeKOD
    del isExisting
    
    if isDummy:
        SSOcatalogueName = catalogueName.replace("_light", "").replace(".fits",
                                                                   "_sso.fits")
        if os.path.exists(SSOcatalogueName) and not overwriteFiles:
            log.warning("{} already exists and will not be re-made."
                        .format(SSOcatalogueName))
            return makeKOD
        
        create_SSOcatalogue(None, runDirectory, SSOcatalogueName, includeComets,
                            timeFunctions, overwriteFiles)
        return madeKOD
    del isDummy
    
    #If makeKOD, create a new known objects database with a reference epoch 
    #corresponding to midnight of the observation night. Also create a symbolic
    #link to the MPC observatory codes list.
    if makeKOD:
        midnight = nightStart + timedelta(days=0.5)
        createKnownObjectsDatabase(midnight, runDirectory, redownloadDatabases,
                                   includeComets, keep_tmp, timeFunctions)
        
        #Check if the run directory contains the proper known objects database
        #files for further processing
        if not checkForDatabaseProducts(runDirectory):
            return madeKOD
        
        madeKOD = True
        
        #Make symbolic link to observatory codes list if it doesn't exist yet
        if not os.path.exists("{}ObsCodes.html".format(runDirectory)):
            log.info("Creating symbolic link to ObsCodes.html")
            os.symlink("{}ObsCodes.html".format(softwareFolder),
                       "{}ObsCodes.html".format(runDirectory))
        
    #Convert the transient catalogue to an MPC-formatted text file
    MPCfileName = "{}{}".format(runDirectory, os.path.basename(catalogueName
                      ).replace("_light", "").replace(".fits","_MPCformat.txt"))
    MPC_code = convertCatalogue2MPCformat(catalogueName, MPCfileName, 
                                          overwriteFiles, timeFunctions)
    if MPC_code == None:
        log.critical("Stop running match2SSO on catalogue because of unknown "
                     +"MPC code.")
        return madeKOD
    
    #Run astcheck on the MPC-formatted transient file
    astcheckOutputFileName = MPCfileName.replace("_MPCformat.txt",
                                                 "_astcheckMatches.txt")
    runAstcheck(MPCfileName, runDirectory, astcheckOutputFileName, timeFunctions,
                overwriteFiles)
    
    #Save matches found by astcheck to an SSO catalogue
    SSOcatalogueName = catalogueName.replace("_light", "").replace(".fits",
                                                                   "_sso.fits")
    create_SSOcatalogue(astcheckOutputFileName, runDirectory, SSOcatalogueName,
                        includeComets, timeFunctions, overwriteFiles)
    
    #Create a submission file that can be used to submit the detections that 
    #were matched to known solar system objects to the MPC
    createSubmissionFile(SSOcatalogueName, MPCfileName, MPC_code, timeFunctions,
                         overwriteFiles)
    
    #Delete temporary files corresponding to the processed transient catalogue. 
    #The other temporary files (the CHK files, the SOF file and the symbolic 
    #links) in the runDirectory are not (yet) removed, as they might be needed 
    #for processing of other data from the same night.
    if not keep_tmp:
        os.remove(MPCfileName)
        log.info("Removed {}".format(MPCfileName))
        os.remove(astcheckOutputFileName)
        log.info("Removed {}".format(astcheckOutputFileName))
    
    if timeFunctions:
        log_timing_memory(t0, label='matchSingleCatalogue')
    
    return madeKOD


# ## Core functions in match2SSO
# * Create known objects catalogue (with asteroids with a max. orbital uncertainty)
# * Convert transient catalogue to MPC format
# * Run astcheck
# * Save astcheck matches to SSO catalogue
# * Create submission file

# In[ ]:


def createKnownObjectsDatabase(midnight, runDirectory, redownloadDatabases,
                               includeComets, keep_tmp, timeFunctions):
    """
    Function downloads the most recent versions of the asteroid database and the
    comet database. It then uses integrat.cpp from the lunar repository to 
    integrate the asteroid orbits to midnight of the observation night, in order
    to optimize the predicted positions of known objects. 
    
    Beware: the current version of integrat.cpp cannot be used on JPL's comet 
    file, as it is not compatible with its format. As a consequence, there might
    be an offset in the predictions of the comet positions, perhaps causing us 
    to miss these objects in the linking routine. Note that the MPC's comet file
    could be integrated to the right epoch using integrat.cpp, but as this file
    is not compatible with astcheck, it cannot be used here.
    
    Parameters:
    -----------
    midnight: datetime object, including time zone
            Local midnight during the observation night.
    runDirectory: string
            Directory in which mpc2sof is run and which the known objects
            catalogue is saved to.
    redownloadDatabases: boolean
            If False, the databases will not be redownloaded. They will only be 
            integrated to the observation epoch (midnight on the observation
            night).
    includeComets: boolean
            Boolean indicating whether comets are to be included in the known 
            objects database. As of 8 dec 2021, there are issues with matching
            to known comets with large orbital uncertainties, so this boolean
            should be False until these issues can be fixed in astcheck.
    keep_tmp: boolean
            Only relevant if a new database version is downloaded. If keep_tmp 
            is False, the old versions of the database will be removed.
    timeFunctions: boolean
            Boolean indicating whether functions need to be (wall-)timed.
    """
    mem_use(label='at start of createKnownObjectsDatabase')
    if timeFunctions:
        t0 = time.time()
    
    def downloadDatabase(SSOtype):
        
        """
        SSOtype: string
                Solar system object type that is in the database. Can be either 
                'asteroid' or 'comet'. (Capitals are allowed as well.)
        """
        if timeFunctions:
            t1 = time.time()
        
        SSOtype = SSOtype.lower()
        if SSOtype == "asteroid":
            databaseURL = settingsFile.URL_asteroidDatabase
        elif SSOtype == "comet":
            databaseURL = settingsFile.URL_cometDatabase
        else:
            errorString = "Database type unknown. Cannot be downloaded."
            log.critical(errorString)
            raise ValueError(errorString)
        
        #Determine whether database needs to be downloaded
        existingDatabases = glob.glob("{}{}DB_*.dat".format(databaseFolder,
                                                            SSOtype))
        existingUnintegratedDatabases = [DB for DB in existingDatabases 
                                         if "epoch" not in DB]
        download = True
        if not redownloadDatabases and len(existingDatabases)>0:
            download = False
        
        #Download database if desired and get database version
        if download:
            log.info("Downloading {} database...".format(SSOtype))
            databaseVersion = datetime.utcnow().strftime("%Y%m%dT%H%M")
            databaseName = "{}{}DB_version{}.dat".format(databaseFolder,SSOtype,
                                                         databaseVersion)
            req = requests.get(databaseURL, allow_redirects=True)
            open(databaseName, "wb").write(req.content)
            log.info("{} database version: {}".format(SSOtype, databaseVersion))
            
            #Remove asteroids with large orbital uncertainties from database
            if SSOtype == "asteroid":
                selectAsteroidsOnUncertainty(databaseName, timeFunctions)
            
            if not keep_tmp and len(existingDatabases)>0:
                log.info("Removing older {} database versions.".format(SSOtype))
                for oldDatabase in existingDatabases:
                    os.remove(oldDatabase)
                    log.info('Removed {}.'.format(oldDatabase))
        else:
            #Retrieve most recent (unintegrated) database version. If there is
            #no unintegrated database, retrieve the most recent integrated one.
            #Empty databases (created when includeComets=False) are not taken 
            #into account, as these are in a different folder.
            databaseFiles = np.sort(existingUnintegratedDatabases)
            if len(databaseFiles)==0:
                databaseFiles = np.sort(existingDatabases)
            databaseName = databaseFiles[-1]
            databaseVersion = os.path.splitext(os.path.basename(databaseName)
                                        )[0].split("_")[1].replace("version","")
            log.info("{} database version: {}".format(SSOtype, databaseVersion))
            
        if timeFunctions:
            log_timing_memory(t1, label='downloadDatabase ({})'
                              .format(SSOtype))
        
        return databaseName, databaseVersion
    
    #Download asteroid database
    asteroidDatabaseName, asteroidDatabaseVersion = downloadDatabase("asteroid")
    
    #Download comet database if requested
    if includeComets:
        _, cometDatabaseVersion = downloadDatabase("comet")
    else:
        log.info("Do not download comet database. Instead, create an empty "
                 +"comet database so that there's no matching to comets.")
        cometDBheader = ("Num  Name                                     Epoch "+
                         "     q           e        i         w        Node   "+
                         "       Tp       Ref\n-------------------------------"+
                         "------------ ------- ----------- ---------- --------"+
                         "- --------- --------- -------------- ------------")
        open("{}ELEMENTS.COMET".format(runDirectory), "w").write(cometDBheader)
    
    #Integrat only accepts UTC midnights. Choose the one closest to local 
    #midnight.
    date_midnight = midnight.date()
    if midnight.hour >= 12.:
        date_midnight = date_midnight + timedelta(days=1)
    midnight_utc = pytz.utc.localize(datetime.strptime(date_midnight.strftime(
                                         "%Y%m%d")+ " 000000", "%Y%m%d %H%M%S"))
    del midnight
    
    #Integrate the asteroid database to the observation date.
    midnight_utc_str = midnight_utc.strftime("%Y%m%dT%H%M")
    filenameIntegratedAsteroidDB = ("{}asteroidDB_version{}_epoch{}.dat"
                                    .format(databaseFolder,
                                            asteroidDatabaseVersion,
                                            midnight_utc_str))
    if not os.path.exists(filenameIntegratedAsteroidDB):
        log.info("Integrating asteroid database to epoch {}..."
                 .format(midnight_utc_str))
        if timeFunctions:
            t1 = time.time()
        subprocess.run(["integrat", asteroidDatabaseName, 
                        filenameIntegratedAsteroidDB,str(Time(midnight_utc).jd),
                        "-f{}".format(settingsFile.JPL_ephemerisFile)],
                       cwd=databaseFolder)
        if timeFunctions:
            log_timing_memory(t1, label='integrat')
        
        #Remove temporary file created by integrat
        if os.path.exists("{}ickywax.ugh".format(databaseFolder)):
            os.remove("{}ickywax.ugh".format(databaseFolder))
    
    #Create the symbolic links in the run directory that mpc2sof needs
    runNameAsteroidDatabase = "{}MPCORB.DAT".format(runDirectory)
    if os.path.exists(runNameAsteroidDatabase):
        log.info("Removing the old MPCORB.DAT symbolic link")
        os.unlink(runNameAsteroidDatabase)
    os.symlink(filenameIntegratedAsteroidDB, runNameAsteroidDatabase)
    log.info("Created symbolic link {}".format(runNameAsteroidDatabase))
    
    if includeComets:
        runNameCometDatabase = "{}ELEMENTS.COMET".format(runDirectory)
        if os.path.exists(runNameCometDatabase):
            log.info("Removing the old ELEMENTS.COMET symbolic link")
            os.unlink(runNameCometDatabase)
        os.symlink("{}cometDB_version{}.dat".format(databaseFolder, 
                                                    cometDatabaseVersion),
                   runNameCometDatabase)
        log.info("Created symbolic link {}".format(runNameCometDatabase))
    mem_use(label='after creating symbolic links to the databases')
    
    #Combine the known comets and asteroids into a SOF file, which astcheck will
    #then use as input
    log.info("Combining asteroids and comets into SOF file.")
    if timeFunctions:
        t2 = time.time()
    subprocess.run("mpc2sof", cwd=runDirectory)
    if timeFunctions:
        log_timing_memory(t2, label='mpc2sof')
    
    log.info("Finished loading and formatting external databases.")
    if timeFunctions:
        log_timing_memory(t0, label='createKnownObjectsDatabase')
    
    return


# In[ ]:


def selectAsteroidsOnUncertainty(asteroidDBname, timeFunctions):
    
    """
    Go through the asteroid database (MPCORB format) and select the asteroids 
    which have orbital uncertainty parameters smaller than maxUncertainty. 
    The MPC uncertainty parameters that we consider are explained here: 
    https://www.minorplanetcenter.net/iau/info/UValue.html
    
    Overwrite the database with just the asteroids selected on their orbital
    uncertainties, so that there will be no matching with poorly known objects.
    
    Parameters:
    -----------
    asteroidDBname: string
            Name of the full asteroid database (MPCORB-formatted text-file).
    timeFunctions: boolean
            Boolean indicating whether functions need to be (wall-)timed.
    """
    mem_use(label='at start of selectAsteroidsOnUncertainty')
    
    if settingsFile.maxUncertainty == None:
        log.info("All known solar system bodies are used in the matching, "
                 +"irrespective of their uncertainty parameter.")
        return
    
    if timeFunctions:
        t0 = time.time()
    log.info("Removing asteroids with too large uncertainties...\n")
    
    #Open asteroid database
    asteroidDB = open(asteroidDBname, "r").readlines()
    
    #Find the size of the header of the asteroid database, assuming that the
    #header ends with a line of dashes.
    headerEndIndex = 0
    lineIndex = 0
    for line in asteroidDB:
        if line[:5] == "-----":
            headerEndIndex = lineIndex
            break
        lineIndex += 1
    
    #Re-write the asteroid database file, including only the header and the 
    #lines corresponding to asteroids that have small orbital uncertainties.
    Nast_preSelection = 0
    Nast_postSelection = 0
    with open(asteroidDBname, "w") as f:
        for lineIndex in range(len(asteroidDB)-1):
            
            #Copy header to file
            if lineIndex <= headerEndIndex:
                f.write(asteroidDB[lineIndex])
                continue
            
            line = asteroidDB[lineIndex]
            
            #Copy empty lines
            if line == "\n":
                f.write(line)
                continue
            
            Nast_preSelection += 1
            
            #Filter on uncertainty parameter. Copy lines of asteroids for
            #which orbits are determined reasonably well.
            uncertainty = line[105]
            if uncertainty.isdigit():
                if float(uncertainty) <= settingsFile.maxUncertainty:
                    f.write(line)
                    Nast_postSelection += 1
    
    log.info("{} out of {} asteroids have U<={}"
             .format(Nast_postSelection, Nast_preSelection,
                     settingsFile.maxUncertainty))
    log.info("Asteroid database now only includes sources with U<={}"
            .format(settingsFile.maxUncertainty))
    if timeFunctions:
        log_timing_memory(t0, label='selectAsteroidsOnUncertainty')
    
    return


# In[ ]:


def convertCatalogue2MPCformat(transientCatalogue, MPCfileName, overwriteFiles,
                               timeFunctions):
    """
    Function converts the transient catalogue to a text file of the MPC 
    80-column format, so that astcheck can run on it. For the asteroid / comet
    identifier used in the MPC file, the transient number is used. This 
    transient number cannot be used for MPC submissions as it is not all-time 
    unique (per telescope). But it is a straight-forward way to link detections 
    to known solar system objects within match2SSO.
    
    Parameters
    ----------
    transientCatalogue: string
            Path to and name of the transient catalogue.
    MPCfileName: string
            Path to and name of the MPC-formatted text file that is made in this
            function.
    overwriteFiles: boolean
            Boolean indicating whether files are allowed to be overwritten.
    timeFunctions: boolean
            Boolean indicating whether functions need to be (wall-)timed.
    """
    mem_use(label='at start of convertCatalogue2MPCformat')
    log.info("Converting transient catalogue to MPC-format.")
    if timeFunctions:
        t0 = time.time()
    
    #Load transient catalogue header
    with fits.open(transientCatalogue) as hdu:
        transientHeader = hdu[1].header
    
    #Get the MPC observatory code from the header
    MPC_code = transientHeader["MPC-CODE"].strip()
    if MPC_code not in list(pd.read_fwf("{}ObsCodes.html"
                                        .format(softwareFolder),widths=[4,2000],
                                        skiprows=1)['Code'])[:-1]:
        log.critical("MPC code {} is not in the MPC list of observatory codes"
                    .format(MPC_code))
        return None
    
    #Check if MPC-formatted file exists and if it should be overwritten or not
    if not overwriteFiles and os.path.exists(MPCfileName):
        log.info("MPC-formatted file already exists and will not re-made.")
        return MPC_code
    
    #Get observation date in the right format
    observationTime = Time(transientHeader['MJD-OBS'], format='mjd').datetime
    decimalDay = observationTime.day + (observationTime.hour + 
                 (observationTime.minute + (observationTime.second + 
                 (observationTime.microsecond/10.**6))/60.)/60.)/24.
    MPC_char16to32 = "{} {:0>2} {:08.5f} ".format(observationTime.year, 
                                                  observationTime.month, 
                                                  decimalDay)
    #Load transient catalogue data
    with fits.open(transientCatalogue) as hdu:
        detections = Table(hdu[1].data)
    
    #Create output file
    MPCfile = open(MPCfileName, "w")
    
    #Loop over the detections and add data to the MPC-formatted file
    for detection in range(len(detections)):
        #Use the source numbers as "temporary designations" in the MPC format. 
        #In this way, we will be able to link the known objects to the right 
        #source.
        MPC_char1to15 = ("     {:0>7}  C"
                         .format(detections[numberColumn][detection]))
        
        #Get the coordinates and magnitude of the source
        coord = SkyCoord(detections[RA_column][detection]*u.deg,                          detections[DEC_column][detection]*u.deg, frame='icrs')
        ra = coord.to_string('hmsdms', sep=' ', precision=2)[:11]
        dec = coord.to_string('hmsdms', sep=' ', precision=1)[-11:]
        mag = detections[magnitudeColumn][detection]
        mag = "{:.1f}".format(mag)
        
        MPC_char33to80 = "{} {}          {} G      {}".format(ra, dec,                                           mag.rjust(4), MPC_code)
        
        #Write the data to the MPC-formatted file
        MPC_line = "{}{}{}".format(MPC_char1to15, MPC_char16to32,                                    MPC_char33to80)
        MPCfile.write("{}\n".format(MPC_line))
    
    MPCfile.close()
    
    log.info("MPC-formatted file saved to {}.".format(MPCfileName))
    if timeFunctions:
        log_timing_memory(t0, label='convertCatalogue2MPCformat')
    
    return MPC_code


# In[ ]:


def runAstcheck(MPC_formattedFile, runDirectory, outputFileName, timeFunctions,
                overwriteFiles, matchingRadius=settingsFile.matchingRadius, 
                limitingMagnitude=settingsFile.limitingMagnitude, 
                maximalNumberOfObjects=settingsFile.maximalNumberOfAsteroids):
    """
    Run astcheck on the input transient catalogue to find matches between 
    transient detections and known solar system objects. Per detection, all
    matches within the matchingRadius are selected and saved to the output text
    file.
    
    Parameters:
    -----------
    MPC_formattedFile: string
            Name of the input text-file that is formatted according to the MPC's
            80-column MPC format. This file can list detections / tracklets for
            astcheck to match, but it can also contain just a single row 
            representing the observation. In the latter use case, the specified
            coordinates should correspond to the centre of the observation.
    runDirectory: string
            Directory in which astcheck is run. This directory should contain
            the mpc2sof catalogue that contains the known solar system objects.
    outputFileName: string
            Path to and name of the output text file in which the matches found
            by astcheck are stored.
    timeFunctions: boolean
            Boolean indicating whether functions need to be (wall-)timed.
    overwriteFiles: boolean
            Boolean indicating whether files are allowed to be overwritten.
    matchingRadius: int or float
            Matching radius in arcsec. The default value is the one specified in
            the settings file [set_match2sso.py].
    limitingMagnitude: int or float
            Limiting V-magnitude up to which asteroids are taken into account. 
            The default value is the one specified in the settings file 
            [set_match2sso.py].
    maximalNumberOfObjects: int
            Maximal number of matches that are returned by astcheck. The default
            value is the one specified in the settings file [set_match2sso.py].
    """
    mem_use(label='at start of runAstcheck')
    log.info("Running astcheck: matching detections to known solar system " +               "bodies.")
    if timeFunctions:
        t0 = time.time()
    
    if not overwriteFiles and os.path.exists(outputFileName):
        log.info("Astcheck output file already exists and will not be re-made.")
        return
    
    #Create a file for storing the output of the astcheck run
    outputFile = open(outputFileName, 'w')
    
    #Run astcheck from folder containing .sof-file
    subprocess.call(["astcheck", MPC_formattedFile, 
                     "-r{}".format(matchingRadius), 
                     "-m{}".format(limitingMagnitude),  
                     "-M{}".format(maximalNumberOfObjects)], 
                    stdout=outputFile, cwd=runDirectory)
    outputFile.close()
    
    log.info("Matches saved to {}.".format(outputFileName))
    if timeFunctions:
        log_timing_memory(t0, label='runAstcheck')
    
    return 


# In[ ]:


def create_SSOcatalogue(astcheckOutputFileName, runDirectory,
                        SSOcatalogueName, includeComets, timeFunctions,
                        overwriteFiles):
    """
    Open the text-file that was produced when running astcheck 
    (astcheckOutputFileName) and save the information to an SSO catalogue.
    
    Parameters:
    -----------
    astcheckOutputFileName: string  
            Name of the file containing astcheck's output (the matches). This
            can be None, in which case we will create a dummy SSO catalogue
            without matches.
    runDirectory: string
            Directory in which astcheck was run. This directory also contains
            symbolic links to the asteroid and comet databases that were used to
            create the known objects catalogue that astcheck used.
    SSOcatalogueName: string
            Name of the SSO catalogue to be created in which the matches
            are stored.
    includeComets: boolean
            Boolean indicating whether comets are included in the known objects 
            catalogue.
    timeFunctions: boolean
            Boolean indicating whether functions need to be (wall-)timed.
    overwriteFiles: boolean
            Boolean indicating whether files are allowed to be overwritten.
    """
    mem_use(label='at start of create_SSOcatalogue')
    log.info("Converting astcheck output into an SSO catalogue.")
    if timeFunctions:
        t0 = time.time()
    
    if not overwriteFiles and os.path.exists(SSOcatalogueName):
        log.info("SSO catalogue already exists and will not be re-made.")
        return
    
    #If the transient catalogue was red-flagged, matching was not performed and
    #an empty SSO catalogue needs to be created.
    if astcheckOutputFileName == None:
        log.info("Creating a dummy SSO catalogue.")
        SSOheader = create_SSOheader(runDirectory, includeComets, True, 
                                     timeFunctions)
        write2fitsFile(Table(), None, [], SSOcatalogueName, 
                       startHeader=SSOheader)
        if timeFunctions:
            log_timing_memory(t0, label='create_SSOcatalogue')
        return
    
    #Read astcheck's output file
    #The footer is variable in terms of the number of lines it spans, but it 
    #always starts with the footerString as defined below and can hence be
    #recognized by this string.
    footerString = "The apparent motion and arc length"
    astcheckOutputFile = open(astcheckOutputFileName, "r").readlines()
    indexFooter = [index for index in range(len(astcheckOutputFile)) if                    footerString in astcheckOutputFile[index]][0]
    astcheckOutputFile = astcheckOutputFile[5:indexFooter]
    
    separator = "\n"
    indicesSeparator = np.where(np.array(astcheckOutputFile)==separator)[0]
    indicesSeparator = np.append(-1, indicesSeparator)
    indicesSeparator = np.append(indicesSeparator, len(astcheckOutputFile))
    
    #Create table to store match information in
    #Dictionary of all output columns with their (numpy) formats and their units
    outputColumns = {
                    numberColumn:   ["i4",  ""],
                    "ID_SSO":       ["12a", ""],
                    "DIST_RA_SSO":  ["i2",  "arcsec"],
                    "DIST_DEC_SSO": ["i2",  "arcsec"],
                    "DIST_SSO":     ["i2",  "arcsec"],
                    "MAG_V_SSO":    ["f4",  ""]
                    }
    outputTable = Table()
    for key in outputColumns.keys():
        outputTable.add_column(Column(name=key, dtype=outputColumns[key][0],
                                     unit=outputColumns[key][1]))
    #Loop over sources
    for index in range(len(indicesSeparator)-1):
        minimalIndex = indicesSeparator[index]+1
        maximalIndex = indicesSeparator[index+1]
        
        if minimalIndex == maximalIndex:
            continue
        
        #Name of the source in the MPC-formatted input file (= transient number)
        transientNumber = astcheckOutputFile[minimalIndex:maximalIndex][0].split(                                                                 ':')[0].split()
        #Lines corresponding to matches in the astcheck output file
        matches = astcheckOutputFile[minimalIndex:maximalIndex][1:]
        
        if len(matches) == 0:
            continue
        
        #Get properties of closest match
        matchProperties = re.split('  +', matches[0])
        if len(matchProperties) == 7:
            identifier, offsetRA, offsetDEC, offset, magnitude,             properMotionRA, properMotionDEC = matchProperties
        elif len(matchProperties) == 8:
            _, identifier, offsetRA, offsetDEC, offset, magnitude,             properMotionRA, properMotionDEC = matchProperties
        else:
            log.critical("Match could not be split into correct parameters:\n{}"
                        .format(matches[0]))
            continue
        
        try:
            magnitude = float(magnitude)
        except:
            log.warning("Magnitude '{}' could not be converted to float."
                        .format(magnitude))
            magnitude = None
            
        #Add match to output table
        outputTableRow = [transientNumber, str(identifier), float(offsetRA), 
                          float(offsetDEC), float(offset), magnitude]
        outputTable.add_row(outputTableRow)
    
    #If a solar system object was matched to multiple transient sources, remove
    #all these matches as they are unreliable.
    uniqueObjects = np.unique(outputTable['ID_SSO'])
    if len(uniqueObjects) != len(outputTable):
        for obj in uniqueObjects:
            objIndices = np.where(outputTable['ID_SSO'] == obj)[0]
            if len(objIndices) > 1:
                log.warning("{} was matched to multiple transients. Removing " 
                            .format(obj) + "these matches from the SSO " +
                            "catalogue due to unreliability!")
                outputTable.remove_rows(objIndices)
    
    #Set dummy parameter (dummy means that there were no matches found)
    dummy = False
    if len(outputTable) == 0:
        dummy = True
    
    #Create header for SSO catalogue
    SSOheader = create_SSOheader(runDirectory, includeComets, dummy,
                                timeFunctions)
    
    write2fitsFile(outputTable, None, [], SSOcatalogueName, 
                   startHeader=SSOheader)
    
    log.info("Matches saved to SSO catalogue: {}".format(SSOcatalogueName))
    if timeFunctions:
        log_timing_memory(t0, label='create_SSOcatalogue')
    return


# In[ ]:


def create_SSOheader(runDirectory, includeComets, dummy, timeFunctions):
    
    """
    Function creates the header for the SSO catalogue.
    Parameters:
    -----------
    runDirectory: string
            name of the folder in which the symbolic links to the databases are
            are stored. These are used to get the version numbers of the 
            databases.
    includeComets: boolean
            Boolean indicating whether comets were included in the known objects
            database.
    dummy: boolean
            Boolean indicating whether the catalogue is a dummy catalogue without
            sources (dummy=True). If False, there are sources in the catalogue.
    timeFunctions: boolean
            Boolean indicating whether functions need to be (wall-)timed.
    """
    mem_use(label='at start of create_SSOheader')
    log.info("Creating SSO header.")
    if timeFunctions:
        t0 = time.time()
    
    #Create empty SSO header
    header = fits.Header()
    
    #Add Python version to SSO header
    header['PYTHON-V'] = (platform.python_version(), "Python version used")
    
    #Get C++ version and add to the SSO header. Based on [https://stackoverflow
    # .com/questions/44734397/which-c-standard-is-the-default-when-compiling-
    # with-g/44735016#44735016]
    proc = subprocess.run("g++ -dM -E -x c++  /dev/null | grep -F __cplusplus", 
                          capture_output=True, shell=True)
    CPP_macro = proc.stdout.decode("utf-8").replace("\n","").split()[-1]
    
    if CPP_macro not in settingsFile.CPPmacro2version.keys():
        log.error("C++ macro unknown: {}".format(CPP_macro))
        CPP_version = "None"
    else:
        CPP_version = settingsFile.CPPmacro2version[CPP_macro]
    header['CPP-V'] = (CPP_version, "C++ version used")
    
    #Get G++ version and add to SSO header
    proc = subprocess.run("g++ --version", capture_output=True, shell=True)
    GPP_version = proc.stdout.decode("utf-8").split("\n")[0].split()[-1]
    header['GPP-V'] = (GPP_version, "G++ version used")
    
    #Add match2SSO & header keyword versions to the SSO header
    header['SSO-V'] = (__version__, "match2SSO version used")
    header['SSOKW-V'] = (keywords_version, 
                         "SSO header keywords version used")
    
    #Get unique strings with git, signifying the latest commit that was made to
    #the lunar & jpl_eph repositories and hence signifynig the versions of these
    #repositories. Save the strings to the SSO header.
    proc = subprocess.run("git rev-parse --short=4 HEAD", capture_output=True, 
                          shell=True, cwd="{}lunar/".format(softwareFolder))
    lunarVersion = proc.stdout.decode("utf-8").replace("\n","")
    header['LUNAR-V'] = (lunarVersion, "lunar repository version used")
    
    proc = subprocess.run("git rev-parse --short=4 HEAD", capture_output=True, 
                          shell=True, cwd="{}jpl_eph/".format(softwareFolder))
    jplEphVersion = proc.stdout.decode("utf-8").replace("\n","")
    header['JPLEPH-V'] = (jplEphVersion, "jpl_eph repository version used")
    
    #Add version of JPL lunar & planetary ephemerides file to SSO header
    header['JPLDE-V'] = ("DE{}".format(settingsFile.JPL_ephemerisFile.split("."
                                       )[-1]),"JPL ephemeris file version used")
    
    #Add asteroid database version & reference epoch to the SSO header.
    #The MPCORB.DAT symbolic link in the run directory refers to the 
    #asteroid database version that was used. The name structure of this 
    #database is: asteroidDB_version[yyyymmddThhmm]_epoch[yyyymmddThhmm].dat
    if os.path.exists("{}MPCORB.DAT".format(runDirectory)):
        asteroidDBname = os.readlink("{}MPCORB.DAT".format(runDirectory))
        asteroidDBdate = os.path.basename(asteroidDBname).split("_"
                                                      )[1].replace("version","")
        asteroidDBversion = "{}-{}-{}T{}:{}".format(asteroidDBdate[0:4],
                                       asteroidDBdate[4:6], asteroidDBdate[6:8],
                                    asteroidDBdate[9:11], asteroidDBdate[11:13])
        
        referenceEpoch = os.path.basename(asteroidDBname).split("_"
                                                        )[2].replace("epoch","")
        asteroidDBepoch = ("{}-{}-{}T{}:{}" .format(referenceEpoch[0:4],
                                       referenceEpoch[4:6], referenceEpoch[6:8],
                                   referenceEpoch[9:11], referenceEpoch[11:13]))
    else:
        asteroidDBversion = "None"
        asteroidDBepoch = "None"
    header['ASTDB-V'] = (asteroidDBversion, 
                         "asteroid database version (date in UTC)")
    header['ASTDB-EP'] = (asteroidDBepoch, "asteroid database epoch in UTC")
    
    #Add comet database version to the SSO header.
    #The ELEMENTS.COMET symbolic link in the run directory refers to the 
    #comet database version that was used. The name structure of this 
    #database is: cometDB_version[yyyymmddThhmm].dat
    cometDBversion = "None"
    if includeComets:
        try:
            cometDBname = os.readlink("{}ELEMENTS.COMET".format(runDirectory))
            cometDBdate = os.path.basename(cometDBname).split("_"
                                                      )[1].replace("version","")
            cometDBversion = "{}-{}-{}T{}:{}".format(cometDBdate[0:4], 
                                          cometDBdate[4:6], cometDBdate[6:8], 
                                          cometDBdate[9:11], cometDBdate[11:13])
        except:
            #Readlink won't work for an empty comet database that was created 
            #when includeComet=False, because ELEMENTS.COMET is a file instead 
            #of a symbolic link.
            log.info("Empty comet database was used.")
    header['COMDB-V'] = (cometDBversion, "comet database version (date in UTC)")
    
    #Add matching radius and maximum orbital uncertainty parameter to header
    header['RADIUS'] = (float(settingsFile.matchingRadius), 
                        "matching radius in arcsec")
    header['U-MAX'] = (settingsFile.maxUncertainty, 
                       "maximum orbital uncertainty parameter")
    
    #Add keyword indicating whether there is
    header['SDUMCAT'] = (bool(dummy), "dummy SSO catalogue without sources?")
    
    log.info("SSO header complete.")
    if timeFunctions:
        log_timing_memory(t0, label='create_SSOheader')
    
    return header


# In[ ]:


def createSubmissionFile(SSOcatalogueName, MPCformattedFileName, MPC_code,
                        timeFunctions, overwriteFiles):
    """
    Make an MPC submission file using the SSO catalogue and the MPC-formatted 
    file that were created within match2SSO to link the transient detections 
    from a single catalogue to known solar system objects. The detections 
    corresponding to matches are grouped in a 'known objects submission file'. 
    The identifiers used in the submission file are the packed designation of 
    the matching objects. These are the packed permanent designations if 
    available. Otherwise, the packed provisional designations are used.
    
    Parameters:
    -----------
    SSOcatalogueName: string
            Path to and name of the SSO catalogue of which the matches need to 
            be converted to a submission file.
    MPCformattedFileName: string
            Name of the MPC formatted file that was made in match2SSO for the
            matching (but does not contain the correct SSO identifiers yet for
            submission to the MPC).
    MPC_code: string
            MPC code corresponding to the telescope.
    timeFunctions: boolean
            Boolean indicating whether functions need to be (wall-)timed.
    overwriteFiles: boolean
            Boolean indicating whether files are allowed to be overwritten.
    """
    mem_use(label='at start of createSubmissionFile')
    log.info("Creating MPC submission file.")
    if timeFunctions:
        t0 = time.time()
    
    #Compose submission file name
    submissionFileVersion = Time.now().strftime("%Y%m%dT%H%M%S")
    submissionFileName = "{}{}_{}.txt".format(submissionFolder, 
                          os.path.basename(SSOcatalogueName).replace(".fits",
                                                                     "_submit"),
                          submissionFileVersion)
    
    #Check if file already exists (will only happen when running this function
    #multiple times in close succession, as the production time is used in the
    #file name)
    if not overwriteFiles and os.path.exists(submissionFileName):
        log.info("Submission file already exists and will not be re-made.")
        return
    
    #Open SSO catalogue
    with fits.open(SSOcatalogueName) as hdu:
        SSOcatalogue = Table(hdu[1].data)
    
    #Check SSO catalogue for matches
    if len(SSOcatalogue)==0:
        log.info("No matches found. Submission file will not be made.")
        return
    
    #Create submission file (in same folder as SSO catalogue)
    log.info("Making submission file {}".format(submissionFileName))
    if os.path.exists(submissionFileName):
        log.warning("MPC submission file {} is overwritten."
                    .format(submissionFileName))
    
    submissionFile = open(submissionFileName, 'w')
    
    #Create header for the submission file
    submissionHeader = createSubmissionHeader(submissionFileName, MPC_code,
                                             timeFunctions)
    submissionFile.write(submissionHeader)
    
    #Open MPC-formatted file
    detections_MPCformat = pd.read_fwf(MPCformattedFileName, widths=[14,66],
                                 names=["char1to14", "char15to80"],
                                 dtype={'char1to14':np.int32, 'char15to80':str})
    
    #For each detection that was matched to a known solar system object, 
    #get the packed designation of the matching object and write the detection
    #to the submission file.
    for matchIndex, transientNumber in enumerate(SSOcatalogue[numberColumn]):
        designation = SSOcatalogue["ID_SSO"][matchIndex].strip()
        #Get packed designation
        if re.match("^[0-9]{4}\s[A-Z]", designation) or "/" in designation:
            #Asteroid or comet with provisional designation (or survey 
            #designation)
            packedDesignation = packProvisionalDesignation(designation)
            if packedDesignation == None:
                continue
            char1to12 = "    {}".format(packedDesignation)
            
        else:
            #Asteroid or comet with permanent designation
            packedDesignation, fragment = packPermanentDesignation(designation)
            if packedDesignation == None:
                continue
            char1to12 = "{}{}".format(packedDesignation, fragment.rjust(7))
            
        #Get detection details
        detectionIndex = np.where(np.array(detections_MPCformat["char1to14"]) ==
                                  int(transientNumber))[0]
        if len(detectionIndex) != 1:
            log.error("{} detections found that correspond to transient number "
                      .format(len(detectionIndex)) + "{}. Should be only one."
                      .format(transientNumber))
            continue
        detectionDetails = detections_MPCformat["char15to80"][detectionIndex[0]]
        
        #Write detection to submission file
        detectionLine = "{}  {}".format(char1to12, detectionDetails)
        if len(detectionLine) != 80:
            log.error("Detection not formatted correctly in 80 columns:\n{}"
                     .format(detectionLine))
        submissionFile.write(detectionLine+"\n")
    
    submissionFile.close()
    log.info("MPC submission file saved to {}".format(submissionFileName))
    if timeFunctions:
        log_timing_memory(t0, label='createSubmissionFile')
    
    return


# In[ ]:


def createSubmissionHeader(submissionFileName, MPC_code, timeFunctions, 
                           comment=None):
    """
    Function composes the header of the MPC submission file corresponding to a
    single transient catalogue.
    
    Parameters:
    -----------
    submissionFileName: string
            Name of the submission file for which the header is composed.
    MPC_code: string
            MPC code of the telescope with which the observation was made.
    timeFunctions: boolean
            Boolean indicating whether functions need to be (wall-)timed.
    comment: string
            Comment to be added to the header in the COM line. By default, this
            is None, meaning that the COM line is not added to the header.
    """
    mem_use(label='at start of createSubmissionHeader')
    log.info("Creating header for submission file.")
    if timeFunctions:
        t0 = time.time()
    
    firstline = "COD {}\n".format(MPC_code)
    defaultHeader = (
       "CON Radboud University, Houtlaan 4, 6525XZ, Nijmegen, The Netherlands\n"
       +"CON [p.groot@astro.ru.nl]\n"
       +"OBS P. J. Groot, S. L. D. Bloemen, L. Townsend\n"
       +"MEA P. M. Vreeswijk, D. L. A. Pieterse, K. Paterson\n"
       +"TEL 0.65-m reflector + CCD\n"
       +"NET Gaia-DR2\n"
       +"AC2 mpc-response@blackgem.org\n")
    

    #Special cases for which a phrase needs to be included in the ACK line 
    #of the header of the MPC submission file:
    #neocand = "NEO CANDIDATE" #submitting new NEO candidate
    #neocp = "NEOCP"           #submitting observations of NEOCP objects
    #Add ACK line to the header of the MPC submission file
    
    ACK_line = "ACK {}\n".format(Path(submissionFileName).stem)
    if len(ACK_line)>82:
        log.error("ACK line in submission file {} is too long!"
                  .format(submissionFileName))
    
    #Add COM line to the header
    COM_line = ''
    if comment != None:
        if len(comment)>76:
            log.warning("COM line is too long and therefore not used. " 
                        +"Use at most 76 characters!")
        else:
            COM_line = "COM {}\n".format(comment)
    
    log.info("Submission file header complete.")
    if timeFunctions:
        log_timing_memory(t0, label='createSubmissionHeader')
    
    return firstline + defaultHeader + ACK_line + COM_line


# ### Helper functions to make an MPC submission file
# Convert asteroid designations to their packed form

# In[ ]:


def abbreviateNumber(num):
    
    """
    Number packing function needed to pack MPC designations.
    """
    
    Ndict = {str(index): letter for index, letter in              enumerate(ascii_uppercase + ascii_lowercase, start=10)}
    
    if int(num)>9:
        return Ndict[str(num)]
    else:
        return num


# In[ ]:


def packProvisionalDesignation(fullDesignation):
    
    """
    Function converts provisional minor planet designation into its packed form,
    using the definitions given in
    https://www.minorplanetcenter.net/iau/info/PackedDes.html#prov .
    As described in https://www.minorplanetcenter.net/iau/info/OpticalObs.html,
    for comets a character is added in front of the provisional designation (at
    column 5), describing the comet type. For asteroids, we add a space in front
    so that the returned string is 8 characters long (spanning columns 5-12 in
    the submission file). 
    
    Parameters
    ----------
    fullDesignation: string
            Unpacked provisional designation assigned to the object by the MPC.
    """
    
    #Remove space before or after designation
    fullDesignation = fullDesignation.strip()
    
    #Their are four special survey designation forms (for surveys that were 
    #undertaken between 1960 and 1977)that should be packed differently
    surveyStrings = ["P-L", "T-1", "T-2", "T-3"]
    for surveyString in surveyStrings:
        if re.match("^[0-9]{4}\s"+"{}$".format(surveyString),fullDesignation):
            packedDesignation = "{}S{}".format(surveyString.replace("-",""), 
                                              fullDesignation[:4])
            return packedDesignation
    
    packYear = {18: "I", 19: "J", 20: "K"}
    
    def packCycleNumber(Ncycle):
        
        """Input parameter Ncycle is a string of 0-3 digits."""
        
        if len(Ncycle) == 0:
            return "00"
        
        elif int(Ncycle)>99:
            return "{}{}".format(abbreviateNumber(Ncycle[0:2]), Ncycle[2])
        
        return "{:0>2}".format(Ncycle)
    
    
    #For a comet
    if "/" in fullDesignation:
        cometType, designation = fullDesignation.split("/")
        
        #In case of a comet fragment, the last character of the packed 
        #designation is the fragment letter. Otherwise, it is zero.
        fragment = "0"
        if '-' in designation:
            designation, fragment = designation.split("-")
            fragment = fragment.lower()
        year, remainder = designation.split(" ")
        
        if int(year[:2]) not in packYear.keys():
            log.error("Provisional designation of comet {} "
                     .format(fullDesignation) +
                     "cannot be packed. Skipping it.")
            return None
        
        packedYear = "{}{}".format(packYear[int(year[:2])], year[2:])
        
        
        #In case there are two letters after the space in the designation. This
        #can be the case if the object was thought to be an asteroid early on.
        if remainder[1].isalpha():
            
            if fragment!="0":
                #A comet with two letters in its provisional designation after
                #the space and a fragment letter cannot be submitted in the old
                #submission format. It can in the new ADES format, but we are
                #not yet using this. Skip detection.
                log.error("Provisional designation of comet {} "
                          .format(fullDesignation) +
                          "cannot be packed. Skipping it.")
                return None
            
            #Although this object is not a fragment, its provisional designation
            #does contain a second letter after the space which should be 
            #written to the same position as the fragment letter.
            fragment = remainder[1]
            remainder = "{}{}".format(remainder[0], remainder[2:])
            
        #There should be at most three digits after the space-letter combination
        #in the provisional designation
        if len(remainder)>4:
            log.error("Unclear how to pack provisional designation of comet "
                      + "{}. Skipping it.".format(fullDesignation))
            return None
        
        if int(year[:2]) not in packYear.keys():
            log.error("Data from before 1800 or after 2099 cannot be assigned a"
                      + " packed provisional designation.")
            return None
        
        packedDesignation = ("{}{}{}{}{}".format(cometType, packedYear,
                        remainder[0], packCycleNumber(remainder[1:]), fragment))
        return packedDesignation
    
    #For an asteroid
    if int(fullDesignation[:2]) not in packYear.keys():
        log.error("Provisional designation of asteroid {} "
                 .format(fullDesignation) +
                 "cannot be packed. Skipping it.")
        return None
    
    packedYear = "{}{}".format(packYear[int(fullDesignation[:2])], 
                               fullDesignation[2:4])
    packedDesignation = (" {}{}{}{}".format(packedYear, fullDesignation[5],
                                           packCycleNumber(fullDesignation[7:]),
                                           fullDesignation[6]))
    
    #Final check
    if len(packedDesignation)!=8:
        log.error("Packed provisional designation is of incorrect length: '{}'"
                 .format(packedDesignation))
        return None
    
    return packedDesignation


# In[ ]:


def packPermanentDesignation(fullDesignation):
    
    """
    Function converts the permanent minor planet designation into its packed 
    form (5 characters), using the definitions given in
    https://www.minorplanetcenter.net/iau/info/PackedDes.html#perm
    Return the packed designation and - if applicable - the letter 
    corresponding to the comet fragment. If the object is not a comet fragment,
    an empty string will be returned for the fragment letter.
    Parameters
    ----------
    fullDesignation: string
            Unpacked permanent designation assigned to the object by the MPC.
    """
    fragment = ''
    
    if not fullDesignation.isdigit():
        #Object is a comet
        if len(fullDesignation.split("-"))==2:
            designation, fragment = fullDesignation.split("-")
        else:
            designation = fullDesignation
            fragment = ''
        packedDesignation = designation.zfill(5)
        fragment = fragment.lower()
        
    elif int(fullDesignation)<99999:
        packedDesignation = "{:0>5}".format(int(fullDesignation))
    
    elif int(fullDesignation)<620000:
        quotient = int(fullDesignation)//10000
        packedDesignation = "{}{:0>4}".format(abbreviateNumber(quotient), 
                                              int(fullDesignation)%10000)
    
    else:
        remainder = int(fullDesignation) - 620000
        quotient3 = remainder//62**3
        remainder -= quotient3*62**3
        quotient2 = remainder//62**2
        remainder -= quotient2*62**2
        quotient1 = remainder//62
        remainder -= quotient1*62
        packedDesignation = "~{}{}{}{}".format(abbreviateNumber(quotient3), 
                                               abbreviateNumber(quotient2), 
                                               abbreviateNumber(quotient1), 
                                               abbreviateNumber(remainder))
    #Final check
    if len(packedDesignation)!=5:
        log.error("Packed permanent designation is of incorrect length: '{}'"
                  .format(packedDesignation))
        return None, ''
    
    return packedDesignation, fragment


# ### General helper functions

# In[ ]:


#BlackBOX function
# from https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse
def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')
        return


# In[ ]:


def checkInputParameters(mode, cat2process, date2process, list2process):
    
    """
    Check if the correct (combination of) input parameters was/were defined for
    run_match2SSO. If so, this function returns True. Otherwise, False is 
    returned.
    """
    if mode not in ["day", "night", "historic"]:
        print("CRITICAL: unknown mode.")
        return False
    
    if mode=="day":
        if cat2process!=None:
            print("CRITICAL: when processing a specified catalog, the mode must"
                  " be 'historic' or 'night'.")
            return False
        if list2process!=None:
            print("CRITICAL: when processing a catalog list, the mode must be "
                  "'historic'.")
            return False
    
    elif mode=="night":
        if cat2process==None:
            print("Critical: --catalog needs to be specified when running "
                  "match2SSO in night mode.")
            return False
        if date2process!=None:
            print("CRITICAL: when processing a specified date, the mode must "
                  "be 'historic' or 'day'.")
            return False
        if list2process!=None:
            print("CRITICAL: when processing a catalog list, the mode must be "
                  "'historic'.")
            return False
    
    elif mode=="historic" and cat2process==None and (date2process==None and
                                                     list2process==None):
        print("CRITICAL: --date, --catalog and --catlist are all None. Nothing "
              "to process.")
        return False
        
    if (date2process!=None and cat2process!=None) or (date2process!=None and
        list2process!=None) or (cat2process!=None and list2process!=None):
        print("CRITICAL: either specify --date, --catalog OR --catlist. A "
              "combination is not allowed.")
        return False
    
    if cat2process!=None:
        if "trans" not in cat2process:
            print("CRITICAL: specified catalog should correspond to a transient"
                  " catalog.")
            return False
        if not os.path.exists(cat2process):
            print("CRITICAL: the specified catalog does not exist.")
            return False
    
    if list2process!=None and not os.path.exists(list2process):
        print("CRITICAL: the specified catalog list does not exist.")
        return False
    
    return True


# In[ ]:


def checkSettingsFile(tel):
    
    #Set up global variables for all folders in the settings file and ensure the
    #paths end with a slash
    global inputDataFolder
    inputDataFolder = get_par(settingsFile.inputDataFolder, tel)
    if inputDataFolder[-1] != "/":
        inputDataFolder += "/"
    
    global softwareFolder
    softwareFolder = get_par(settingsFile.softwareFolder, tel)
    if softwareFolder[-1] != "/":
        softwareFolder += "/"
    
    global databaseFolder
    databaseFolder = get_par(settingsFile.databaseFolder, tel)
    if databaseFolder[-1] != "/":
        databaseFolder += "/"
    
    global logFolder
    logFolder = get_par(settingsFile.logFolder, tel)
    if logFolder[-1] != "/":
        logFolder += "/"
    
    global submissionFolder
    submissionFolder = get_par(settingsFile.submissionFolder, tel)
    if submissionFolder[-1] != "/":
        submissionFolder += "/"
    
    #Check if database folder and submission folder exists (it will be needed / 
    #prepared in any mode) and create the folders if they don't exist
    if not os.path.isdir(databaseFolder):
        os.makedirs(databaseFolder)
    if not os.path.exists(submissionFolder):
        os.makedirs(submissionFolder)
    
    #Check if software folder exists. If it doesn't, exit match2SSO
    if not os.path.isdir(softwareFolder):
        print("CRITICAL: software folder given in settings file doesn't exist.")
        return False
    
    #Check that astcheck parameters are numbers
    if not (isinstance(settingsFile.matchingRadius, float) or 
            isinstance(settingsFile.matchingRadius, int)):
        print("CRITICAL: incorrectly specified matching radius in settings "
              +"file. Must be float or integer.")
        return False
    
    if not (isinstance(settingsFile.limitingMagnitude, float) or
            isinstance(settingsFile.limitingMagnitude, int)):
        print("CRITICAL: incorrectly specified limiting mag. in settings file.")
        return False
    
    if not isinstance(settingsFile.maximalNumberOfAsteroids, int):
        print("CRITICAL: incorrectly specified max. number of asteroids in "
              +"settings file.")
        return False
    
    if (not isinstance(settingsFile.maxUncertainty, int) and 
        settingsFile.maxUncertainty != None):
        print("CRITICAL: incorrectly specified max. uncertainty in settings "
              +"file. Must be 0-9 or None.")
        return False
        
    #Check if JPL ephemeris file exists
    if not os.path.exists(settingsFile.JPL_ephemerisFile):
        print("CRITICAL: JPL ephemeris file specified in settings file doesn't "
              +"exist.")
        return False
    
    return True


# In[ ]:


def setUpLogFile(logName):
    
    """
    Function creates log file and configures the log handler.
    """
    
    if logName==None:
        return
    
    logDir, logFileName = os.path.split(logName)
    
    #If no folder is specified, use the log folder from the settings file.
    if len(logDir)==0:
        logDir = logFolder
    
    #Create folder to store log in, if it does not yet exist
    if not os.path.isdir(logDir):
        os.makedirs(logDir)
    
    #Configure log handling
    logFile = "{}/{}" .format(logDir, logFileName)
    if os.path.exists(logFile):
        filePathAndName, extension = os.path.splitext(logFile)
        logFile = "{}_{}{}".format(filePathAndName, 
                                Time.now().strftime("%Y%m%d_%H%M%S"), extension)
        print("Log file already exists. Creating a new log named {}"
              .format(logFile))
        
    fileHandler = logging.FileHandler(logFile, 'a')
    fileHandler.setFormatter(logFormatter)
    fileHandler.setLevel('INFO')
    log.addHandler(fileHandler)
    
    return


# In[ ]:


#ZOGY function
def get_par(par, tel):
    
    """Function to check if [par] is a dictionary with one of the keys
       being [tel] or the alphabetic part of [tel] (e.g. 'BG'), and if
       so, return the corresponding value. Otherwise just return the
       parameter value."""

    par_val = par
    if type(par) is dict:
        if tel in par:
            par_val = par[tel]
        else:
            # cut off digits from [tel]
            tel_base = ''.join([char for char in tel if char.isalpha()])
            if tel_base in par:
                par_val = par[tel_base]
        
    return par_val


# In[ ]:


#ZOGY function
def mem_use(label=''):

    # ru_maxrss is in units of kilobytes on Linux; however, this seems
    # to be OS dependent as on mac os it is in units of bytes; see
    # manpages of "getrusage"
    if sys.platform=='darwin':
        norm = 1024**3
    else:
        norm = 1024**2
        
    mem_max = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/norm
    mem_now = psutil.Process().memory_info().rss / 1024**3
    mem_virt = psutil.Process().memory_info().vms / 1024**3
    
    log.info ('memory use [GB]: rss={:.3f}, maxrss={:.3f}, vms={:.3f} in {}'
              .format(mem_now, mem_max, mem_virt, label))

    return


# In[ ]:


#ZOGY function
def log_timing_memory(t_in, label=''):
  
    log.info('wall-time spent in {}: {:.3f} s'.format(label, time.time()-t_in))
    mem_use(label=label)

    return


# In[ ]:


def getTransientFileNames(minimalDate, maximalDate, tel, timeFunctions,
                         excludeFlagged=False):
    
    """
    Function returns a list with the transient file names that were taken
    between the minimal and maximal specified dates. The input dates should be
    Time objects. If excludeFlagged is True, the dummy transient catalogues
    (which are red-flagged) are excluded.
    
    Parameters
    ----------
    minimalDate: datetime object, incl time zone
            Minimal observation date of the time block for which the 
            observations are selected.
    maximalDate: datetime object, incl time zone
            Maximal observation date of the time block for which the
            observations are selected.
    tel: string
            Telescope abbreviation.
    timeFunctions: boolean
            Boolean indicating whether functions need to be (wall-)timed.
    excludeFlagged: boolean
            Boolean indicating whether red-flagged (dummy) catalogues should be
            excluded or not.
    """
    mem_use(label='at start of getTransientFileNames')
    log.info("Selecting transient catalogues between {} and {}."              .format(minimalDate.strftime("%Y-%m-%d %H:%M:%S"), 
                     maximalDate.strftime("%Y-%m-%d %H:%M:%S")))
    if timeFunctions:
        t0 = time.time()
    
    #Convert to local time
    localTimeZone = timezone(get_par(settingsFile.timeZoneTelescope, tel))
    minimalDate = minimalDate.astimezone(localTimeZone)
    maximalDate = maximalDate.astimezone(localTimeZone)
    
    #Select the transient files by observation date
    year, month, day = "*", "*", "*"
    if minimalDate.year == maximalDate.year:
        year = "%d"%(minimalDate.year)
        if minimalDate.month == maximalDate.month:
            month = "{:0>2}".format(minimalDate.month)
            if minimalDate.day == maximalDate.day:
                day = "{:0>2}".format(maximalDate.day)
    transientFiles = glob.glob(os.path.join(inputDataFolder, 
                                            "%s/%s/%s/*_trans_light.fits"
                                            %(year, month, day)))
    if len(transientFiles)==0:
        return []
    
    files2process = []
    for transientFile in transientFiles:
        #Parse date encoded in filename and compare with our limits
        #(e.g. ML1_20200517_034221_red_trans_light.fits)
        splittedFileName = os.path.basename(transientFile).split("_")
        date_obs = splittedFileName[1]
        time_obs = splittedFileName[2]
        observationTime = Time.strptime(date_obs+time_obs, "%Y%m%d%H%M%S").mjd
        if (observationTime >= Time(minimalDate).mjd and 
            observationTime <= Time(maximalDate).mjd):
            with fits.open(transientFile) as hdu:
                header = hdu[1].header
            
            if not excludeFlagged:
                files2process.append(transientFile)
            else:
                log.info("Excluding red-flagged (dummy) catalogues.")
                
                if dummyColumn not in header.keys():
                    log.critical("{} not in the header!".format(dummyColumn))
                    return []
                
                if header[dummyColumn]==False:
                    files2process.append(transientFile)
    
    log.info("{} transient catalogues have been selected."
             .format(len(files2process)))
    if timeFunctions:
        log_timing_memory(t0, label='getTransientFileNames')
    return files2process


# In[ ]:


def getNightStartFromCatalogueName(catalogueName, tel, noonType="local"):
    
    """
    This function returns the noon corresponding to the start of the
    observation night, as a datetime object. This is either the local noon or
    the noon in UTC, as specified. The noon is deduced from the information in
    the catalogue name (e.g. ML1_yyyymmdd_hhmmss_red_trans_light.fits).
    
    Parameters
    ----------
    catalogueName: string
            Name of the catalogue corresponding to an observation that took 
            place on the observation night for which the noon that signifies the
            start of the night must be determined.
    tel: string
            Telescope abbreviation.
    noonType: string
            Must be either "local" or "utc". If "utc", this function will return
            the noon corresponding to the start of the night in UTC. This can be
            different from the local noon.
    """
    noonType = noonType.lower()
    
    #Get observation time from filename and define as being in UTC
    splittedFileName = os.path.basename(catalogueName).split("_")
    observationDate = splittedFileName[1]
    observationHour = splittedFileName[2]
    observationTime = datetime.strptime(observationDate+" "+observationHour, 
                                   "%Y%m%d %H%M%S").replace(tzinfo=pytz.utc)
    
    if noonType == "local":
        localTimeZone = timezone(get_par(settingsFile.timeZoneTelescope, tel))
        
        #Get local noon corresponding to the start of the observing night
        local_noon = localTimeZone.localize(datetime.strptime(observationDate+
                                                     " 120000","%Y%m%d %H%M%S"))
        #Get date of observing night
        if observationTime < local_noon:
            date = (observationTime - timedelta(days=1)).strftime("%Y%m%d")
        else:
            date = observationTime.strftime("%Y%m%d")
        
        #Make local noon variable
        startNight = localTimeZone.localize(datetime.strptime(date+" 120000",
                                                              "%Y%m%d %H%M%S"))
    else:
        if noonType != "utc":
            log.error("Noon type not understood. Assuming noon in utc.")
        
        startNight = pytz.utc.localize(datetime.strptime(
                                    observationDate+" 120000","%Y%m%d %H%M%S"))
        if int(observationHour[:2])<12.:
            startNight = startNight - timedelta(days=1)
    
    return startNight


# In[ ]:


def checkInputCatalogue(catalogueName):
    
    """
    Check if the input catalogue exists and if it is a dummy (red-flagged)
    catalogue or not. If a light version of the catalogue is available, use that
    version. This function returns a boolean for "does the catalogue exist?", a
    boolean for "is the catalogue a dummy?" and the catalogue name is returned,
    as the light version might have been selected instead of the transient
    catalogue that includes the thumbnails.
    """
    
    #Check whether the (light) catalogue exists and ensure the use of the light
    #version of the catalogue if it is available (better in terms of memory 
    #usage & processing speed)
    if "_light" not in catalogueName:
        lightCatalogueName = catalogueName.replace("_trans", "_trans_light")
        if os.path.exists(lightCatalogueName):
            catalogueName = lightCatalogueName
            
    if not os.path.exists(catalogueName):
        log.critical("The specified catalog does not exist:\n{}"
                     .format(catalogueName))
        return False, None, catalogueName
            
    #Check quality control flag of the catalogue
    with fits.open(catalogueName) as hdu:
        header = hdu[1].header
            
    if dummyColumn not in header.keys():
        log.critical("{} not in the header of {}!"
                     .format(dummyColumn, catalogueName))
        return False, None, catalogueName

    if header[dummyColumn]:
        log.info("{} is a dummy catalogue.".format(catalogueName))
        return True, True, catalogueName
    
    return True, False, catalogueName


# In[ ]:


def checkForDatabaseProducts(runDirectory):
    
    """
    This function checks if the database products that astcheck needs in order
    to process transient catalogues are located in the runDirectory. These are
    the known objects catalogue, the symbolic link to the asteroid catalogue
    (needed for reading out the asteroid database version) and ELEMENTS.COMET,
    which is either an empty comet database (if comets should not be included in
    the matching) or a symbolic link to the comet database used (again needed to
    read out the version).
    The function also indirectly checks for the existence of the runDirectory.
    A boolean is returned: True if all is well and False if something is
    missing.
    
    Parameters
    ----------
    runDirectory: string
            Directory in which astcheck will be run.
    """
    
    #Check for known objects catalogue
    if not os.path.exists("{}mpcorb.sof".format(runDirectory)):
        log.critical("The known objects catalogue (SOF format) could not be "
                     + "found.")
        return False
    
    #Check for symbolic links pointing to the used version of the SSO databases
    if not os.path.exists("{}MPCORB.DAT".format(runDirectory)):
        log.critical("MPCORB.DAT could not be found")
        return False
    if not os.path.exists("{}ELEMENTS.COMET".format(runDirectory)):
        log.critical("ELEMENTS.COMET could not be found")
        return False
    
    return True


# In[ ]:


def write2fitsFile(data, header, headerKeys, outputFileName, startHeader=None):
    
    """
    Function formats the output data, composes the header and combines the two
    into a hdu table. The table is then written to a fits file.
    
    Parameters
    ----------
    data : table data
            Table data which is to be used as the data for the output fits
            table.
    header: header
            Header from which certain keywords are copied to the header of the
            output catalogue.
    headerKeys: list of strings
            Contains names of header keywords from the header mentioned above, 
            that will be included in the output catalogue header.
    outputFileName: string
            File name (including path) under which the output binary fits table 
            will be stored.
    startHeader: header
            Header which will be included in the header of the output catalogue.
            startHeader can be None.
    """
    mem_use(label='at start of write2fitsFile')
    
    #Format fits table data
    columns = []
    for columnName in data.columns:
        
        columnFormat = data[columnName].dtype
        
        #Converting bytestring format to fits format does not work properly for
        #strings, as the length is not taken into account properly. Manually
        #correct this.
        if 'S' in columnFormat.str:
            stringLength = columnFormat.str.split("S")[-1]
            columnFormat = "{}A".format(stringLength)
        
        columnUnit = str(data[columnName].unit)
        if columnUnit=="None":
            columnUnit = ""
        
        column = fits.Column(name=columnName, format=columnFormat,
                             unit=columnUnit, array=data[columnName])
        columns.append(column)
    
    #Compose fits table header
    if startHeader != None:
        header = startHeader
    else:
        header = fits.Header()
        
    for key in headerKeys:
        header[key] = (header[key], header.comments[key])
   
    #Combine formatted fits columns and header into output binary fits table
    fitsTable = fits.BinTableHDU.from_columns(columns, header=header)
    fitsTable.writeto(outputFileName, overwrite=True)
    
    mem_use(label='at end of write2fitsFile')
    return


# In[ ]:


#Based on BlackBOX function clean_tmp
def removeTemporaryFolder(tmp_path):
    
    """
    Function that removes the specified folder and its contents.
    """

    if os.path.isdir(tmp_path):
        shutil.rmtree(tmp_path)
        log.info ('Removing temporary folder: {}'.format(tmp_path))
    
    else:
        log.warning ('tmp folder {} does not exist'.format(tmp_path))
    
    mem_use(label='after removing temporary folder')
    return


# ### Run match2SSO using command line parameters

# In[ ]:


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="User parameters")
    
    parser.add_argument("--telescope", type=str, default="ML1", 
                        help="Telescope name (ML1, BG2, BG3 or BG4); "
                        "default='ML1'")
    
    parser.add_argument("--mode", type=str, default="historic", 
                        help="Day, night or historic mode of pipeline; "
                        "default='historic'")
    
    parser.add_argument("--catalog", type=str, default=None, 
                        help="Only process this particular transient catalog. "
                        "Requires full path and requires mode to be 'historic' "
                        "or 'night'; default=None")
    
    parser.add_argument("--date", type=str, default=None,
                        help="Date to process (yyyymmdd, yyyy-mm-dd, yyyy/mm/dd"
                        " or yyyy.mm.dd). Mode is required to be 'historic'; "
                        "default=None")
    
    parser.add_argument("--catlist", type=str, default=None, 
                        help="Process all transient catalogs in the input list."
                        " List entries require full path. Mode  must be "
                        "'historic'; default=None")
    
    parser.add_argument("--logname", type=str, default=None,
                        help="Name of log file to save. Requires full path; "
                        "default of None will not create a logfile")
    
    parser.add_argument("--keep_tmp", type=str, default="False",
                        help="Boolean to indicate if temporary directories / "
                        "files need to be kept; default=False")
    
    parser.add_argument("--newdatabases", type=str, default="True",
                        help="Boolean to indicate if the asteroid and comet "
                        "databases need to be redownloaded. If False and "
                        "downloaded versions already exist, the newest "
                        "downloaded versions are used; default=True")
    
    parser.add_argument("--includecomets", type=str, default="False",
                        help="Boolean to indicate if comets are to be "
                        "included in the matching; default=False")
    
    parser.add_argument("--overwrite", type=str, default="True",
                        help="Boolean to indicate if files are allowed to be "
                        "overwritten; default=True")
    
    parser.add_argument("--timing", type=str, default="False",
                        help="Boolean to indicate if functions need to be "
                        "(wall-)timed; default=False")
    
    args = parser.parse_args()
    run_match2SSO(tel=args.telescope, mode=args.mode, cat2process=args.catalog, 
                  date2process=args.date, list2process=args.catlist, 
                  logName=args.logname, keep_tmp=args.keep_tmp, 
                  redownloadDatabases=args.newdatabases, 
                  includeComets=args.includecomets, overwriteFiles=args.overwrite,
                  timeFunctions=args.timing)

