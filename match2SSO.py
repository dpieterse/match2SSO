#!/usr/bin/env python
# coding: utf-8

# # match2SSO
# * Running instructions are given at the start of function run_match2SSO
# * Originally written as a Jupyter Notebook using Python 3.8.10
# * Compatible with MeerLICHT / BlackGEM observations
# * Compatible with BlackBOX / ZOGY version 1.0.0 and up (image processing pipeline of MeerLICHT & BlackGEM)
# * Compatible with Unix file systems and Google cloud storage systems
# 
# Output (SSO) catalogue columns and header keywords are listed here:
# https://www.overleaf.com/read/zrhqwcbkfqns
# 
# <i>match2SSO</i> makes grateful use of the <i>lunar</i> and 
# <i>jpl_eph</i> repositories that were written by Bill Gray under 
# Project Pluto. The core of <i>match2SSO</i> is <i>astcheck</i>: a C++
# script in the <i>lunar</i> repository that matches detections to known
# solar system objects. More information on <i>astcheck</i> can be found
# at: https://www.projectpluto.com/astcheck.htm
# 
# <b>Dependencies on scripts and files</b>
# * <i>lunar</i> package (https://github.com/Bill-Gray/lunar)
# * <i>jpl_eph</i> package (https://github.com/Bill-Gray/jpl_eph)
# * JPL DE ephemeris file (ftp://ssd.jpl.nasa.gov/pub/eph/planets/Linux/)
# * MPC's Observatory codes list
#   (https://www.minorplanetcenter.net/iau/lists/ObsCodes.html)
# 
# In addition, <i>match2SSO</i> uses MPCORB.DAT (MPC's asteroid database) and
# COMET.ELEMENTS (JPL's comet database), but these are downloaded when
# running the script and hence do not need to be pre-downloaded.

# ## Python packages and settings

# In[ ]:


__version__ = "1.2.0"
__author__ = "Danielle Pieterse"
KEYWORDS_VERSION = "1.1.0"


# In[ ]:


# Standard imports
import os
import re
import sys
import glob
import time
import argparse
import logging
import resource
import platform
import shutil
import subprocess
from math import sqrt
from datetime import datetime, timedelta
from pathlib import Path
from string import ascii_lowercase, ascii_uppercase
import psutil
import requests
import numpy as np

# Third party imports
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table, Column
from astropy.time import Time
import pandas as pd
import pytz
from pytz import timezone

# In case a Google cloud storage system is being used
from google.cloud import storage

# Local imports
import set_match2SSO as settingsFile # Load match2SSO settings file


# In[ ]:


# Set log format and global logging constants
LOG_FORMAT = ("%(asctime)s.%(msecs)03d [%(levelname)s, %(process)s] "
              "%(message)s [%(funcName)s, line %(lineno)d]")
DATE_FORMAT = "%Y-%m-%dT%H:%M:%S"
logging.basicConfig(level="INFO", format=LOG_FORMAT, datefmt=DATE_FORMAT)
LOG_FORMATTER = logging.Formatter(LOG_FORMAT, DATE_FORMAT)
logging.Formatter.converter = time.gmtime # Convert time in logger to UTC
LOG = logging.getLogger()


# In[ ]:


# Set global constants

# Relevant data columns
NUMBER_COLUMN = settingsFile.colNumber
RA_COLUMN = settingsFile.colRA
DEC_COLUMN = settingsFile.colDec
MAG_COLUMN = settingsFile.colMag
SNR_COLUMN = settingsFile.colSNR

# Relevant header keywords of detection catalogue
DUMMY_KEYWORD = settingsFile.keyDummy
DATE_KEYWORD = settingsFile.keyDate
MPC_CODE_KEYWORD = settingsFile.keyMPCcode
CENTRAL_RA_KEYWORD = settingsFile.keyRACentre
CENTRAL_DEC_KEYWORD = settingsFile.keyDecCentre
LIMMAG_KEYWORD = settingsFile.keyLimmag

# Load switches from settings file
INCLUDE_COMETS = bool(settingsFile.include_comets)
KEEP_TMP = bool(settingsFile.keep_tmp)
TIME_FUNCTIONS = bool(settingsFile.time_functions)


# ## Main routines to call match2SSO functions in the correct order

# In[ ]:


def run_match2SSO(tel, mode, cat2process, date2process, list2process,
                  logname, redownload, overwrite):
    """
    Run match2SSO on the input catalogue(s)/date. match2SSO can be run in
    different mode / date2process / cat2process / list2process combinations.
    Allowed combinations are: (if not mentioned, the variable is None)
      * Day mode
      * Day mode + date2process
      * Night mode + cat2process
      * Historic mode + cat2process
      * Historic mode + date2process
      * Historic mode + list2process
    
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
        Formatted as yyyymmdd, yyyy-mm-dd, yyyy/mm/dd or yyyy.mm.dd. When used
        with day mode: date for which the known objects database needs to be
        prepared. When used with historic mode: date for which all (light)
        transient catalogues need to be processed.
    list2process: string
        Path to and name of the text file that contains the paths to and names
        of transient catalogues (one per line) that need to be processed.
    logname: string
        Path to and name of the log file in which comments about the run are
        stored.
    redownload: boolean
        Boolean to indicate whether the asteroid (and comet) databases should be
        redownloaded. Alternatively the most recently downloaded databases will
        be used. Parameter is not relevant for the night mode.
    overwrite: boolean
        Boolean to indicate whether files will be remade and overwritten.
        Alternatively existing files will be  used.
    """
    t_glob = time.time()
    
    # Format input parameters
    mode = mode.lower()
    if date2process:
        date2process = date2process.replace(".", "").replace(
            "/", "").replace("-", "")
    
    # Set global variables
    global FILE_JPLEPH
    FILE_JPLEPH = get_par(settingsFile.JPL_ephemerisFile, tel)
    global OVERWRITE_FILES
    OVERWRITE_FILES = overwrite
    
    # Perform checks on input parameter combinations and setting file
    # parameters
    if not check_input_parameters(mode, cat2process, date2process,
                                  list2process):
        return
    if not check_settings():
        return
    folders = load_and_check_folders(tel)
    if not folders: # Empty tuple
        return
    input_folder, tmp_folder, log_folder, report_folder = folders
    
    # Logging
    setup_logfile(logname, log_folder)
    LOG.info("Mode: {}".format(mode))
    #mem_use(label="at start of run_match2SSO")
    
    # Get local noon corresponding to the night start in case date2process or
    # cat2process are specified.
    night_start = None
    local_timezone = timezone(get_par(settingsFile.timeZoneTelescope, tel))
    if cat2process:
        night_start = get_night_start_from_date(cat2process, tel)
    elif date2process:
        night_start = local_timezone.localize(datetime.strptime(" ".join([
            date2process, "120000"]), "%Y%m%d %H%M%S"))
    
    # Run match2SSO
    if mode == "day":
        day_mode(night_start, tel, tmp_folder, redownload)
        
    elif mode == "night":
        night_mode(cat2process, night_start, tel, tmp_folder, report_folder)
    
    elif mode == "historic" or "hist":
        hist_mode(cat2process, date2process, list2process, night_start, tel,
                  input_folder, tmp_folder, report_folder)
    
    LOG.info("Finished running match2SSO.")
    log_timing_memory(t_glob, label="run_match2SSO")
    logging.shutdown()
    
    return


# In[ ]:


def day_mode(night_start, tel, tmp_folder, redownload_db):
    
    """
    Run match2SSO in day mode to prepare for the night mode. The day mode
    creates the known objects database & CHK files for the upcoming night, or
    - in case [date] is specified - for the specified night.
      0) Creates a run directory in preparation of the nightly processing.
      1) Downloads asteroid and comet databases to the database folder.
      2) Integrates the asteroid database to midnight of the observation night.
      3) Combines the comet and integrated asteroid databases into a SOF-
         formatted known objects database.
      4) Creates symbolic links to the used databases and the observatory codes
         list in the run directory.
      5) Runs astcheck on a fake detection in order to create the CHK files
         that astcheck will need for faster processing when running on
         observations.
      6) Removes the fake detection in- & output, excluding the CHK files.
    Products of steps 1-2 are saved to the database folder, those of steps 3-5
    to the run directory.
    
    Parameters:
    -----------
    night_start: datetime object (incl timezone) or None
        Noon corresponding to the start of the observation night.
    tel: string
        Abbreviated telescope name. Can be either ML1, BG2, BG3 or BG4.
    tmp_folder: string
        Name of the folder where the known objects databases will be downloaded
        to and in which the run directory will be placed.
    redownload_db: boolean
        Boolean to indicate whether the asteroid (and comet) databases should be
        redownloaded. Alternatively the most recently downloaded databases will
        be used.
    """
    
    LOG.info("Running the day mode.")
    
    # Use the upcoming local night if no observation night is specified
    if not night_start:
        local_timezone = timezone(get_par(settingsFile.timeZoneTelescope, tel))
        night_start = (datetime.now(local_timezone)).strftime("%Y%m%d 120000")
        night_start = local_timezone.localize(datetime.strptime(
            night_start, "%Y%m%d %H%M%S"))
    
    # Delete the temporary products made in the night mode of the previous
    # night.
    #if not KEEP_TMP:
    #    night_previous = (night_start - timedelta(days=1)).strftime("%Y%m%d")
    #    rundir_previous = "{}{}/".format(tmp_folder, night_previous)
    #    
    #    if isfile("{}MPCORB.DAT".format(rundir_previous)):
    #        asteroid_database = os.readlink("{}MPCORB.DAT"
    #                                        .format(rundir_previous))
    #        if "epoch" in asteroid_database:
    #            os.remove(asteroid_database)
    #            LOG.info("Removed {}".format(asteroid_database))
    #    remove_tmp_folder(rundir_previous)
    
    # Create a run directory corresponding to the observation night
    rundir = ("{}{}/".format(tmp_folder, night_start.strftime("%Y%m%d")))
    LOG.info("Run directory: {}".format(rundir))
    if not os.path.isdir(rundir):
        os.makedirs(rundir)
    
    # Create symbolic link to the observatory codes list
    if not isfile("{}ObsCodes.html".format(rundir)):
        LOG.info("Creating symbolic link to ObsCodes.html")
        os.symlink(settingsFile.obsCodesFile,
                   "{}ObsCodes.html".format(rundir))
    
    # Download and integrate known object databases
    midnight = night_start + timedelta(days=0.5)
    create_known_objects_database(midnight, rundir, tmp_folder, redownload_db)
    
    # Create CHK files that astcheck needs in advance, to allow parallelisation
    # in the night mode
    create_chk_files(night_start, "night_start", rundir)
    create_chk_files(night_start + timedelta(days=1), "night_end", rundir)
    
    # Check for known object database products
    if not find_database_products(rundir):
        LOG.error("No database products found...")
        return
        
    LOG.info("Day mode finished.")
    
    return


# In[ ]:


def night_mode(cat_name, night_start, tel, tmp_folder, report_folder):
    
    """
    Run match2SSO on a single transient catalogue. The day mode should have been
    run once before the night mode. This allows the night mode to run in
    parallel on multiple transient catalogues of the same night, as steps that
    cannot be parallelised (making known objects database and CHK files) have
    already been executed in the day mode. The night mode:
      1) Converts the transient catalogue into an MPC-formatted text file.
      2) Runs astcheck on that file, to find matches between the transient
         detections and known solar system objects.
      3) Makes an SSO catalogue containing the matches.
      4) Makes an MPC report of the matches.
    Running the night mode in parallel on multiple catalogues can be done by
    calling match2SSO (with --mode night) multiple times in parallel.
    
    Parameters:
    -----------
    cat_name: string
        Path to and name of the transient catalogue that is to be processed.
    night_start: datetime object (with timezone)
        Noon corresponding to the start of the observation night.
    tel: string
        Abbreviated telescope name. Can be either ML1, BG2, BG3 or BG4.
    tmp_folder: string
        Name of the folder that contains the known objects databases and the run
        directory.
    report_folder: string
        Name of the folder in which the MPC reports will be stored.
    """
    
    LOG.info("Running the night mode on transient catalogue: \n{}"
             .format(cat_name))
    
    # Get name of run directory
    rundir = "{}{}/".format(tmp_folder, night_start.strftime("%Y%m%d"))
    LOG.info("Run directory: {}".format(rundir))
    
    # Check for known object database products. Stop processing if it doesn't
    # exist.
    if not find_database_products(rundir):
        logging.shutdown()
        return
    
    # Check for CHK files
    night_start_utc = get_night_start_from_date(cat_name, tel, "utc")
    night_end_utc = night_start_utc + timedelta(days=1)
    night_start_utc = night_start_utc.strftime("%Y%m%d")
    night_end_utc = night_end_utc.strftime("%Y%m%d")
    if not isfile("{}{}.chk".format(rundir, night_start_utc)):
        LOG.critical("Missing {}.chk!".format(night_start_utc))
        logging.shutdown()
        return
    if not isfile("{}{}.chk".format(rundir, night_end_utc)):
        LOG.critical("Missing {}.chk!".format(night_end_utc))
        logging.shutdown()
        return
    del night_start_utc
    del night_end_utc
    
    # Check for observatory codes list. Stop processing if it doesn't exist
    if not isfile("{}ObsCodes.html".format(rundir)):
        LOG.critical("{}ObsCodes.html doesn't exist.".format(rundir))
        logging.shutdown()
        return
    
    _ = match_single_catalogue(cat_name, rundir, tmp_folder, report_folder,
                               night_start, make_kod=False, redownload_db=False)
    
    # Beware that the run directory created for the processing of the
    # catalogue is not removed. This is the case because a single parallel
    # process does not know about the rest. A cleaning function can be run
    # in the day mode.
    
    return


# In[ ]:


def hist_mode(cat_name, date, catlist, night_start, tel, input_folder,
              tmp_folder, report_folder, redownload_db):
    
    """
    The historic mode does the entire processing of match2SSO, including the
    preparation of the known objects databases. The historic mode can be run on
    a single transient catalogue, an entire night of observations or a list of
    transient catalogues (possibly spanning multiple nights).
    
    For the first catalogue that is processed of each observation night, the
    code:
      0) Creates a run directory in preparation of the processing.
      1) Downloads asteroid and comet databases to the database folder.
      2) Integrates the asteroid database to midnight of the observation night.
      3) Combines the comet and integrated asteroid databases into a SOF-
         formatted known objects database.
      4) Creates symbolic links to the used databases and the observatory codes
         list in the run directory.
    The asteroid and comet databases used for this are only downloaded once per
    historic mode run (and only if [redownload_db] is True). For the
    remaining files (or if [redownload_db] is False), the most recently
    downloaded versions are used.
    
    For all transient catalogues, the code then:
      5) Converts the transient catalogue into an MPC-formatted text file.
      6) Runs astcheck on that file, to find matches between the transient
         detections and known solar system objects.
      7) Makes an SSO catalogue containing the matches.
      8) Makes an MPC report of the matches.
    
    The historic mode can only be parallelised if there is no overlap in 
    observing nights.
    
    Parameters:
    -----------
    cat_name: string or None
        Path to and name of the transient catalogue that is to be processed.
    date: string or None
        Date for which the known objects database needs to be prepared.
        Formatted as yyyymmdd.
    catlist: string or None
        Path to and name of the text file that contains the paths to and names
        of transient catalogues (one per line) that need to be processed.
    night_start: datetime object (with timezone) or None
        Noon corresponding to the start of the observation night.
    tel: string
        Abbreviated telescope name. Can be either ML1, BG2, BG3 or BG4.
    input_folder: string
        Folder which contains the yyyy/mm/dd/ folders in which the transient
        catalogues are stored.
    tmp_folder: string
        Name of the folder that contains the known objects databases and the run
        directory.
    report_folder: string
        Name of the folder in which the MPC reports will be stored.
    redownload_db: boolean
        Boolean to indicate whether the asteroid (and comet) databases should be
        redownloaded. Alternatively the most recently downloaded databases will
        be used.
    Beware that exactly two of the variables (cat_name, date, catlist) need to
    be None!
    """
    
    if catlist:
        LOG.info("Running historic mode on catalogue list: \n{}"
                 .format(catlist))
        
        # Open catalogue list
        with open(catlist, "r") as catalogue_list:
            catalogues2process = [name.strip() for name in catalogue_list                                   if name[0] != "#"]
        
        # Order by observation date (noon that equals the start of the
        # observation day)
        noons = [get_night_start_from_date(catname, tel).strftime(
            "%Y%m%d %H%M%S") for catname in catalogues2process]
        
        # Process files per night
        LOG.info("Catalogue list spans {} nights".format(len(np.unique(noons))))
        first_night = True
        for noon in np.unique(noons):
            LOG.info("Processing night that starts at {}".format(noon))
            night_index = np.where(np.array(noons) == noon)[0]
            catalogues2process_1night = np.array(
                catalogues2process)[night_index]
            
            local_timezone = timezone(get_par(settingsFile.timeZoneTelescope,
                                              tel))
            noon = local_timezone.localize(datetime.strptime(noon,
                                                             "%Y%m%d %H%M%S"))
            if first_night:
                match_catalogues_single_night(
                    catalogues2process_1night, noon, redownload_db, tmp_folder,
                    report_folder)
            else:
                match_catalogues_single_night(catalogues2process_1night, noon,
                                              False, tmp_folder, report_folder)
            first_night = False
        return
    
    if cat_name:
        LOG.info("Running historic mode on transient catalogue: \n{}"
                 .format(cat_name))
        catalogues2process = (cat_name)
    
    elif date:
        LOG.info("Running historic mode on night {}".format(date))
        catalogues2process = get_transient_filenames(
            input_folder, night_start,
            night_start+timedelta(days=1)-timedelta(minutes=1), tel)
        
        if not catalogues2process:
            LOG.critical("No (light) transient catalogues exist for night {}"
                         .format(date))
            return
    
    match_catalogues_single_night(catalogues2process, night_start,
                                  redownload_db, tmp_folder, report_folder)
    
    return


# In[ ]:


def match_catalogues_single_night(catalogues_single_night, night_start,
                                  redownload_db, tmp_folder, report_folder):
    """
    Wrapper function that calls the match_single_catalogue function for
    each catalogue in a list that contains catalogues corresponding to
    observations taken on the same night.
    
    Parameters:
    -----------
    catalogues_single_night: list of strings
        Names of the catalogues that will need to be processed by match2SSO.
        The catalogues are expected to correspond to observations taken on the
        same night.
    night_start: datetime object
        Start of the observing night.
    redownload_db: boolean
        Boolean indicating whether the known object databases need to be
        redownloaded, or alternatively if existing downloaded versions of the
        databases can be used for the matching.
    tmp_folder: string
        Name of the folder containing the known objects databases, or where
        these databases can be downloaded.
    report_folder: string
        Name of the folder in which the MPC reports will be stored.
    """
    if TIME_FUNCTIONS:
        t_func = time.time()
    
    LOG.info("{} catalogues to process for the night around {}.".format(
        len(catalogues_single_night),
        night_start.strftime("%Y-%m-%d %H:%M:%S")))
    
    #Create run directory
    rundir = "{}{}/".format(tmp_folder, night_start.strftime("%Y%m%d"))
    LOG.info("Run directory: {}".format(rundir))
    if not os.path.isdir(rundir):
        os.makedirs(rundir)
    
    #Run matching per catalogue
    make_kod = True
    for cat_name in catalogues_single_night:
        made_kod = match_single_catalogue(
            cat_name, rundir, tmp_folder, report_folder, night_start,
            make_kod, redownload_db)
        if made_kod:
            make_kod = False #Only make known objects database once
    
    #Remove the run directory after processing the last catalogue of the night
    if not KEEP_TMP:
        #Remove integrated database made for this night
        if isfile("{}MPCORB.DAT".format(rundir)):
            asteroid_database = os.readlink("{}MPCORB.DAT".format(rundir))
            if "epoch" in asteroid_database:
                os.remove(asteroid_database)
                LOG.info("Removed {}".format(asteroid_database))
        
        #Remove temporary folder made for this night
        remove_tmp_folder(rundir)
    
    if TIME_FUNCTIONS:
        log_timing_memory(t_func, label="match_catalogues_single_night")
    
    return


# In[ ]:


def match_single_catalogue(cat_name, rundir, tmp_folder, report_folder,
                           night_start, make_kod, redownload_db):
    """
    Run matching routine on a single transient catalogue. Optionally, a new
    known objects database is created where the reference epoch corresponds to
    midnight on the observation night. The detections in the transient 
    catalogue are then matched to the positions of the solar system bodies in
    the known objects catalogue. Matches are saved to an SSO catalogue.
    Function returns a boolean indicating whether a known objects catalogue was
    (successfully) made.
    
    Parameters:
    -----------
    cat_name: string
        Name of the transient catalogue of which the detections are to be
        matched to known solar system objects.
    rundir: string
        Directory corresponding to observation night where all temporary
        products will be stored during the running of match2SSO.
    tmp_folder: string
        Name of the folder containing the known objects databases, or where
        these databases can be downloaded.
    report_folder: string
        Name of the folder in which the MPC reports will be stored.
    night_start: datetime object, including time zone
        Noon corresponding to the start of the local night during which the
        observation corresponding to the transient catalogue was made.
    make_kod: boolean
        Boolean indicating whether a new known objects database needs to be
        made. In night mode, this should be False. In historic mode, this is
        only True for the first catalogue that is processed per observation
        night, as the database will need to be integrated to that observation
        night.
    redownload_db: boolean
        Only used when make_kod is True. This boolean indicates whether the
        asteroid and comet databases will need to be redownloaded before making
        the known objects database. Alternatively, the most recent, previously
        downloaded version of the databases are used.
    """
    #mem_use(label="at start of match_single_catalogue")
    if TIME_FUNCTIONS:
        t_func = time.time()
    LOG.info("Running match2SSO on {}".format(cat_name))
    
    # Check if input catalogue exists and is not flagged red
    is_existing, is_dummy, cat_name = check_input_catalogue(cat_name)
    if not is_existing:
        return made_kod
    del is_existing
    
    # File names
    mpcformat_file = "{}{}".format(rundir, os.path.basename(cat_name).replace(
        "_light.fits", ".fits").replace(".fits", "_MPCformat.txt"))
    sso_cat = cat_name.replace("_light.fits", ".fits").replace(".fits",
                                                               "_sso.fits")
    predict_cat = sso_cat.replace("_trans", "").replace("_sso.fits",
                                                        "_sso_predict.fits")
    reportname = "{}{}.txt".format(
        report_folder, os.path.basename(sso_cat).replace(".fits", "_report"))
    
    # Keep track of whether known object database has been made
    made_kod = False
    
    # Create predictions and SSO catalogues in case of a dummy (empty) detection
    # catalogue
    if is_dummy:
        _ = predictions(None, rundir, predict_cat, "")
        create_sso_catalogue(None, rundir, sso_cat, 0)
        return made_kod
    
    # Check if output catalogues exist, in which case the known objects database
    # will not need to be made
    if isfile(predict_cat) and isfile(sso_cat) and not OVERWRITE_FILES:
        LOG.info("Prediction and SSO catalogues exist and won't be remade.\n")
        
        # Check for any version of an MPC report for this observation
        reportnames = list_files(reportname.replace(".txt", ""), end_str=".txt")
        if reportnames:
            LOG.info("There is at least one version of an MPC report for this "
                     "observation. No new one will be made.\n")
            return made_kod
        
        # Check if MPC-formatted file that the create_MPC_report function needs
        # exists and get the MPC code. If it doesn't exist yet / anymore, remake
        # this file first.
        mpc_code, create_dummy = convert_fits2mpc(cat_name, mpcformat_file)
        if mpc_code is None:
            LOG.critical("Unknown MPC code - MPC report will not be made.")
            return made_kod
        
        if create_dummy:
            # No MPC report will need to be made, as the SSO catalogue is a
            # dummy (empty) catalogue.
            return made_kod
        
        # Create an report that can be used to submit the detections that were
        # matched to known solar system objects to the MPC
        create_MPC_report(sso_cat, mpcformat_file, reportname, rundir, mpc_code)
        
        if not KEEP_TMP:
            os.remove(mpcformat_file)
            LOG.info("Removed {}".format(mpcformat_file))
        return made_kod
    
    # Convert the transient catalogue to an MPC-formatted text file
    mpc_code, create_dummy = convert_fits2mpc(cat_name, mpcformat_file)
    if mpc_code is None:
        LOG.critical("Matching cannot be done because of unknown MPC code.")
        LOG.info("Creating dummy catalogues.")
        _ = predictions(None, rundir, predict_cat, "")
        create_sso_catalogue(None, rundir, sso_cat, 0)
        return made_kod
    
    if create_dummy:
        _ = predictions(None, rundir, predict_cat, "")
        create_sso_catalogue(None, rundir, sso_cat, 0)
        return made_kod
    
    # If make_kod, create a new known objects database with a reference epoch
    # corresponding to midnight of the observation night. Also create a
    # symbolic link to the MPC observatory codes list.
    if make_kod:
        midnight = night_start + timedelta(days=0.5)
        create_known_objects_database(midnight, rundir, tmp_folder,
                                      redownload_db)
        made_kod = True
        
        # Make symbolic link to observatory codes list if it doesn't exist yet
        if not isfile("{}ObsCodes.html".format(rundir)):
            LOG.info("Creating symbolic link to ObsCodes.html")
            os.symlink(settingsFile.obsCodesFile,
                       "{}ObsCodes.html".format(rundir))
    
    # Make predictions catalogue
    N_sso = predictions(cat_name, rundir, predict_cat, mpc_code)
    
    # Check if the SSO catalogue and MPC report already exist
    reportnames = list_files(reportname.replace(".txt", ""), end_str=".txt")
    if reportnames and isfile(sso_cat) and not OVERWRITE_FILES:
        LOG.info("SSO catalogue and MPC report exist and won't be remade.\n")
        if not KEEP_TMP:
            os.remove(mpcformat_file)
            LOG.info("Removed {}".format(mpcformat_file))
        return made_kod
    
    # Run astcheck on the MPC-formatted transient file
    astcheck_file = mpcformat_file.replace("_MPCformat.txt",
                                           "_astcheckMatches.txt")
    run_astcheck(mpcformat_file, rundir, astcheck_file)
    
    # Save matches found by astcheck to an SSO catalogue
    create_sso_catalogue(astcheck_file, rundir, sso_cat, N_sso)
    
    # Create a report that can be used to submit the detections that were
    # matched to known solar system objects to the MPC
    create_MPC_report(sso_cat, mpcformat_file, reportname, rundir, mpc_code)
    
    # Delete temporary files corresponding to the processed transient
    # catalogue. The other temporary files (the CHK files, the SOF file and the
    # symbolic links) in the run directory are not (yet) removed, as they
    # might be needed for processing of other data from the same night.
    if not KEEP_TMP:
        os.remove(mpcformat_file)
        LOG.info("Removed {}".format(mpcformat_file))
        os.remove(astcheck_file)
        LOG.info("Removed {}".format(astcheck_file))
    
    if TIME_FUNCTIONS:
        log_timing_memory(t_func, label="match_single_catalogue")
    
    return made_kod


# ## Core routines in match2SSO
# * Create known objects catalogue (with asteroids with a max. orbital
#   uncertainty)
# * Create catalogue with predictions of asteroids in the FOV
# * Convert transient catalogue to MPC format
# * Run matching algorithm (astcheck) on file made in previous step
# * Save solar system object detections (matches) to SSO catalogue
# * Create MPC report

# In[ ]:


def create_chk_files(noon, noon_type, rundir):
    
    """
    Function that creates the CHK files that astcheck uses (and produces if
    they don't exist yet) when matching, by running astcheck on a fake
    detection and subsequently removing the fake data products excluding the
    CHK files. These files describe the positions of all asteroids at the start
    and end of the night (both at noon) in UTC. This function is only needed
    when running match2SSO in the day mode, as this will then allow
    parallelisation in the night mode.
    
    As the 24-hour local observing night (between local noons) can overlap with
    two UTC nights (between UTC noons) if there is a large time difference
    between the timezone of the telescope and UTC, we will have to try
    producing CHK files both for a time close to the start of the local night
    (we take 1 min after) as well as for a time close to the end of it (1 min
    before). This function should hence be run twice. This will produce a total
    of 2 or 3 CHK files that astcheck will subsequently use to match any
    observation taken during the local observing night to the known solar
    system objects.
    
    Parameters:
    -----------
    noon: datetime object
        Noon that corresponds to either the start or the end of the observing
        night.
    noon_type: string
        Should be either "night_start" or "night_end", depending on which noon
        was given as input.
    rundir: string
        Directory in which astcheck is run. This directory should contain
        the mpc2sof catalogue that contains the known solar system objects.
    """
    
    # Check noon_type parameter and set observation time for fake
    # detection
    noon_type = noon_type.lower()
    if noon_type not in ["night_start", "night_end"]:
        LOG.error("Unknown noon type!")
        return
    if noon_type == "night_start":
        obstime = noon + timedelta(minutes=1)
    elif noon_type == "night_end":
        obstime = noon - timedelta(minutes=1)
    
    # Convert observation time of fake detection to UTC
    obstime = obstime.astimezone(pytz.utc)
    
    # Create MPC-formatted file with fake detection
    mpcformat_file_fake = "{}fakedetection_{}_MPCformat.txt".format(rundir,
                                                                    noon_type)
    LOG.info("Creating fake detection: {}".format(mpcformat_file_fake))
    mpcformat_file_fake_content = open(mpcformat_file_fake, "w")
    fake_detection = "".join([
        "     0000001  C{} {:0>2} {:08.5f} "
        .format(obstime.year, obstime.month, obstime.day),
        "00 00 00.00 +00 00 00.0          0.00 G      L66"
        ])
    mpcformat_file_fake_content.write(fake_detection)
    mpcformat_file_fake_content.close()
    
    # Run astcheck on fake observation to create CHK files
    LOG.info("Running astcheck on fake detection")
    astcheck_file_fake = mpcformat_file_fake.replace("_MPCformat.txt",
                                                     "_astcheckMatches.txt")
    run_astcheck(mpcformat_file_fake, rundir, astcheck_file_fake,
                 matching_radius=0)
    
    # Remove MPC-formatted file and astcheck output related to the fake
    # detection
    os.remove(mpcformat_file_fake)
    LOG.info("Removed {}".format(mpcformat_file_fake))
    os.remove(astcheck_file_fake)
    LOG.info("Removed {}".format(astcheck_file_fake))
    
    return


# In[ ]:


def use_existing_database(databaselist, sso_type):
    
    """
    Select the most recently downloaded unintegrated asteroid/comet database
    from the input list [databaselist], by the database's name. If there is no
    unintegrated database, select the most recent integrated one. Empty
    databases (created when INCLUDE_COMETS = False) should not be in the input
    list.
    The function returns the database name & version.
    
    Parameters:
    -----------
    databaselist: list of strings
        List of the asteroid / comet databases that were downloaded previously.
        This can include different database versions and/or integrated to
        different epochs.
    sso_type: string
        Type of database: asteroid or comet
    """
    LOG.info("Selecting {} database that was downloaded in a previous run."
             .format(sso_type))
    
    # Select unintegrated databases only if available and sort on name
    unintegrated_databases = [DB for DB in databaselist if "epoch" not in DB]
    databases_sorted = sorted(unintegrated_databases)
    if not databases_sorted:
        databases_sorted = sorted(databaselist)
    
    # Select most recently downloaded database
    database_name = databases_sorted[-1]
    
    # Retrieve database version from database name
    database_version = os.path.splitext(
        os.path.basename(database_name))[0].split("_")[1].replace("version", "")
    LOG.info("{} database version {}".format(sso_type, database_version))
    
    return database_name, database_version


# In[ ]:


def download_database(sso_type, redownload_db, tmp_folder):
    
    """
    Function downloads the asteroid or comet database if desired. After
    downloading, asteroids with too large uncertainties will be removed from
    the downloaded asteroid database copy.
    Alternatively, the function will load the latest downloaded version of the
    database. The function returns the database name and version number.
    
    Parameters:
    -----------
    sso_type: string
        Solar system object type that is in the database. Can be either
        'asteroid' or 'comet'. Capitals are allowed as well.
    redownload_db: boolean
        If False, the databases will not be redownloaded. Instead, the name and
        version number of the latest downloaded database version are returned.
    tmp_folder: string
        Folder to save the asteroid and comet databases to that are downloaded
        in this function.
    """
    if TIME_FUNCTIONS:
        t_subfunc = time.time()
        
    sso_type = sso_type.lower()
    if sso_type == "asteroid":
        database_url = settingsFile.URL_asteroidDatabase
    elif sso_type == "comet":
        database_url = settingsFile.URL_cometDatabase
    else:
        error_string = "Database type unknown. Cannot be downloaded."
        LOG.critical(error_string)
        raise ValueError(error_string)
    
    # Determine whether database needs to be downloaded
    existing_databases = list_files("{}{}DB_".format(tmp_folder, sso_type),
                                    end_str=".dat")
    if not redownload_db and existing_databases:
        database_name, database_version = use_existing_database(
            existing_databases, sso_type)
        return database_name, database_version
    
    # Download database and get database version
    LOG.info("Downloading {} database...".format(sso_type))
    database_version = datetime.utcnow().strftime("%Y%m%dT%H%M")
    database_name = "{}{}DB_version{}.dat".format(tmp_folder, sso_type,
                                                  database_version)
    req = requests.get(database_url, allow_redirects=True)
    open(database_name, "wb").write(req.content)
    LOG.info("{} database version: {}".format(sso_type, database_version))
    
    # Remove asteroids with large orbital uncertainties from database
    if sso_type == "asteroid":
        select_asteroids_on_uncertainty(database_name)
    
    # Compare database size to that of a previously downloaded database, as a
    # crude check whether the database was downloaded fully. As both databases
    # have objects with large orbital uncertainties removed, there may be a
    # small size difference. We therefore round to the nearest 10,000 bytes
    # before comparing.
    if existing_databases:
        oldDBname = sorted(existing_databases)[-1]
        size_oldDB = round(os.path.getsize(oldDBname)/10000.)
        size_newDB = round(os.path.getsize(database_name)/10000.)
        if size_newDB < size_oldDB:
            LOG.error("Downloaded database is {}e4 bytes, whereas the old "
                      "database {} is {}e4 bytes in size. Something went wrong!"
                      .format(size_newDB, oldDBname, size_oldDB))
            LOG.info("Removing just downloaded database")
            os.remove(database_name)
            
            # Use older database version (downloaded before)
            database_name, database_version = use_existing_database(
                existing_databases, sso_type)
            
            if TIME_FUNCTIONS:
                log_timing_memory(t_subfunc, label="downloadDatabase ({})"
                                  .format(sso_type))
            return database_name, database_version  
    
    #if not KEEP_TMP and len(existing_databases) > 0:
    #    LOG.info("Removing older {} database versions.".format(sso_type))
    #    for old_database in existing_databases:
    #        os.remove(old_database)
    #        LOG.info("Removed {}.".format(old_database))
    
    if TIME_FUNCTIONS:
        log_timing_memory(t_subfunc, label="downloadDatabase ({})"
                          .format(sso_type))
    
    return database_name, database_version


# In[ ]:


def create_known_objects_database(midnight, rundir, tmp_folder, redownload_db):
    
    """
    Function downloads the most recent versions of the asteroid database and
    the comet database. It then uses integrat.cpp from the lunar repository to
    integrate the asteroid orbits to midnight of the observation night, in
    order to optimize the predicted positions of known objects.
    
    Beware: the current version of integrat.cpp cannot be used on JPL's comet
    file, as it is not compatible with its format. As a consequence, there
    might be an offset in the predictions of the comet positions, perhaps
    causing us to miss these objects in the linking routine. Note that the
    MPC's comet file could be integrated to the right epoch using integrat.cpp,
    but as this file is not compatible with astcheck, it cannot be used here.
    
    Parameters:
    -----------
    midnight: datetime object, including time zone
        Local midnight during the observation night.
    rundir: string
        Directory in which mpc2sof is run and which the known objects catalogue
        is saved to.
    tmp_folder: string
        Folder to save the known objects databases to that are created in this
        function.
    redownload_db: boolean
        If False, the databases will not be redownloaded. They will only be
        integrated to the observation epoch (midnight on the observation
        night).
    """
    #mem_use(label="at start of create_known_objects_database")
    if TIME_FUNCTIONS:
        t_func = time.time()
    
    # Download asteroid database
    asteroid_database, asteroid_database_version = download_database(
        "asteroid", redownload_db, tmp_folder)
    
    # Download comet database if requested
    if INCLUDE_COMETS:
        _, comet_database_version = download_database("comet", redownload_db,
                                                      tmp_folder)
    else:
        LOG.info("Do not download comet database. Instead, create an empty "
                 "comet database so that there's no matching to comets.")
        open("{}ELEMENTS.COMET".format(rundir), "w").write("".join([
            "Num  Name                                     Epoch      q      ",
            "     e        i         w        Node          Tp       Ref\n---",
            "---------------------------------------- ------- ----------- ---",
            "------- --------- --------- --------- -------------- ------------"
            ]))
    
    # Integrat only accepts UTC midnights. Choose the one closest to local
    # midnight.
    date_midnight = midnight.date()
    if midnight.hour >= 12.:
        date_midnight = date_midnight + timedelta(days=1)
    midnight_utc = pytz.utc.localize(datetime.strptime(" ".join([
        date_midnight.strftime("%Y%m%d"), "000000"]), "%Y%m%d %H%M%S"))
    midnight_utc_str = midnight_utc.strftime("%Y%m%dT%H%M")
    
    # Integrate the asteroid database to the observation date
    integrated_asteroid_database = (
        "{}asteroidDB_version{}_epoch{}.dat"
        .format(tmp_folder, asteroid_database_version, midnight_utc_str))
    if not isfile(integrated_asteroid_database):
        LOG.info("Integrating asteroid database to epoch {}..."
                 .format(midnight_utc_str))
        if TIME_FUNCTIONS:
            t_subtiming = time.time()
        subprocess.run(["integrat", asteroid_database,
                        integrated_asteroid_database,
                        str(Time(midnight_utc).jd), "-f{}".format(FILE_JPLEPH)],
                       cwd=tmp_folder, stdout=subprocess.DEVNULL, check=True)
        if TIME_FUNCTIONS:
            log_timing_memory(t_subtiming, label="integrat")
        
        # Check file size (rounded to nearest 10,000 bytes) as a crude check to
        # see if integration step was performed without issues
        size_DB = round(os.path.getsize(asteroid_database)/10000.)
        size_intDB = round(os.path.getsize(integrated_asteroid_database)/10000.)
        if size_intDB < size_DB:
            LOG.error("Integrated database is smaller than original database! "
                      "({}e4 byte < {}e4 byte)".format(size_intDB, size_DB))
            LOG.warning("Using unintegrated database instead!")
            integrated_asteroid_database = asteroid_database
        
        # Remove temporary file created by integrat
        if isfile("{}ickywax.ugh".format(tmp_folder)):
            os.remove("{}ickywax.ugh".format(tmp_folder))
    
    # Create the symbolic links in the run directory that mpc2sof needs
    symlink_asteroid_database = "{}MPCORB.DAT".format(rundir)
    if isfile(symlink_asteroid_database):
        LOG.info("Removing the old MPCORB.DAT symbolic link")
        os.unlink(symlink_asteroid_database)
    os.symlink(integrated_asteroid_database, symlink_asteroid_database)
    LOG.info("Created symbolic link {}".format(symlink_asteroid_database))
    
    if INCLUDE_COMETS:
        symlink_comet_database = "{}ELEMENTS.COMET".format(rundir)
        if isfile(symlink_comet_database):
            LOG.info("Removing the old ELEMENTS.COMET symbolic link")
            os.unlink(symlink_comet_database)
        os.symlink("{}cometDB_version{}.dat".format(
            tmp_folder, comet_database_version), symlink_comet_database)
        LOG.info("Created symbolic link {}".format(symlink_comet_database))
    #mem_use(label="after creating symbolic links to the databases")
    
    # Combine the known comets and asteroids into a SOF file, which astcheck
    # will then use as input
    LOG.info("Combining asteroids and comets into SOF file.")
    if TIME_FUNCTIONS:
        t_subtiming = time.time()
    
    subprocess.run("mpc2sof", cwd=rundir, check=True)
    
    if TIME_FUNCTIONS:
        log_timing_memory(t_subtiming, label="mpc2sof")
        log_timing_memory(t_func, label="create_known_objects_database")
    LOG.info("Finished loading and formatting external databases.")
    
    return


# In[ ]:


def select_asteroids_on_uncertainty(asteroid_database):
    
    """
    Go through the asteroid database (MPCORB format) and select the asteroids
    which have orbital uncertainty parameters smaller than maxUncertainty.
    The MPC uncertainty parameters that we consider are explained here:
    https://www.minorplanetcenter.net/iau/info/UValue.html
    
    Overwrite the database with just the asteroids selected on their orbital
    uncertainties, so that there will be no matching with poorly known objects.
    
    Parameters:
    -----------
    asteroid_database: string
        Name of the full asteroid database (MPCORB-formatted text-file).
    """
    #mem_use(label="at start of select_asteroids_on_uncertainty")
    
    if settingsFile.maxUncertainty is None:
        LOG.info("All known solar system bodies are used in the matching, "
                 "irrespective of their uncertainty parameter.")
        return
    
    if TIME_FUNCTIONS:
        t_func = time.time()
    LOG.info("Removing asteroids with too large uncertainties...\n")
    
    #Open asteroid database
    asteroid_database_content = open(asteroid_database, "r").readlines()
    
    #Find the size of the header of the asteroid database, assuming that the
    #header ends with a line of dashes.
    header_end_index = 0
    line_index = 0
    for line in asteroid_database_content:
        if line.startswith("-----"):
            header_end_index = line_index
            break
        line_index += 1
    
    #Re-write the asteroid database file, including only the header and the
    #lines corresponding to asteroids that have small orbital uncertainties.
    number_asteroids_pre_selection = 0
    number_asteroids_post_selection = 0
    with open(asteroid_database, "w") as asteroid_database_content_new:
        for line_index in range(len(asteroid_database_content)-1):
            
            #Copy header to file
            if line_index <= header_end_index:
                asteroid_database_content_new.write(
                    asteroid_database_content[line_index])
                continue
            
            line = asteroid_database_content[line_index]
            
            #Copy empty lines
            if line == "\n":
                asteroid_database_content_new.write(line)
                continue
            
            number_asteroids_pre_selection += 1
            
            #Filter on uncertainty parameter. Copy lines of asteroids for
            #which orbits are determined reasonably well.
            uncertainty = line[105]
            if uncertainty.isdigit():
                if float(uncertainty) <= settingsFile.maxUncertainty:
                    asteroid_database_content_new.write(line)
                    number_asteroids_post_selection += 1
    
    LOG.info("{} out of {} asteroids have U <= {}".format(
        number_asteroids_post_selection, number_asteroids_pre_selection,
        settingsFile.maxUncertainty))
    LOG.info("Asteroid database now only includes sources with U <= {}"
             .format(settingsFile.maxUncertainty))
    if TIME_FUNCTIONS:
        log_timing_memory(t_func, label="select_asteroids_on_uncertainty")
    
    return


# In[ ]:


def predictions(transient_cat, rundir, predict_cat, mpc_code,
                is_FOV_circle=False):
    """
    Use astcheck to predict which asteroids are in the FOV during the
    observation. Predictions can be made for a circular or square FOV. The
    function returns the number of asteroids that are estimated to be bright
    enough to be detected (V mag <= limiting AB magnitude). This is a crude
    estimate, as no correction for the different type of magnitudes is taken
    into account.
    
    Parameters:
    -----------
    transient_cat: string or None
        Path to and name of the transient catalogue. If None, we will make a
        dummy (empty) predictions catalogue.
    rundir: string
        Directory corresponding to observation night where the temporary output
        file that astcheck makes will be stored.
    predict_cat: string
        Path to and name of the output catalogue to which the predictions will
        be saved.
    mpc_code: string
        MPC code corresponding to the telescope.
    is_FOV_circle: boolean
        Boolean indicating if the FOV is circular. If False, a square FOV is
        assumed where the sides are aligned with RA & Dec.
    """
    #mem_use(label="at start of predictions")
    
    if TIME_FUNCTIONS:
        t_func = time.time()
    
    if not OVERWRITE_FILES and isfile(predict_cat):
        LOG.info("Prediction catalogue already exists and won't be re-made.")
        with fits.open(predict_cat) as hdu:
            hdr = hdu[1].header
        return hdr["N-SSO"]
    
    # If the transient catalogue was red-flagged, we will not make a predictions
    # catalogue as there is little use for it.
    if transient_cat is None:
        LOG.info("Creating a dummy predictions catalogue.")
        sso_header = create_sso_header(rundir, 0, 0, True, False)
        fitstable = format_cat(Table(), start_header=sso_header)
        save_fits(fitstable, predict_cat, rundir=rundir)
        if TIME_FUNCTIONS:
            log_timing_memory(t_func, label="predictions")
        return 0
    
    LOG.info("Making predictions...")
    
    # Open header of transient catalogue
    with fits.open(transient_cat) as hdu:
        hdr = hdu[1].header
    
    # Get the observation date (isot format) and central field coordinates (in
    # deg) from the header
    date = hdr[DATE_KEYWORD]
    if CENTRAL_RA_KEYWORD in hdr.keys() and CENTRAL_DEC_KEYWORD in hdr.keys():
        ra_field = hdr[CENTRAL_RA_KEYWORD]
        dec_field = hdr[CENTRAL_DEC_KEYWORD]
    else:
        ra_field = hdr["RA"]
        dec_field = hdr["DEC"]
    limmag = hdr[LIMMAG_KEYWORD]
    
    # Create temporary output file for astcheck results
    output_file = "{}{}".format(rundir, os.path.basename(transient_cat).replace(
        "_light", "").replace("_trans.fits", "_sso_predictions.txt"))
    output_file_content = open(output_file, "w")
    
    if is_FOV_circle:
        field_radius = 3600.*settingsFile.FOV_width/2.
    else:
        #Square FOV. Use the half diagonal of the FOV as the radius.
        field_radius = 3600.*sqrt(2)*settingsFile.FOV_width/2.
    
    # Run astcheck from folder containing .sof-file
    subprocess.call(["astcheck", "-c", str(date), str(ra_field), str(dec_field),
                     str(mpc_code), "-r{}".format(field_radius), "-h",
                     "-m{}".format(settingsFile.limitingMagnitude),
                     "-M{}".format(settingsFile.maximalNumberOfAsteroids)],
                    stdout=output_file_content, cwd=rundir)
    output_file_content.close()
    
    # Read in astcheck output
    astcheck_file_content = open(output_file, "r").readlines()
    astcheck_file_content = remove_astcheck_header_and_footer(
        astcheck_file_content)
    
    # Create table to store solar system bodies and their properties
    output_columns = {
        "ID_SSO":       ["12a", "", ""],
        "RA_SSO":       ["f4", "deg", "%.6f"],
        "DEC_SSO":      ["f4", "deg", "%.6f"],
        "V_RA_SSO":     ["f4", "arcsec/hour", "%.4f"],
        "V_DEC_SSO":    ["f4", "arcsec/hour", "%.4f"],
        "MAG_V_SSO":    ["f4", "", "%.2f"]
        }
    output_table = Table()
    for key in output_columns.keys():
        output_table.add_column(Column(name=key, dtype=output_columns[key][0],
                                       unit=output_columns[key][1]))
        if output_columns[key][2]:
            output_table[key].format = output_columns[key][2]
        
    # Loop through SSOs and store their properties in the output table
    for source in astcheck_file_content:
        source = re.sub("\n$", "", source) # Remove line end character
        source_properties = re.split(" +", source)
        if source_properties[0]:
            identifier = " ".join([source_properties[0], source_properties[1]])
        else:
            identifier = source_properties[1]
        ra_source, dec_source, magnitude, v_ra, v_dec = source_properties[2:]
        
        output_row = (identifier, float(ra_source), float(dec_source),
                      float(v_ra), float(v_dec), float(magnitude))
        output_table.add_row(output_row)
    
    # Remove temporary astcheck output file
    if not KEEP_TMP:
        os.remove(output_file)
        LOG.info("Removed {}".format(output_file))
    
    # In case of a square FOV, disregard SSOs outside the FOV
    if not is_FOV_circle and len(output_table)>0:
        center = SkyCoord(ra_field, dec_field, unit="deg", frame="icrs")
        sources = SkyCoord(output_table["RA_SSO"],
                           output_table["DEC_SSO"], unit="deg", frame="icrs")
        dra, ddec = center.spherical_offsets_to(sources)
        mask_dist = ((abs(dra.deg) <= settingsFile.FOV_width/2.) &
                     (abs(ddec.deg) <= settingsFile.FOV_width/2.))
        output_table = output_table[mask_dist]
    
    # Create header for predicted asteroids in FOV
    i_bright = np.where(output_table["MAG_V_SSO"] <= limmag)[0]
    N_sso = len(i_bright)
    if N_sso > 0:
        dummy = False
    else:
        dummy = True
    sso_header = create_sso_header(rundir, 0, N_sso, dummy, False)
    
    # Save to table
    fitstable = format_cat(output_table, start_header=sso_header)
    save_fits(fitstable, predict_cat, rundir=rundir)
    
    LOG.info("Predictions saved to {}.".format(predict_cat))
    if TIME_FUNCTIONS:
        log_timing_memory(t_func, label="predictions")
    
    return N_sso


# In[ ]:


def create_sso_header(rundir, N_det, N_sso, dummy, incl_detections):
    
    """
    Function creates the headers for the SSO catalogue and the SSO predictions
    catalogue.
    
    Parameters:
    -----------
    rundir: string
        Name of the folder in which the symbolic links to the databases are
        stored. These are used to get the version numbers of the databases.
    N_det: int
        Number of detected solar system objects. Only one matched object is
        counted per transient detection and if an object is matched to multiple
        transients, it is also only counted once.
    N_sso: int
        Number of solar system objects in the FOV that are supposedly bright
        enough for a detection (V magnitude < T-LMAG). The difference between V
        and AB magnitudes is ignored here. This number is therefore a rough
        prediction for the number of detections.
    dummy: boolean
        Boolean indicating whether the catalogue is a dummy catalogue without
        sources (dummy=True). If False, there are sources in the catalogue.
    incl_detections: boolean
        Boolean indicating whether the catalogue is the SSO catalogue -
        corresponding to detected solar system objects - or not. If not, it is
        the catalogue containing the predicted objects in the FOV and some of
        the header keywords will not be included in the header.
    """
    #mem_use(label="at start of create_sso_header")
    LOG.info("Creating SSO header.")
    if TIME_FUNCTIONS:
        t_func = time.time()
    
    # Create empty SSO header
    header = fits.Header()
    
    # Add Python version to SSO header
    header["PYTHON-V"] = (platform.python_version(), "Python version used")
    
    # Get C++ version and add to the SSO header. Based on 
    # [https://stackoverflow.com/questions/44734397/which-c-standard-is-the-
    # default-when-compiling-with-g/44735016#44735016]
    proc = subprocess.run("g++ -dM -E -x c++  /dev/null | grep -F __cplusplus",
                          capture_output=True, shell=True, check=True)
    cpp_macro = proc.stdout.decode("utf-8").replace("\n", "").split()[-1]
    
    if cpp_macro not in settingsFile.CPPmacro2version.keys():
        LOG.error("C++ macro unknown: {}".format(cpp_macro))
        cpp_version = "None"
    else:
        cpp_version = settingsFile.CPPmacro2version[cpp_macro]
    header["CPP-V"] = (cpp_version, "C++ version used")
    
    # Get G++ version and add to SSO header
    proc = subprocess.run("g++ --version", capture_output=True, shell=True,
                          check=True)
    gpp_version = proc.stdout.decode("utf-8").split("\n")[0].split()[-1]
    header["GPP-V"] = (gpp_version, "G++ version used")
    
    # Add match2SSO & header keyword versions to the SSO header
    header["SSO-V"] = (__version__, "match2SSO version used")
    header["SSOKW-V"] = (KEYWORDS_VERSION,
                         "SSO header keywords version used")
    
    # Get the versions of the lunar & jpl_eph repositories. The versions are
    # given by unique strings signifying the latest commit that was made to the
    # repositories. Save the strings to the SSO header.
    lunar_version = retrieve_version("lunar")
    header["LUNAR-V"] = (lunar_version, "lunar repository version used")
    
    jpl_eph_version = retrieve_version("jpl_eph")
    header["JPLEPH-V"] = (jpl_eph_version, "jpl_eph repository version used")
    
    # Add version of JPL lunar & planetary ephemerides file to SSO header
    header["JPLDE-V"] = ("DE{}".format(FILE_JPLEPH.split(".")[-1]),
                         "JPL ephemeris file version used")
    
    # Add asteroid database version & reference epoch to the SSO header.
    # The MPCORB.DAT symbolic link in the run directory refers to the
    # asteroid database version that was used. The name structure of this
    # database is: 
    # asteroid_database_version[yyyymmddThhmm]_epoch[yyyymmddThhmm].dat
    asteroid_database_version, asteroid_database_epoch = "None", "None"
    if isfile("{}MPCORB.DAT".format(rundir)):
        asteroid_database = os.readlink("{}MPCORB.DAT".format(rundir))
        asteroid_database_date = os.path.basename(
            asteroid_database).split("_")[1].replace("version", "")
        asteroid_database_version = "{}-{}-{}T{}:{}".format(
            asteroid_database_date[0:4], asteroid_database_date[4:6],
            asteroid_database_date[6:8], asteroid_database_date[9:11],
            asteroid_database_date[11:13])
        reference_epoch = os.path.basename(asteroid_database).split(
            "_")[2].replace("epoch", "")
        asteroid_database_epoch = ("{}-{}-{}T{}:{}".format(
            reference_epoch[0:4], reference_epoch[4:6], reference_epoch[6:8],
            reference_epoch[9:11], reference_epoch[11:13]))
    
    header["ASTDB-V"] = (asteroid_database_version,
                         "asteroid database version (date in UTC)")
    header["ASTDB-EP"] = (asteroid_database_epoch,
                          "asteroid database epoch in UTC")
    
    # Add comet database version to the SSO header.
    # The ELEMENTS.COMET symbolic link in the run directory refers to the
    # comet database version that was used. The name structure of this
    # database is: cometDB_version[yyyymmddThhmm].dat
    comet_database_version = "None"
    if INCLUDE_COMETS:
        comet_database = os.readlink("{}ELEMENTS.COMET".format(rundir))
        comet_database_date = os.path.basename(comet_database).split(
            "_")[1].replace("version", "")
        comet_database_version = "{}-{}-{}T{}:{}".format(
            comet_database_date[0:4], comet_database_date[4:6],
            comet_database_date[6:8], comet_database_date[9:11],
            comet_database_date[11:13])
    
    header["COMDB-V"] = (comet_database_version,
                         "comet database version (date in UTC)")
    
    # Add matching radius and maximum orbital uncertainty parameter to header
    if incl_detections:
        header["RADIUS"] = (float(settingsFile.matchingRadius),
                            "matching radius in arcsec")
    header["U-MAX"] = (settingsFile.maxUncertainty,
                       "maximum orbital uncertainty parameter")
    
    # Add number of (predicted) detections to header
    header["N-SSO"] = (N_sso, "predicted number of bright SSOs in FOV")
    if incl_detections:
        header["N-SSODET"] = (N_det, "number of detected solar system objects")
    
    # Add keyword indicating whether there is
    header["SDUMCAT"] = (bool(dummy), "dummy SSO catalogue without sources?")
    
    LOG.info("SSO header complete.")
    if TIME_FUNCTIONS:
        log_timing_memory(t_func, label="create_sso_header")
    
    return header


# In[ ]:


def convert_fits2mpc(transient_cat, mpcformat_file):
    
    """
    Function converts the transient catalogue to a text file of the MPC
    80-column format, so that astcheck can run on it. For the asteroid / comet
    identifier used in the MPC file, the transient number is used. This
    transient number cannot be used for MPC reports as it is not all-time unique
    (per telescope). But it is a straight-forward way to link detections to
    known solar system objects within match2SSO.
    
    Negative transients are excluded from the MPC-formatted file. In the case of
    moving objects, these are sources that are present in the reference image
    but not in the new image. They can be recognised as having a negative
    signal-to-noise ratio value in the significance (Scorr) image for MeerLICHT
    and BlackGEM. As reference images can be stacked images (that are not
    centered on the asteroid) for which the observation date and time is
    unclear, and as we use the date of the new image as the observation date,
    asteroids in the reference images are not useful.
    
    Function returns the MPC Observatory code of the telescope with which the
    observation was made, and a boolean indicating whether there were detections
    to be converted (dummy = False) or not (dummy = True).
    
    Parameters:
    -----------
    transient_cat: string
        Path to and name of the transient catalogue.
    mpcformat_file: string
        Path to and name of the MPC-formatted text file that is made in this
        function.
    """
    #mem_use(label="at start of convert_fits2mpc")
    LOG.info("Converting transient catalogue to MPC-format.")
    if TIME_FUNCTIONS:
        t_func = time.time()
    
    # Load transient catalogue header
    with fits.open(transient_cat) as hdu:
        transient_header = hdu[1].header
    
    # Get the MPC observatory code from the header
    mpc_code = transient_header[MPC_CODE_KEYWORD].strip()
    if mpc_code not in list(pd.read_fwf(settingsFile.obsCodesFile,
                                        widths=[4, 2000],
                                        skiprows=1)["Code"])[:-1]:
        LOG.critical("MPC code {} is not in the MPC list of observatory codes"
                     .format(mpc_code))
        return None, True
    
    # Check if MPC-formatted file exists and if it should be overwritten or not
    if not OVERWRITE_FILES and isfile(mpcformat_file):
        LOG.info("MPC-formatted file already exists and will not re-made.")
        return mpc_code, False
    
    # Get observation date in the right format
    date_obs = Time(transient_header[DATE_KEYWORD], format="isot").datetime
    decimal_day = date_obs.day + 1./24.*(
        date_obs.hour + 1./60.*(date_obs.minute + 1./60.*(
            date_obs.second + date_obs.microsecond/10.**6)))
    mpc_char16to32 = "{} {:0>2} {:08.5f} ".format(date_obs.year, 
                                                  date_obs.month, decimal_day)
    # Load transient catalogue data
    with fits.open(transient_cat) as hdu:
        detections = Table(hdu[1].data)
    
    # Remove negative transients
    index_positives = np.where(detections[SNR_COLUMN]>=0)[0]
    detections = detections[index_positives]
    
    # Check if there are positive transients to include
    if len(detections) == 0:
        dummy = True
        LOG.info("No (positive) sources available for linking.")
        return mpc_code, dummy
    dummy = False
    
    # Create output file
    mpcformat_file_content = open(mpcformat_file, "w")
    
    # Loop over the detections and add data to the MPC-formatted file
    for detection_index in range(len(detections)):
        # Use the source numbers as "temporary designations" in the MPC format.
        # In this way, we will be able to link the known objects to the right
        # source.
        line = ("     {:0>7}  C"
                .format(detections[NUMBER_COLUMN][detection_index]))
        line = "".join([line, mpc_char16to32])
        
        # Get the coordinates and magnitude of the source
        coord = SkyCoord(detections[RA_COLUMN][detection_index],
                         detections[DEC_COLUMN][detection_index],
                         unit="deg", frame="icrs") 
        mag = "{:.1f}".format(detections[MAG_COLUMN][detection_index])
        
        line = "".join([line, "{} {}          {} G      {}"
                        .format(coord.to_string("hmsdms", sep=" ",
                                                precision=2)[:11],
                                coord.to_string("hmsdms", sep=" ",
                                                precision=1)[-11:],
                                mag.rjust(4), mpc_code)])
        
        # Write the data to the MPC-formatted file
        mpcformat_file_content.write("{}\n".format(line))
    
    mpcformat_file_content.close()
    
    LOG.info("MPC-formatted file saved to {}.".format(mpcformat_file))
    if TIME_FUNCTIONS:
        log_timing_memory(t_func, label="convert_fits2mpc")
    
    return mpc_code, dummy


# In[ ]:


def run_astcheck(mpcformat_file, rundir, output_file,
                 matching_radius=settingsFile.matchingRadius):
    """
    Run astcheck on the input transient catalogue to find matches between
    transient detections and known solar system objects. Per detection, all
    matches within the matching_radius are selected and saved to the output text
    file.
    
    Parameters:
    -----------
    mpcformat_file: string
        Name of the input text-file that is formatted according to the MPC's
        80-column MPC format. This file can list detections / tracklets for
        astcheck to match, but it can also contain just a single row
        representing the observation. In the latter use case, the specified
        coordinates should correspond to the centre of the observation.
    rundir: string
        Directory in which astcheck is run. This directory should contain
        the mpc2sof catalogue that contains the known solar system objects.
    output_file: string
        Path to and name of the output text file in which the matches found
        by astcheck are stored.
    matching_radius: int or float
        Matching radius in arcsec. The default value is the one specified in
        the settings file [set_match2SSO.py].
    """
    #mem_use(label="at start of run_astcheck")
    LOG.info("Running astcheck: matching detections to known solar system "
             "bodies.")
    if TIME_FUNCTIONS:
        t_func = time.time()
    
    if not OVERWRITE_FILES and isfile(output_file):
        LOG.info("Astcheck output file already exists and won't be re-made.")
        return
    
    # Create a file for storing the output of the astcheck run
    output_file_content = open(output_file, "w")
    
    # Run astcheck from folder containing .sof-file
    subprocess.call(["astcheck", mpcformat_file, "-h",
                     "-r{}".format(matching_radius),
                     "-m{}".format(settingsFile.limitingMagnitude),
                     "-M{}".format(settingsFile.maximalNumberOfAsteroids)],
                    stdout=output_file_content, cwd=rundir)
    output_file_content.close()
    
    LOG.info("Matches saved to {}.".format(output_file))
    if TIME_FUNCTIONS:
        log_timing_memory(t_func, label="run_astcheck")
    
    return


# In[ ]:


def create_sso_catalogue(astcheck_file, rundir, sso_cat, N_sso):
    
    """
    Open the text-file that was produced when running astcheck [astcheck_file]
    and save the information to an SSO catalogue.
    
    Parameters:
    -----------
    astcheck_file: string
        Name of the file containing astcheck's output (the matches). This can
        be None, in which case we will create a dummy SSO catalogue without
        matches.
    rundir: string
        Directory in which astcheck was run. This directory also contains
        symbolic links to the asteroid and comet databases that were used to
        create the known objects catalogue that astcheck used.
    sso_cat: string
        Name of the SSO catalogue to be created in which the matches
        are stored.
    N_sso: int
        Number of solar system objects in the FOV that are supposedly bright
        enough for a detection (V magnitude < T-LMAG). The difference between V
        and AB magnitudes is ignored here. This number is therefore a rough
        prediction for the number of detections.
    """
    #mem_use(label="at start of create_sso_catalogue")
    LOG.info("Converting astcheck output into an SSO catalogue.")
    if TIME_FUNCTIONS:
        t_func = time.time()
    
    if not OVERWRITE_FILES and isfile(sso_cat):
        LOG.info("SSO catalogue already exists and will not be re-made.")
        return
    
    # If the transient catalogue was red-flagged, matching was not performed
    # and an empty SSO catalogue needs to be created.
    if astcheck_file is None:
        LOG.info("Creating a dummy SSO catalogue.")
        sso_header = create_sso_header(rundir, 0, N_sso, True, True)
        fitstable = format_cat(Table(), start_header=sso_header)
        save_fits(fitstable, sso_cat, rundir=rundir)
        if TIME_FUNCTIONS:
            log_timing_memory(t_func, label="create_sso_catalogue")
        return
    
    # Remove astcheck header and footer if needed
    astcheck_file_content = open(astcheck_file, "r").readlines()
    astcheck_file_content = remove_astcheck_header_and_footer(
        astcheck_file_content)
    
    # Find empty lines
    separator = "\n"
    indices_separator = np.where(np.array(astcheck_file_content)
                                 == separator)[0]
    indices_separator = np.append(-1, indices_separator)
    indices_separator = np.append(indices_separator,
                                  len(astcheck_file_content))
    
    # Create table to store match information in
    output_columns = {
        NUMBER_COLUMN:  ["i4", ""],
        "ID_SSO":       ["12a", ""],
        "DIST_RA_SSO":  ["i2", "arcsec"],
        "DIST_DEC_SSO": ["i2", "arcsec"],
        "DIST_SSO":     ["i2", "arcsec"],
        "MAG_V_SSO":    ["f4", ""],
        "FLAGS_SSO":    ["i2", ""]
        }
    output_table = Table()
    for key in output_columns.keys():
        output_table.add_column(Column(name=key, dtype=output_columns[key][0],
                                       unit=output_columns[key][1]))
    # Loop over sources
    N_det = 0
    for index in range(len(indices_separator)-1):
        minimal_index = indices_separator[index]+1
        maximal_index = indices_separator[index+1]
        
        if minimal_index == maximal_index:
            continue
        
        # Name of the source in the MPC-formatted input file (transient number)
        transient_number = astcheck_file_content[minimal_index:maximal_index][
            0].split(":")[0].split()
        # Lines corresponding to matches in the astcheck output file
        matches = astcheck_file_content[minimal_index:maximal_index][1:]
        
        if not matches: #Empty list
            continue
        
        N_det += 1
        
        # If a source is matched to multiple solar system objects, assign the
        # matches a flag of 1
        if len(matches)>1:
            initial_flag = 1
        else:
            initial_flag = 0
        
        # Get properties of matches
        for i_match in range(len(matches)):
            match_properties = re.split(" +", matches[i_match])
            match_properties = [x for x in match_properties if len(x) > 0]
            identifier = match_properties[0]
            try:
                int(match_properties[1])
                offset_ra, offset_dec, offset, magnitude = match_properties[1:5]
            except ValueError:
                identifier = " ".join([identifier, match_properties[1]])
                offset_ra, offset_dec, offset, magnitude = match_properties[2:6]
            
            try:
                magnitude = float(magnitude)
            except ValueError:
                LOG.warning("Magnitude '{}' could not be converted to float."
                            .format(magnitude))
                magnitude = None
        
            # Add match to output table
            output_row = (transient_number, str(identifier), float(offset_ra),
                          float(offset_dec), float(offset), magnitude, initial_flag)
            output_table.add_row(output_row)
    
    # If a solar system object was matched to multiple transient sources in the
    # image, assign it a flag of 2
    unique_objects = np.unique(output_table["ID_SSO"])
    if len(unique_objects) != len(output_table):
        for obj in unique_objects:
            obj_indices = np.where(output_table["ID_SSO"] == obj)[0]
            if len(obj_indices) > 1:
                output_table["FLAGS_SSO"][obj_indices] += 2
    
    # Set dummy parameter (dummy means that there were no matches found)
    dummy = False
    if not output_table: #Empty table
        dummy = True
    
    # Create header for SSO catalogue and save catalogue
    sso_header = create_sso_header(rundir, N_det, N_sso, dummy, True)
    fitstable = format_cat(output_table, start_header=sso_header)
    save_fits(fitstable, sso_cat, rundir=rundir)
    
    LOG.info("Matches saved to SSO catalogue: {}".format(sso_cat))
    if TIME_FUNCTIONS:
        log_timing_memory(t_func, label="create_sso_catalogue")
    
    return


# In[ ]:


def create_MPC_report(sso_cat, mpcformat_file, reportname, rundir, mpc_code):
    
    """
    Make an MPC report using the SSO catalogue and the MPC-formatted file that
    were created within match2SSO to link the transient detections from a single
    catalogue to known solar system objects. The detections corresponding to
    matches are grouped in a 'known objects report'. The identifiers used in the
    report are the packed designations of the matching objects. These are the
    packed permanent designations if available. Otherwise, the packed
    provisional designations are used.
    
    The MPC report will be compiled in a temporary folder and the complete file
    will be moved to the reports folder at the end of the function.
    
    Parameters:
    -----------
    sso_cat: string
        Path to and name of the SSO catalogue of which the matches need to be
        converted to a MPC report.
    mpcformat_file: string
        Name of the MPC formatted file that was made in match2SSO for the
        matching (but does not contain the correct SSO identifiers yet for
        reporting to the MPC).
    reportname: string
        Name of the MPC report that will be made from the SSO catalogue. The
        report should have the extension ".txt".
    rundir: string
        Run directory in which the MPC report will be created, before being
        moved to its proper destination as given in the [sso_cat] path.
    mpc_code: string
        MPC code corresponding to the telescope.
    """
    #mem_use(label="at start of create_MPC_report")
    LOG.info("Creating MPC report...")
    if TIME_FUNCTIONS:
        t_func = time.time()
    
    # Compose temporary report name using run directory as path
    destination = os.path.dirname(reportname)
    report_basename = os.path.basename(reportname).replace(
        ".txt", "_{}.txt".format(Time.now().strftime("%Y%m%dT%H%M%S")))
    reportname = rundir + report_basename
    destination_file = "/".join([destination, report_basename])
    
    # Check if file already exists (will only happen when running this function
    # multiple times in close succession, as the production time is used in the
    # file name)
    if not OVERWRITE_FILES and isfile(destination_file):
        LOG.info("MPC report already exists and will not be re-made.")
        return
    
    # Open SSO catalogue
    with fits.open(sso_cat) as hdu:
        sso_cat_content = Table(hdu[1].data)
    
    # Check SSO catalogue for matches
    if not sso_cat_content:
        LOG.info("No matches found. MPC report will not be made.")
        return
    
    # Create MPC report
    LOG.info("Writing report {}".format(reportname))
    report_content = open(reportname, "w")
    
    # Write header to the report
    report_content.write(create_report_header(reportname, mpc_code))
    
    # Open MPC-formatted file as the MPC report will be very similar, only with
    # MPC designations rather than transient numbers as the first column. For a
    # large part, the MPC-formatted file content will therefore be copied over.
    detections_mpcformat = pd.read_fwf(mpcformat_file, widths=[14, 66],
                                       names=["char1to14", "char15to80"],
                                       dtype={"char1to14":np.int32, 
                                              "char15to80":str})
    
    # For each detection that was matched to a known solar system object,
    # get the packed designation of the matching object and write the detection
    # to the report
    for match_index in range(len(sso_cat_content[NUMBER_COLUMN])):
        designation = sso_cat_content["ID_SSO"][match_index].strip()
        
        # Start creating the line of the report corresponding to the detection,
        # by adding the packed designation of the object to the line
        if re.match(r"^[0-9]{4}\s[A-Z]", designation) or "/" in designation:
            # Provisional or survey designation
            packed_designation = wrapper_pack_provisional_designation(
                designation)
            if packed_designation is None:
                continue
            detection_line = "    {}  ".format(packed_designation)
        
        else:
            # Asteroid or comet with permanent designation
            packed_designation, fragment = pack_permanent_designation(
                designation)
            if packed_designation is None:
                continue
            detection_line = "{}{}  ".format(packed_designation,
                                             fragment.rjust(7))
        
        # Get the detection details from the MPC-formatted file and add to the
        # line of the MPC report corresponding to the detection
        detection_index = np.where(
            np.array(detections_mpcformat["char1to14"])
            == int(sso_cat_content[NUMBER_COLUMN][match_index]))[0]
        if len(detection_index) != 1:
            LOG.error("{} detections found that correspond to transient number"
                      " {}. Should be only one.".format(len(detection_index),
                      int(sso_cat_content[NUMBER_COLUMN][match_index])))
            continue
        detection_index = detection_index[0]
        detection_line = "".join([detection_line, detections_mpcformat[
            "char15to80"][detection_index]])
        
        # Check line corresponding to detection and write to the MPC report if
        # all is well
        if len(detection_line) != 80:
            LOG.error("Detection not formatted correctly in 80 columns:\n{}"
                      .format(detection_line))
        report_content.write(detection_line+"\n")
    
    report_content.close()
    
    # Move report from run directory to final destination
    LOG.info("Moving report to {}".format(destination))
    copy_file(reportname, destination_file, move=True)
    if isfile(destination_file):
        LOG.warning("MPC report {} is overwritten.".format(destination_file))
    else:
        LOG.info("MPC report saved to {}".format(destination_file))
    
    if TIME_FUNCTIONS:
        log_timing_memory(t_func, label="create_MPC_report")
    
    return


# In[ ]:


def create_report_header(reportname, mpc_code, comment=None):
    
    """
    Function composes the header of the MPC report corresponding to a single
    transient catalogue.
    
    Parameters:
    -----------
    reportname: string
        Name of the MPC report for which the header is composed.
    mpc_code: string
        MPC code of the telescope with which the observation was made.
    comment: string
        Comment to be added to the header in the COM line. By default, this is
        None, meaning that the COM line is not added to the header.
    """
    #mem_use(label="at start of create_report_header")
    LOG.info("Creating header for MPC report...")
    if TIME_FUNCTIONS:
        t_func = time.time()
    
    firstline = "COD {}\n".format(mpc_code)
    mainheader = get_par(settingsFile.MPCreportHeader, mpc_code)
    
    # Special cases for which a phrase needs to be included in the ACK line
    # of the header of the MPC report:
    # neocand = "NEO CANDIDATE" #submitting new NEO candidate
    # neocp = "NEOCP"           #submitting observations of NEOCP objects
    
    # Add ACK line to the header of the MPC report.
    ack_line = "ACK {}\n".format(Path(reportname).stem)
    if len(ack_line) > 82:
        LOG.error("ACK line in report {} is too long!".format(reportname))
    
    # Add COM line to the header
    com_line = ""
    if comment is not None:
        if len(comment) > 76:
            LOG.warning("COM line is too long and therefore not used. "
                        "Use at most 76 characters!")
        else:
            com_line = "COM {}\n".format(comment)
    
    LOG.info("MPC report header complete.")
    if TIME_FUNCTIONS:
        log_timing_memory(t_func, label="create_report_header")
    
    return "".join([firstline, mainheader, ack_line, com_line])


# ## Subroutines
# The functions below are subroutines, which we split up in
# - subroutines to support selecting input catalogues
# - subroutines to pack SSO designations
# - general subroutines for checks and set-up
# - subroutines to support Google Cloud bucket systems
# - subroutines for formatting & saving output and removing temporary files

# ### Subroutines to support selecting input catalogues
# Support selecting the correct input transient catalogues to run match2SSO on

# In[ ]:


def get_transient_filenames(input_folder, minimal_date, maximal_date, tel,
                            exclude_flagged=False):
    """
    Function returns a list with the transient file names that were taken
    between the minimal and maximal specified dates with the specified [tel]
    telescope. This function works for MeerLICHT and BlackGEM data and runs on
    the lighter version of the transient catalogues if available. Otherwise it
    runs on the 'standard' transient catalogues. If exclude_flagged is True, the
    dummy transient catalogues (which are red-flagged) are excluded.
    
    Function assumes a directory and filename structure and hence might not be
    directly applicable to other telescopes than MeerLICHT & BlackGEM.
    
    Parameters:
    -----------
    input_folder: string
        Folder which contains the yyyy/mm/dd/ folders in which the transient
        catalogues are stored.
    minimal_date: datetime object, incl time zone
        Minimal observation date of the time block for which the observations
        are selected.
    maximal_date: datetime object, incl time zone
        Maximal observation date of the time block for which the observations
        are selected.
    tel: string
        Telescope abbreviation.
    exclude_flagged: boolean
        Boolean indicating whether red-flagged (dummy) catalogues should be
        excluded or not.
    """
    #mem_use(label="at start of get_transient_filenames")
    LOG.info("Selecting transient catalogues between {} and {}."
             .format(minimal_date.strftime("%Y-%m-%d %H:%M:%S"),
                     maximal_date.strftime("%Y-%m-%d %H:%M:%S")))
    if TIME_FUNCTIONS:
        t_func = time.time()
    
    # Convert to local time
    local_timezone = timezone(get_par(settingsFile.timeZoneTelescope, tel))
    minimal_date = minimal_date.astimezone(local_timezone)
    maximal_date = maximal_date.astimezone(local_timezone)
    
    # Select the transient directories from observation date
    search_path = input_folder
    if minimal_date.year == maximal_date.year:
        yeardir = "%d/"%(minimal_date.year)
        search_path += yeardir
        if minimal_date.month == maximal_date.month:
            monthdir = "{:0>2}/".format(minimal_date.month)
            search_path += monthdir
            
            # Get observing date (defined to start at noon and end at noon the
            # next day)
            minday = minimal_date.day
            if minday.hour < 12.:
                minday -= 1
            maxday = maximal_date.day
            if maxday.hour < 12.:
                maxday -= 1
            
            if minday == maxday:
                daydir = "{:0>2}/".format(maxday)
                search_path += daydir
    
    # Infer folder depth to use
    if 'gs://' in input_folder:
        bucket_name, __ = get_bucket_name(input_folder)
        folder_tmp = input_folder.split(bucket_name)[-1]
        depth = folder_tmp.rstrip('/').count('/')
    else:
        depth = input_folder.rstrip('/').count('/')
    depth += 3 # Adds 3 to the depth for the yyyy/mm/dd folders
    
    # Use function [list_folders] to list all reduced data folders that will be
    # searched for transient catalogues
    list_paths = list_folders(search_path, depth=depth)
    
    # List transient catalogues in the selected folders. List the lighter
    # version of the catalogue if it exists. Otherwise list the 'standard'
    # version.
    transient_files = []
    for reddir in list_paths:
        lighttransnames = list_files(reddir, end_str="_trans_light.fits")
        transient_files.extend(lighttransnames)
        
        transnames = list_files(reddir, end_str="_trans.fits")
        for transname in transnames:
            if transname.replace(".fits", "_light.fits") not in lighttransnames:
                transient_files.append(transname)
    transient_files = sorted(transient_files)
    
    if not transient_files:
        return []
    
    files2process = []
    for transient_cat in transient_files:
        # Parse date encoded in filename and compare with our limits
        # (e.g. ML1_20200517_034221_red_trans_light.fits)
        splitted_filename = os.path.basename(transient_cat).split("_")
        date_obs = splitted_filename[1]
        time_obs = splitted_filename[2]
        observation_time = Time.strptime(date_obs+time_obs, "%Y%m%d%H%M%S").mjd
        if (observation_time >= Time(minimal_date).mjd and
                observation_time <= Time(maximal_date).mjd):
            with fits.open(transient_cat) as hdu:
                header = hdu[1].header
            
            if not exclude_flagged:
                files2process.append(transient_cat)
            else:
                LOG.info("Excluding red-flagged (dummy) catalogues.")
                
                if DUMMY_KEYWORD not in header.keys():
                    LOG.critical("{} not in the header!".format(DUMMY_KEYWORD))
                    return []
                
                if not header[DUMMY_KEYWORD]:
                    files2process.append(transient_cat)
    
    LOG.info("{} transient catalogues have been selected."
             .format(len(files2process)))
    if TIME_FUNCTIONS:
        log_timing_memory(t_func, label="get_transient_filenames")
    return files2process


# In[ ]:


def get_night_start_from_date(cat_name, tel, noon_type="local"):
    
    """
    This function returns the noon corresponding to the start of the
    observation night, as a datetime object. This is either the local noon or
    the noon in UTC, as specified. The noon is deduced from the catalogue
    header.
    
    Parameters:
    -----------
    cat_name: string
        Name of the catalogue corresponding to an observation that took place
        on the observation night for which the noon that signifies the start of
        the night must be determined.
    tel: string
        Telescope abbreviation.
    noon_type: string
        Must be either "local" or "utc". If "utc", this function will return
        the noon corresponding to the start of the night in UTC. This can be
        different from the local noon.
    """
    noon_type = noon_type.lower()
    
    # Get observation time from catalogue header and define as being in UTC
    with fits.open(cat_name) as hdu:
        hdr = hdu[1].header
    
    observation_time = pytz.utc.localize(Time(hdr[DATE_KEYWORD],
                                              format="isot").datetime)
    observation_date = str(observation_time.date())
    
    if noon_type == "local":
        local_timezone = timezone(get_par(settingsFile.timeZoneTelescope, tel))
        
        # Get local noon corresponding to the start of the observing night
        local_noon = local_timezone.localize(datetime.strptime(" ".join([
            observation_date, "120000"]), "%Y-%m-%d %H%M%S"))
        # Get date of observing night
        if observation_time < local_noon:
            date = (observation_time - timedelta(days=1)).strftime("%Y%m%d")
        else:
            date = observation_time.strftime("%Y%m%d")
        
        # Make local noon variable
        night_start = local_timezone.localize(datetime.strptime(" ".join([
            date, "120000"]), "%Y%m%d %H%M%S"))
    else:
        if noon_type != "utc":
            LOG.error("Noon type not understood. Assuming noon in utc.")
        
        night_start = pytz.utc.localize(datetime.strptime(" ".join([
            observation_date, "120000"]), "%Y-%m-%d %H%M%S"))
        if int(observation_time.hour) < 12.:
            night_start -= timedelta(days=1)
    
    return night_start


# ### Subroutines for SSO designation packing
# Pack asteroid & comet designations for the MPC reports

# In[ ]:


def pack_permanent_designation(full_designation):
    
    """
    Function converts the permanent minor planet designation into its packed
    form (5 characters), using the definitions given in
    https://www.minorplanetcenter.net/iau/info/PackedDes.html#perm
    Return the packed designation and - if applicable - the letter
    corresponding to the comet fragment. If the object is not a comet fragment,
    an empty string will be returned for the fragment letter.
    
    Parameters:
    -----------
    full_designation: string
        Unpacked permanent designation assigned to the object by the MPC.
    """
    fragment = ""
    
    if not full_designation.isdigit():
        # Object is a comet
        if len(full_designation.split("-")) == 2:
            designation, fragment = full_designation.split("-")
        else:
            designation = full_designation
            fragment = ""
        packed_designation = designation.zfill(5)
        fragment = fragment.lower()
    
    elif int(full_designation) < 99999:
        packed_designation = "{:0>5}".format(int(full_designation))
    
    elif int(full_designation) < 620000:
        quotient = int(full_designation)//10000
        packed_designation = "{}{:0>4}".format(abbreviate_number(quotient),
                                               int(full_designation)%10000)
    else:
        remainder = int(full_designation) - 620000
        quotient3 = remainder//62**3
        remainder -= quotient3*62**3
        quotient2 = remainder//62**2
        remainder -= quotient2*62**2
        quotient1 = remainder//62
        remainder -= quotient1*62
        packed_designation = "~{}{}{}{}".format(abbreviate_number(quotient3),
                                                abbreviate_number(quotient2),
                                                abbreviate_number(quotient1),
                                                abbreviate_number(remainder))
    # Final check
    if len(packed_designation) != 5:
        LOG.error("Packed permanent designation is of incorrect length: '{}'"
                  .format(packed_designation))
        return None, ""
    
    return packed_designation, fragment


# In[ ]:


def wrapper_pack_provisional_designation(full_designation):
    
    """
    Wrapper function for packing provisional minor planet designations into
    their packed forms. The function first checks whether [full_designation] is
    a survey designation. If so, it is packed using the definitions in
    https://www.minorplanetcenter.net/iau/info/PackedDes.html#prov
    
    If the designation is not a survey designation, the function determines
    whether it belongs to a comet or an asteroid. For comets, we call the
    function pack_provisional_designation_comet. For asteroids, we call
    pack_provisional_designation_asteroid.
    
    The function returns an 8-character long string, spanning columns 5-12 in
    the MPC report, or None in case of an issue.
    
    Parameters:
    -----------
    full_designation: string
        Unpacked provisional designation assigned to the object by the MPC.
    """
    
    # Remove space before or after designation
    full_designation = full_designation.strip()
    
    # There are four special survey designation forms (for surveys that were
    # undertaken between 1960 and 1977) that should be packed differently
    survey_strings = ("P-L", "T-1", "T-2", "T-3")
    for survey_string in survey_strings:
        if (re.match(r"^[0-9]{4}\s[A-Z]", full_designation)
                and full_designation.endswith(survey_string)):
            packed_designation = "{}S{}".format(survey_string.replace("-", ""),
                                                full_designation[:4])
            return packed_designation
    
    pack_year = {18: "I", 19: "J", 20: "K"}
    
    # For a comet
    if "/" in full_designation:
        packed_designation = pack_provisional_designation_comet(
            full_designation, pack_year)
        return packed_designation
    
    # For an asteroid
    packed_designation = pack_provisional_designation_asteroid(
        full_designation, pack_year)
    
    return packed_designation


# In[ ]:


def pack_provisional_designation_asteroid(full_designation, pack_year):
    
    """
    Function converts the provisional asteroid designation into its packed form
    using the definitions given in
    https://www.minorplanetcenter.net/iau/info/PackedDes.html#prov
    We add a space in front so that the returned string is 8 characters long,
    spanning columns 5-12 in the MPC report. The function returns None in
    case of an issue.
    
    Parameters:
    -----------
    full_designation: string
        Unpacked provisional designation assigned to the asteroid by the MPC.
    pack_year: dict
        Dictionary needed to pack the first two numbers of the year (indicating
        the century) into a single letter, according to the MPC standard.
    """
    
    if int(full_designation[:2]) not in pack_year.keys():
        LOG.error("Provisional designation of asteroid {} cannot be packed. "
                  "Skipping it.".format(full_designation))
        return None
    
    packed_year = "{}{}".format(pack_year[int(full_designation[:2])],
                                full_designation[2:4])
    packed_designation = (" {}{}{}{}"
                          .format(packed_year, full_designation[5],
                                  pack_cycle_number(full_designation[7:]),
                                  full_designation[6]))
    # Final check
    if len(packed_designation) != 8:
        LOG.error("Packed provisional designation is of incorrect length: "
                  "'{}'".format(packed_designation))
        return None
    
    return packed_designation


# In[ ]:


def pack_provisional_designation_comet(full_designation, pack_year):
    
    """
    Function converts the provisional comet designation into its packed form
    using the definitions given in
    https://www.minorplanetcenter.net/iau/info/PackedDes.html#prov
    As described in https://www.minorplanetcenter.net/iau/info/OpticalObs.html,
    for comets a character is added in front of the provisional designation (at
    column 5), describing the comet type.
    The function returns an 8-character long string, spanning columns 5-12 in
    the MPC report, or None in case of an issue.
    
    Parameters:
    -----------
    full_designation: string
        Unpacked provisional designation assigned to the comet by the MPC.
    pack_year: dict
        Dictionary needed to pack the first two numbers of the year (indicating
        the century) into a single letter, according to the MPC standard.
    """
    
    comet_type, designation = full_designation.split("/")
    
    # In case of a comet fragment, the last character of the packed designation
    # is the fragment letter. Otherwise, it is zero.
    fragment = "0"
    if "-" in designation:
        designation, fragment = designation.split("-")
        fragment = fragment.lower()
    year, remainder = designation.split(" ")
        
    if int(year[:2]) not in pack_year.keys():
        LOG.error("Provisional designation of comet {} cannot be packed. "
                  "Skipping it.".format(full_designation))
        return None
    
    packed_year = "{}{}".format(pack_year[int(year[:2])], year[2:])
    
    # In case there are two letters after the space in the designation. This
    # can be the case if the object was thought to be an asteroid early on.
    if remainder[1].isalpha():
        
        if fragment != "0":
            # A comet with two letters in its provisional designation after
            # the space and a fragment letter cannot be submitted in the old
            # report format. It can in the new ADES format, but we are
            # not yet using this. Skip detection.
            LOG.error("Provisional designation of comet {} cannot be packed."
                      "Skipping it.".format(full_designation))
            return None
        
        # Although this object is not a fragment, its provisional designation
        # does contain a second letter after the space which should be written
        # to the same position as the fragment letter.
        fragment = remainder[1]
        remainder = "{}{}".format(remainder[0], remainder[2:])
    
    # There should be at most three digits after the space-letter combination
    # in the provisional designation.
    if len(remainder) > 4:
        LOG.error("Unclear how to pack provisional designation of comet {}. "
                  "Skipping it.".format(full_designation))
        return None
    
    if int(year[:2]) not in pack_year.keys():
        LOG.error("Data from before 1800 or after 2099 cannot be assigned a "
                  "packed provisional designation.")
        return None
    
    packed_designation = ("{}{}{}{}{}"
                          .format(comet_type, packed_year, remainder[0],
                                  pack_cycle_number(remainder[1:]), fragment))
    # Final check
    if len(packed_designation) != 8:
        LOG.error("Packed provisional designation is of incorrect length: "
                  "'{}'".format(packed_designation))
        return None
    
    return packed_designation


# In[ ]:


def abbreviate_number(num):
    
    """
    Number packing function needed to pack MPC designations.
    """
    
    num_dict = {str(index): letter for index, letter in 
                enumerate("".join([ascii_uppercase, ascii_lowercase]), start=10)}
    
    if int(num) > 9:
        return num_dict[str(num)]
    
    return num


# In[ ]:


def pack_cycle_number(number_of_cycles):
    
    """Input parameter number_of_cycles is a string of 0-3 digits."""
    
    if not number_of_cycles: #Empty string
        return "00"
    
    if int(number_of_cycles) > 99:
        return "{}{}".format(abbreviate_number(number_of_cycles[0:2]),
                             number_of_cycles[2])
    
    return "{:0>2}".format(number_of_cycles)


# ### General subroutines for checks and set-up
# Subroutines that support match2SSO's functionality by checking parameters, files, folders and software versions; setting up logging and time & memory tracking and more

# In[ ]:


def check_input_parameters(mode, cat2process, date2process, list2process):
    
    """
    Check if the correct (combination of) input parameters was/were defined for
    run_match2SSO. If so, this function returns True. Otherwise, False is
    returned.
    """
    all_good = True
    
    # General checks on the input parameters
    param_list = (date2process, cat2process, list2process)
    param_num_none = sum([isinstance(par, type(None)) for par in param_list])
    
    if  param_num_none < len(param_list)-1:
        print("CRITICAL: either specify --date, --catalog OR --catlist. A "
              "combination is not allowed.")
        all_good = False
    
    if mode not in ["day", "night", "historic", "hist"]:
        print("CRITICAL: unknown mode.")
        all_good = False
    
    # Checks per mode
    if mode == "day":
        if cat2process is not None or list2process is not None:
            print("CRITICAL: the day mode cannot be combined with the "
                  "--catalog or --catlist arguments.")
            all_good = False
    
    elif mode == "night":
        if cat2process is None:
            print("CRITICAL: --catalog needs to be specified when running "
                  "match2SSO in night mode.")
            all_good = False
        if date2process is not None or list2process is not None:
            print("CRITICAL: the night mode cannot be combined with the "
                  "--date or --catlist arguments.")
            all_good = False
    
    elif (mode == "historic" or "hist") and param_num_none == len(param_list):
        print("CRITICAL: --date, --catalog and --catlist are all None. Nothing"
              " to process.")
        all_good = False
    
    # Check on the existence of the specified input
    if cat2process is not None and not isfile(cat2process):
        print("CRITICAL: the specified catalog does not exist.")
        all_good = False
    
    if list2process is not None and not isfile(list2process):
        print("CRITICAL: the specified catalog list does not exist.")
        all_good = False
    
    return all_good


# In[ ]:


def load_and_check_folders(tel):
    
    """
    Function loads the folders specified in the settings file and checks
    whether they end with a slash. In addition, checks on the existence of the
    folders are performed. A tuple of the folder names is returned. The returned
    tuple is empty if there was an issue.
    
    Parameters:
    -----------
    tel: string
        Telescope abbreviation.
    """
    
    # Load folders
    input_folder = get_par(settingsFile.inputFolder, tel)
    tmp_folder = get_par(settingsFile.tmpFolder, tel)
    log_folder = get_par(settingsFile.logFolder, tel)
    report_folder = get_par(settingsFile.MPCreportFolder, tel)
    
    # Check if folder names end with a slash
    input_folder = check_folder_name(input_folder)
    tmp_folder = check_folder_name(tmp_folder)
    log_folder = check_folder_name(log_folder)
    report_folder = check_folder_name(report_folder)
    
    # Check if critical folders exists. If not, return an empty list.
    if not isdir(input_folder):
        print("CRITICAL: input folder given in settings file doesn't exist")
        return ()
    
    # Create the other folders if they don't exist, except for the log folder
    # as that is taken care of in the setup_logfile function.
    if not os.path.isdir(tmp_folder):
        os.makedirs(tmp_folder)
    if report_folder[0:5] != 'gs://' and not isdir(report_folder):
        os.makedirs(report_folder)
    
    folders = (input_folder, tmp_folder, log_folder, report_folder)
    
    return folders


# In[ ]:


def check_folder_name(directory_name):
    
    """
    Function checks if directory name ends with a slash. If not, it is added.
    """
    
    # If name is None or an empty string, don't check folder name
    if not directory_name:
        return directory_name
    
    if directory_name[-1] != "/":
        directory_name = "".join([directory_name, "/"])
    
    return directory_name


# In[ ]:


def check_settings():
    
    """Function checks the parameters in the settings file for validity."""
    
    # Check that astcheck parameters are numbers
    if not isinstance(settingsFile.matchingRadius, (float, int)):
        print("CRITICAL: incorrectly specified matching radius in settings "
              "file. Must be float or integer.")
        return False
    
    if not isinstance(settingsFile.limitingMagnitude, (float, int)):
        print("CRITICAL: incorrectly specified limiting mag. in settings "
              "file.")
        return False
    
    if not isinstance(settingsFile.maximalNumberOfAsteroids, int):
        print("CRITICAL: incorrectly specified max. number of asteroids in "
              "settings file.")
        return False
    
    if (not isinstance(settingsFile.maxUncertainty, int) and
            settingsFile.maxUncertainty is not None):
        print("CRITICAL: incorrectly specified max. uncertainty in settings "
              "file. Must be 0-9 or None.")
        return False
    
    # Check if JPL ephemeris file exists
    if not isfile(FILE_JPLEPH):
        print("CRITICAL: JPL ephemeris file specified in settings file doesn't"
              " exist.")
        return False
    
    return True


# In[ ]:


def setup_logfile(logname, log_folder):
    
    """
    Function creates log file and configures the log handler.
    """
    
    if logname is None:
        return
    
    log_dir, log_filename = os.path.split(logname)
    
    # If no folder is specified, use the log folder from the settings file
    if not log_dir:
        if not log_folder:
            print("Warning: no log folder specified. Log will be written to"
                  " directory from which match2SSO was run!")
            log_folder="."
        log_dir = log_folder
    
    # Create folder to store log in, if it does not yet exist
    if log_dir[0:5] != 'gs://' and not isdir(log_dir):
        os.makedirs(log_dir)
    
    # Configure log handling
    log_file = "{}/{}".format(log_dir, log_filename)
    if isfile(log_file):
        file_path_and_name, extension = os.path.splitext(log_file)
        log_file = "{}_{}{}".format(file_path_and_name,
                                    Time.now().strftime("%Y%m%d_%H%M%S"),
                                    extension)
        print("Log file already exists. Creating a new log named {}"
              .format(log_file))
    
    file_handler = logging.FileHandler(log_file, "a")
    file_handler.setFormatter(LOG_FORMATTER)
    file_handler.setLevel("INFO")
    LOG.addHandler(file_handler)
    
    return


# In[ ]:


# ZOGY function
def log_timing_memory(t_in, label=""):
    
    """Function to report the time and memory spent in a function."""
    
    LOG.info("wall-time spent in {}: {:.3f} s".format(label, time.time()-t_in))
    #mem_use(label=label)
    
    return


# In[ ]:


# ZOGY function
def mem_use(label=""):
    
    """Function keeps track of the memory usage."""
    
    # ru_maxrss is in units of kilobytes on Linux; however, this seems
    # to be OS dependent as on mac os it is in units of bytes; see
    # manpages of "getrusage"
    if sys.platform == "darwin":
        norm = 1024**3
    else:
        norm = 1024**2
    
    mem_max = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/norm
    mem_now = psutil.Process().memory_info().rss / 1024**3
    mem_virt = psutil.Process().memory_info().vms / 1024**3
    
    LOG.info("memory use [GB]: rss={:.3f}, maxrss={:.3f}, vms={:.3f} in {}"
             .format(mem_now, mem_max, mem_virt, label))
    
    return


# In[ ]:


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


# In[ ]:


# ZOGY function
def get_par(par, tel):
    
    """
    Function to check if [par] is a dictionary with one of the keys being [tel]
    or the alphabetic part of [tel] (e.g. 'BG'), and if so, return the
    corresponding value. Otherwise just return the parameter value.
    """
    
    par_val = par
    if isinstance(par, dict):
        if tel in par:
            par_val = par[tel]
        else:
            # cut off digits from [tel]
            tel_base = "".join([char for char in tel if char.isalpha()])
            if tel_base in par:
                par_val = par[tel_base]
    
    return par_val


# In[ ]:


def check_input_catalogue(cat_name):
    
    """
    Check if the input catalogue exists and if it is a dummy (red-flagged)
    catalogue or not. If a light version of the catalogue is available, use
    that version. This function returns a boolean for "does the catalogue 
    exist?", a boolean for "is the catalogue a dummy?" and the catalogue name
    is returned, as the light version might have been selected instead of the
    transient catalogue that includes the thumbnails.
    """
    
    # Check whether the (light) catalogue exists and ensure the use of the
    # light version of the catalogue if it is available (better in terms of
    # memory usage & processing speed)
    if "_light" not in cat_name:
        light_cat = cat_name.replace(".fits", "_light.fits")
        if isfile(light_cat):
            cat_name = light_cat
    
    if not isfile(cat_name):
        LOG.critical("The specified catalog does not exist:\n{}"
                     .format(cat_name))
        return False, None, cat_name
    
    # Check quality control flag of the catalogue
    with fits.open(cat_name) as hdu:
        header = hdu[1].header
    
    if DUMMY_KEYWORD not in header.keys():
        LOG.critical("{} not in the header of {}!".format(DUMMY_KEYWORD,
                                                          cat_name))
        return False, None, cat_name
    
    if header[DUMMY_KEYWORD]:
        LOG.info("{} is a dummy catalogue.".format(cat_name))
        return True, True, cat_name
    
    return True, False, cat_name


# In[ ]:


def find_database_products(rundir):
    
    """
    This function checks if the database products that astcheck needs in order
    to process transient catalogues are located in the run directory. These are
    the known objects catalogue, the symbolic link to the asteroid catalogue
    (needed for reading out the asteroid database version) and ELEMENTS.COMET,
    which is either an empty comet database (if comets should not be included
    in the matching) or a symbolic link to the comet database used (again 
    needed to read out the version). The function also indirectly checks for
    the existence of the run directory. A boolean is returned: True if all is
    well and False if something is missing.
    
    Parameters:
    -----------
    rundir: string
        Directory in which astcheck will be run.
    """
    
    # Check for known objects catalogue
    if not isfile("{}mpcorb.sof".format(rundir)):
        LOG.critical("The known objects catalogue (SOF format) could not be "
                     "found.")
        return False
    
    # Check for symbolic links pointing to the used version of the SSO
    # databases
    if not isfile("{}MPCORB.DAT".format(rundir)):
        LOG.critical("MPCORB.DAT could not be found")
        return False
    if not isfile("{}ELEMENTS.COMET".format(rundir)):
        LOG.critical("ELEMENTS.COMET could not be found")
        return False
    
    return True


# In[ ]:


def remove_astcheck_header_and_footer(astcheck_file_content):
    
    """
    Before the -h switch was implemented in astcheck, the header and footer
    needed to be removed manually. Check if astcheck has done this already, or
    if manual removal is required. Return the content of the astcheck file
    excluding the header and footer.
    """
    
    # The footer is variable in terms of the number of lines it spans, but it
    # always starts with the footer_string as defined below and can hence be
    # recognized by this string.
    footer_string = "The apparent motion and arc length"
    header_size = 5 # Number of header lines
    footer_index = [index for index in range(len(astcheck_file_content)) if                    footer_string in astcheck_file_content[index]]
    if footer_index:
        return astcheck_file_content[header_size:footer_index[0]]
    
    else:
        return astcheck_file_content[1:] # Just remove first empty line


# In[ ]:


def retrieve_version(package_name):
    
    """
    Function retrieves the string corresponding to the version of the software
    package, from the versions.txt file listed in the settings file. In the
    versions file, each line should correspond to a different package and should
    contain the package name and the version string, separated by a space.
    
    Parameters:
    -----------
    package_name: string
        Name of the software package, as listed in versions.txt
    """
    
    try:
        versions_file = open(settingsFile.versionsFile)
        for line in versions_file:
            if package_name in line:
                line = line.replace("\n","")
                vers = line.split(" ")[1]
                break
        versions_file.close()
    except:
        vers = None
    
    return vers


# ### Subroutines that support the use of Google Cloud systems
# The subroutines below make match2SSO compatible with a Google Cloud bucket system,
# whilst also keeping the compatibility with a Unix type file system.

# In[ ]:


# ZOGY function (for compatibility with Google Cloud system)
def get_bucket_name(path):
    
    """
    Infer bucket- and filename from [path], which is expected to be
    gs://[bucket name]/some/path/file or [bucket name]/some/path/file;
    if [path] starts with a forward slash, empty strings will be returned
    """
    
    bucket_name = path.split('gs://')[-1].split('/')[0]
    if len(bucket_name) > 0:
        # N.B.: returning filename without the starting '/'
        bucket_file = path.split(bucket_name)[-1][1:]
    else:
        bucket_file = ''
    
    return bucket_name, bucket_file


# In[ ]:


# ZOGY function (for compatibility with Google Cloud system)
def isdir(folder):
    
    if folder[0:5] == 'gs://':
        if folder[-1] != '/':
            folder = '{}/'.format(folder)
        storage_client = storage.Client()
        bucket_name, bucket_file = get_bucket_name(folder)
        blobs = storage_client.list_blobs(bucket_name, prefix=bucket_file,
                                          max_results=1)
        nblobs = len(list(blobs))
        return nblobs > 0
    
    else:
        return os.path.isdir(folder)


# In[ ]:


# ZOGY function (for compatibility with Google Cloud system)
def isfile(filename):
    
    if filename[0:5] == 'gs://':
        storage_client = storage.Client()
        bucket_name, bucket_file = get_bucket_name(filename)
        # N.B.: bucket_file should not start with '/'
        bucket = storage_client.bucket(bucket_name)
        blob = bucket.blob(bucket_file)
        return blob.exists()

    else:
        return os.path.isfile(filename)


# In[ ]:


# ZOGY function (for compatibility with Google Cloud system)
def list_files(path, search_str='', end_str='', start_str=None,
               recursive=False):
    """
    Function to return a list of files starting with [path] (can be a folder or
    google cloud bucket name and path; this does not have to be a precise or
    complete folder/path, e.g. [path] can be some_path/some_file_basename), with
    possible [end_str] and containing [search_str] without any wildcards. If
    [path] is an exact folder, then [start_str] can be used as the beginning of
    the filename.
    Beware: if [end_str] is an empty string, subfolders will also be listed
    among the files.
    """
    
    # Is google cloud being used?
    google_cloud = (path[0:5] == 'gs://')
    
    # Split input [path] into folder_bucket and prefix; if path is a folder, we
    # need to add a slash at the end. Otherwise the prefix will be the name of
    # the deepest folder
    if isdir(path) and path[-1] != '/':
        path = '{}/'.format(path)
    
    # If path is indeed a folder, then prefix will be an empty string, which is
    # fine below
    folder_bucket, prefix = os.path.split(path.split('gs://')[-1])
    
    # If path consists of just the bucket name including gs:// or just a path
    # without any slashes at all, the above will lead to an empty folder_bucket
    # and prefix=path; turn these around
    if len(folder_bucket) == 0:
        folder_bucket = prefix
        prefix = ''
    
    # If prefix is empty and [start] is defined, use
    # that as the prefix
    if prefix == '' and start_str is not None:
        prefix = start_str
    
    # If not dealing with google cloud buckets, use glob
    if not google_cloud:
        if recursive:
            files = glob.glob('{}/**/{}*{}*{}'.format(folder_bucket, prefix,
                                                      search_str, end_str),
                              recursive=True)
            if path in files:
                files.remove(path)
        else:
            files = glob.glob('{}/{}*{}*{}'.format(folder_bucket, prefix,
                                                   search_str, end_str))
    else:
        # For buckets, use storage.Client().list_blobs; see
        # https://cloud.google.com/storage/docs/samples/storage-list-files-with-prefix#storage_list_files_with_prefix-python
        # setting delimiter to '/' restricts the results to only the
        # files in a given folder
        if recursive:
            delimiter = None
        else:
            delimiter = '/'
        
        # Bucket name and prefix (e.g. gs://) to add to output files
        bucket_name, bucket_file = get_bucket_name(path)
        bucket_prefix = path.split(bucket_name)[0]
        
        if False:
            LOG.info('folder_bucket: {}'.format(folder_bucket))
            LOG.info('prefix: {}'.format(prefix))
            LOG.info('path: {}'.format(path))
            LOG.info('bucket_name: {}'.format(bucket_name))
            LOG.info('bucket_file: {}'.format(bucket_file))
        
        # Get the blobs
        storage_client = storage.Client()
        blobs = storage_client.list_blobs(bucket_name, prefix=bucket_file,
                                          delimiter=delimiter)
        
        # Collect blobs' names in list of files
        files = []
        for blob in blobs:
            filename = blob.name
            
            # Check for search string; if search string is empty, the following
            # if statement will be False
            if search_str not in filename:
                continue
            
            # Check if filename ends with [end_str]. If not, continue with next
            # blob
            len_ext = len(end_str)
            if len_ext > 0 and filename[-len_ext:] != end_str:
                continue
            
            # After surviving above checks, append filename including the bucket
            # prefix and name
            files.append('{}{}/{}'.format(bucket_prefix, bucket_name, filename))
    
    #LOG.info('number of files returned by [list_files]: {}'.format(len(files)))
    
    return files


# In[ ]:


# ZOGY function (for compatibility with Google Cloud system)
def list_folders(path, search_str='', depth=None):

    """
    Function to return list of existing folders starting with [path] (can be a
    folder or google cloud bucket name and path; this does not have to be a
    precise or complete folder/path, e.g. [path] can be 
    some_path/some_file_basename), optionally containing [search_str] without
    any wildcards. [depth] determines how many folder depths are returned,
    starting from the root directory; for google cloud buckets the root folder
    (depth=1) is considered to start after the bucket name,
    i.e. gs://[bucket name]/[root folder]/etc. If [depth] is not defined, all
    existing folders are returned.
    """
    
    # Use function [list_files] to list all files containing 'search_str'
    file_list = list_files(path, search_str=search_str, recursive=True)
    
    # Running in google cloud?
    google_cloud = (path[0:5]=='gs://')
    
    # Remove potential 'gs://' and bucket name
    if google_cloud:
        bucket_name, __ = get_bucket_name(path)
        
        # Chop off 'gs://[bucket_name]'
        file_list = [fn.split(bucket_name)[-1] for fn in file_list]
    
    if depth is None:
        # List all files and folders in path
        file_list = ['/'.join(fn.split('/')[:-1]) for fn in file_list]
    else:
        # Create list of files and folders with depth [depth]
        file_list = ['/'.join(fn.split('/')[:depth+1]) for fn in file_list
                     if len(fn.rstrip('/').split('/')) > depth]
    
    # Remove duplicate folders and sort
    file_list = sorted(list(set(file_list)))
    
    # Append 'gs://[bucket_name]' if needed
    if google_cloud:
        file_list = ['gs://{}{}'.format(bucket_name, fn) for fn in file_list]
    
    # Only folders, not files
    isfolder = [isdir(fn) for fn in file_list]
    folder_list = np.array(file_list)[isfolder].tolist()
    
    LOG.info('folder_list: {}'.format(folder_list))
    
    return folder_list


# In[ ]:


# BlackBOX function (for compatibility with Google Cloud system)
def copy_file(src_file, dest, move=False):

    """
    Function to copy or move a file [src_file] to [dest], which may be a file or
    folder; [src_file] and/or [dest] may be part of the usual filesystem or in a
    google cloud bucket; in the latter case the argument(s) should start with
    gs://[bucket_name]
    """
    
    if move:
        label = 'moving'
    else:
        label = 'copying'
    LOG.info('{} {} to {}'.format(label, src_file, dest))
    
    # If not dealing with google cloud buckets, use shutil.copy2 or shutil.move
    if not (src_file[0:5] == 'gs://' or dest[0:5] == 'gs://'):
        if not move:
            shutil.copy2(src_file, dest)
        else:
            shutil.move(src_file, dest)
    
    else:
        # This could be done in python, but much easier with gsutil from the
        # shell
        if move:
            cp_cmd = 'mv'
        else:
            cp_cmd = 'cp'
        
        cmd = ['gsutil', '-q', cp_cmd, src_file, dest]
        result = subprocess.run(cmd)
    
    return


# ### Subroutines for formatting & saving output and removing temporary files
# Also includes functionality for Google cloud file systems

# In[ ]:


def format_cat(data, header=None, header_keys=[], start_header=None):
    
    """
    Function formats the output data, composes the header and combines the two
    into a hdu fits table. The fits table is returned.
    
    Parameters:
    -----------
    data : table data
        Astropy table data which is to be used as the data for the output fits
        table.
    header: header or None
        Header from which certain keywords are copied to the header of the
        output catalogue.
    header_keys: list of strings
        Contains names of header keywords from the header mentioned above, that
        will be included in the output catalogue header. List can be empty, in
        which case the fits header is either the start_header or empty.
    start_header: header
        Header which will be included in the header of the output catalogue.
        start_header can be None.
    """
    #mem_use(label="at start of format_cat")
    
    # Format fits table data
    columns = []
    for column_name in data.columns:
        
        column_format = data[column_name].dtype
        
        # Converting bytestring format to fits format does not work properly
        # for strings, as the length is not taken into account properly.
        # Manually correct this.
        if "S" in column_format.str:
            string_length = column_format.str.split("S")[-1]
            column_format = "{}A".format(string_length)
        
        column_unit = str(data[column_name].unit)
        if column_unit == "None":
            column_unit = ""
        
        column = fits.Column(name=column_name, format=column_format,
                             unit=column_unit, array=data[column_name])
        columns.append(column)
    
    # Compose fits table header
    if start_header is not None:
        finalheader = start_header
    else:
        finalheader = fits.Header()
    
    for key in header_keys:
        finalheader[key] = (header[key], header.comments[key])
    
    # Combine formatted fits columns and header into output binary fits table
    fitstable = fits.BinTableHDU.from_columns(columns, header=finalheader)
    
    LOG.info("{} sources in catalogue".format(len(fitstable)))
    #mem_use(label="at end of format_cat")
    return fitstable


# In[ ]:


# Based on BlackBOX function [write_fits] (for compatibility with Google Cloud system)
def save_fits(fitstable, outputname, rundir="", overwrite=True):
    
    """
    Function saves fits table in the directory specified in the string in
    [outputname]. For the Google cloud, writing the file directly to a bucket
    isn't possible. Instead, it's written to a folder [rundir] on a virtual
    machine before being moved over to the bucket.
    
    Parameters:
    -----------
    fitstable: binary table HDU object
        Fits table (data + header) that needs to be saved to file.
    outputname: string
        Full path to and name of the output file.
    rundir: string
        Only used in case of a Google cloud file system. Name of the folder on
        the VM where the fits table is temporarily written to.
    overwrite: boolean
        Determines whether fits table may be overwritten if it exists already.
    """
    #mem_use(label='save_fits at start')
    
    # Dealing with google cloud bucket?
    google_cloud = (outputname[0:5] == 'gs://')
    
    # Add time stamp of file creation to header
    header = fitstable.header
    header['DATEFILE'] = (Time.now().isot, 'UTC date of writing file')
    data = fitstable.data
    
    if not google_cloud:
        
        # Make dir for output file if it doesn't exist yet
        os.makedirs(os.path.dirname(outputname), exist_ok=True)
        
        # Write fits directly to the output [outputname]
        fits.writeto(outputname, data, header, overwrite=overwrite)
    
    else:
        
        # Write the tmp fits file
        fits_tmp = '{}{}'.format(rundir, outputname.split('/')[-1])
        LOG.info('writing tmp fits file {}'.format(fits_tmp))
        fits.writeto(fits_tmp, data, header, overwrite=overwrite)
        
        # Move fits_tmp to [dest_folder] in bucket
        dest_folder = os.path.dirname(outputname)
        copy_file(fits_tmp, dest_folder+'/', move=True)
    
    #mem_use(label='save_fits at end')
    return


# In[ ]:


# Based on BlackBOX function clean_tmp
def remove_tmp_folder(tmp_path):
    
    """
    Function that removes the specified folder and its contents.
    """
    
    if isdir(tmp_path):
        shutil.rmtree(tmp_path)
        LOG.info("Removing temporary folder: {}".format(tmp_path))
    else:
        LOG.warning("tmp folder {} does not exist".format(tmp_path))
    
    #mem_use(label="after removing temporary folder")
    return


# ## Run match2SSO using command line parameters

# In[ ]:


if __name__ == "__main__":
    
    PARSER = argparse.ArgumentParser(description="User parameters")
    
    PARSER.add_argument("--telescope", type=str, default="ML1",
                        help="Telescope name (ML1, BG2, BG3 or BG4); "
                        "default='ML1'")
    
    PARSER.add_argument("--mode", type=str, default="historic",
                        help="Day, night or historic mode of pipeline; "
                        "default='historic'")
    
    PARSER.add_argument("--catalog", type=str, default=None,
                        help="Only process this particular transient catalog. "
                        "Requires full path and requires mode to be 'historic'"
                        " or 'night'; default=None")
    
    PARSER.add_argument("--date", type=str, default=None,
                        help="Date to process (yyyymmdd, yyyy-mm-dd, yyyy/mm/d"
                        "d or yyyy.mm.dd). Mode must be 'historic' or 'day'; "
                        "default=None")
    
    PARSER.add_argument("--catlist", type=str, default=None,
                        help="Process all transient catalogs in the input list"
                        ". List entries require full path. Mode must be "
                        "'historic'; default=None")
    
    PARSER.add_argument("--logname", type=str, default=None,
                        help="Name of log file to save. Requires full path; "
                        "default of None will not create a log file")
    
    PARSER.add_argument("--redownload", type=str2bool, default=True,
                        help="Boolean to indicate whether the asteroid (and "
                        "comet) databases should be redownloaded. Alternatively"
                        " the most recently downloaded databases will be used.")
    
    PARSER.add_argument("--overwrite", type=str2bool, default=False,
                        help="Boolean to indicate whether files will be remade "
                        "and overwritten. Alternatively existing files will be "
                        "used.")
    
    ARGS = PARSER.parse_args()
    run_match2SSO(tel=ARGS.telescope, mode=ARGS.mode, cat2process=ARGS.catalog,
                  date2process=ARGS.date, list2process=ARGS.catlist,
                  logname=ARGS.logname, redownload=ARGS.redownload,
                  overwrite=ARGS.overwrite)

