import os

#==============================================================================
#Directory structure
#==============================================================================
runFolderBase = os.environ['DATAHOME']

runFolder = {}; inputFolder={}; databaseFolder={}; logFolder={}; submissionFolder={}
for tel in ['ML1', 'BG2', 'BG3', 'BG4']:
    runFolder[tel] = '{}/{}'.format(runFolderBase, tel)
    inputFolder[tel] = '{}/red/'.format(runFolder[tel])
    databaseFolder[tel] = '{}/tmp/match2SSO/'.format(runFolder[tel])
    logFolder[tel] = '{}/log/match2SSO/'.format(runFolder[tel])
    submissionFolder[tel] = '{}/mpc/'.format(runFolder[tel])

softwareFolder = "/Software/match2SSO/"

#==============================================================================
#Text file listing the software versions used
#==============================================================================
versionsFile = "/Software/versions.txt"


#==============================================================================
#Astcheck parameters
#==============================================================================
matchingRadius = 20 #matching radius in arcsec
limitingMagnitude = 25 #limiting V-magnitude

#Maximal number of asteroids that are returned as a match by astcheck. Set to a
#value large enough so that it won't restrict the output.
maximalNumberOfAsteroids = 1000


#==============================================================================
#Telescope parameters (see pytz.all_timezones within Python for possible time
#zones).
#==============================================================================
timeZoneTelescope = {"ML": "Africa/Johannesburg", "BG": "America/Santiago"}
FOV_width = 1.6544 # Size of the FOV in degrees
#For a square FOV, the FOV_width corresponds to the width and height of the FOV.
#For a circular FOV, the FOV_width is the diameter of the circle.

#==============================================================================
#Links to databases of known asteroids (MPC) and known comets (JPL)
#==============================================================================
URL_asteroidDatabase = "https://www.minorplanetcenter.net/iau/MPCORB/MPCORB.DAT"
URL_cometDatabase = "https://ssd.jpl.nasa.gov/dat/ELEMENTS.COMET"

#Maximal uncertainty parameter allowed for the asteroids that are used for the
#matching (see https://www.minorplanetcenter.net/iau/info/UValue.html). Allowed
#values for the maximum uncertainty parameter are 0 to 9 or None. If None, all 
#objects will be taken into account, including those with letter uncertainties.
maxUncertainty = 2


#==============================================================================
#JPL DE ephemeris file (describing planetary and lunar ephemerides) needed to
#integrate MPCORB to the observation epoch
#==============================================================================
JPL_ephemerisFile = "{}linux_m13000p17000.441".format(softwareFolder)


#==============================================================================
#Libraries
#==============================================================================
#CPP library to convert from __cplusplus macro to the language version string
CPPmacro2version = {"202002L":"C++20", "201703L":"C++17", "201402L":"C++14",
                    "201103L":"C++11", "199711L":"C++98"}


#==============================================================================
#Switches
#==============================================================================
redownload_databases = True
include_comets = False
keep_tmp = False
overwrite_files = False
time_functions = False

"""
Description of switches:
------------------------
redownload_databases: string
    Boolean indicating whether the asteroid and comet databases will need to be
    redownloaded when making the known objects database. Alternatively, the
    most recently downloaded version of the databases are used.

include_comets:
    Boolean indicating whether comets should be included in the known objects
    database. There have been issues with matching to comets with large orbital
    uncertainties in astcheck, so as long as this has not been solved, comet
    matching should be avoided. (Time of writing: 23-12-2021)

keep_tmp:
    Boolean indicating whether the temporary files made during the processing
    should be kept or removed at the end of the processing.

overwrite_files:
    Boolean indicating whether files are allowed to be overwritten. If False
    and the SSO catalogue and MPC submission files both already exist, the
    observation will be skipped.

time_functions:
    Boolean indicating whether functions need to be (wall-)timed.
"""
