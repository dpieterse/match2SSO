import os


# Links to databases of known asteroids (MPC) and known comets (JPL)
URL_asteroidDatabase = "https://www.minorplanetcenter.net/iau/MPCORB/MPCORB.DAT"
URL_cometDatabase = "https://ssd.jpl.nasa.gov/dat/ELEMENTS.COMET"


# CPP library to convert from __cplusplus macro to the language version string
CPPmacro2version = {"202002L":"C++20", "201703L":"C++17", "201402L":"C++14",
                    "201103L":"C++11", "199711L":"C++98"}

# Switches
"""
Description of switches:
------------------------
include_comets:
    Boolean indicating whether comets should be included in the known objects
    database. There have been issues with matching to comets with large orbital
    uncertainties in astcheck, so as long as this has not been solved, comet
    matching should be avoided. (Time of writing: 23-12-2021)
keep_tmp:
    Boolean indicating whether the temporary files made during the processing
    should be kept or removed at the end of the processing.
time_functions:
    Boolean indicating whether functions need to be (wall-)timed.
"""
include_comets = True
keep_tmp = False
time_functions = True


#==============================================================================
# All parameters below may (allowed, not required) be specified as dictionaries
# where the value depends on the used telescope
#==============================================================================


# Directory structure
"""
Beware: in the Google cloud, the tmpFolder should be on a VM, NOT in a bucket
"""
runFolderBase = {}; runFolder = {}; inputFolder={}; tmpFolder={}; logFolder={};
MPCreportFolder={}
for tel in ["ML1"]:
    runFolderBase["ML"] = "/idia/projects/meerlicht"
    runFolder[tel] = "{}/{}".format(runFolderBase["ML"], tel)
    inputFolder[tel] = "{}/red/".format(runFolder[tel])
    tmpFolder[tel] = "{}/tmp/match2SSO/".format(runFolder[tel])
    logFolder[tel] = "{}/log/match2SSO/".format(runFolder[tel])
    MPCreportFolder[tel] = "{}/mpc/".format(runFolder[tel])
for tel in ["BG2", "BG3", "BG4"]:
    runFolderBase["BG"] = "/home/sa_105685508700717199458"
    runFolder[tel] = "{}/Slurm/{}".format(runFolderBase["BG"], tel)
    inputFolder[tel] = "gs://blackgem-red/{}/".format(tel)
    tmpFolder[tel] = "{}/tmp/match2SSO/".format(runFolder[tel])
    logFolder[tel] = None
    MPCreportFolder[tel] = "{}/mpc/{}/".format(runFolderBase["BG"], tel)


# Text file listing the software versions used
versionsFile = "/Software/versions.txt"


# File listing the Minor Planet Center Observatory Codes
obsCodesFile = "/Software/match2SSO/ObsCodes.html"


# JPL DE ephemeris file (describing planetary and lunar ephemerides) needed to
# integrate MPCORB to the observation epoch
JPL_ephemerisFile = {}
for tel in ["ML", "BG"]:
    JPL_ephemerisFile[tel] = ("{}/CalFiles/linux_m13000p17000.441"
                              .format(runFolderBase[tel]))


# Relevant data columns from detection catalogue
colNumber = "NUMBER" # Source number, unique within the catalogue
colRA = "RA_PSF_D"   # [deg]
colDec = "DEC_PSF_D" # [deg]
colMag = "MAG_ZOGY"
colFlux = "FNU_ZOGY" # [micro Jy], only needed if magnitude column doesn't exist
colSNR = "SNR_ZOGY"  # Negative values are negative transients (to reject)


# Relevant header keywords from detection catalogue
keyDummy = "TDUMCAT"      # boolean (if True, the catalogue is empty)
keyDate = "DATE-OBS"      # isot format
keyMPCcode = "MPC-CODE"   # MPC observatory code
keyRACentre = "RA-CNTR"   # RA of the field center, [deg]
keyDecCentre = "DEC-CNTR" # Dec of the field center, [deg]
keyLimmag = "T-LMAG"      # transient limiting magnitude


# Astcheck parameters
matchingRadius = 20 #matching radius in arcsec
limitingMagnitude = 25 #limiting V-magnitude


# Maximal number of asteroids that are returned as a match by astcheck.
maximalNumberOfAsteroids = 1000


# Telescope parameters
"""
See pytz.all_timezones within Python for possible time zones.
For a square FOV, the FOV_width corresponds to the width and height of the
FOV. For a circular FOV, the FOV_width is the diameter of the circle.
"""
timeZoneTelescope = {"ML": "Africa/Johannesburg", "BG": "America/Santiago"}
FOV_width = 1.6544 # Size of the FOV in degrees
mpc_code = {"ML1": "L66", "BG": "X17"} # only used in day mode

# Maximal orbital uncertainty parameter
"""
Maximum U allowed for the asteroids that are used for the matching (see
https://www.minorplanetcenter.net/iau/info/UValue.html). Allowed values for
U are 0 to 9 or None. If None, all objects will be taken into account, including
those with letter uncertainties.
"""
maxUncertainty = 2


# Header for the MPC report (listed per MPC observatory code)
MPCreportHeader = "".join([
        "CON Radboud University, Houtlaan 4, 6525XZ, Nijmegen, The Netherlands\n",
        "CON [d.pieterse@astro.ru.nl]\n",
        "OBS P. J. Groot, S. L. D. Bloemen\n",
        "MEA D. L. A. Pieterse, P. M. Vreeswijk\n",
        "TEL 0.65-m reflector + CCD\n",
        "NET Gaia-DR3\n",
        "AC2 mpc-response@blackgem.org\n"])
