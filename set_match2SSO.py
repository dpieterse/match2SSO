import os
import numpy as np


proc_env = ""
#==============================================================================
# Get ML/BG processing environment (test/staging/production) from BlackBOX
# settings file. This is only needed for the definition of the paths in this
# settings file. If you're using a telescope that is NOT BlackGEM or MeerLICHT,
# comment out the next four lines:
import sys
sys.path.append("/Software/BlackBOX/Settings")
import set_blackbox as set_bb
proc_env = set_bb.proc_env
#==============================================================================


# Links to databases of known asteroids (MPC) and known comets (JPL)
URL_asteroidDatabase = "https://www.minorplanetcenter.net/iau/MPCORB/MPCORB.DAT"
URL_cometDatabase = "https://ssd.jpl.nasa.gov/dat/ELEMENTS.COMET"


# CPP library to convert from __cplusplus macro to the language version string
CPPmacro2version = {"202302L":"C++23", "202002L":"C++20", "201703L":"C++17",
                    "201402L":"C++14", "201103L":"C++11", "199711L":"C++98"}


# Dictionary to convert MPC uncertainty parameter to runoff in arcsec/decade
"""
By default, match2SSO only keeps solar system objects which have well-
determined orbits to prevent spurious matches. This is done by selecting
objects on their maximal orbital runoff (in arcsec), using the dictionary
below (see www.minorplanetcenter.net/iau/info/UValue.html). The U2maxRunoff
dictionary gives the RUNOFF in arcsec/decade. The runoff in arcsec is then
calculated as RUNOFF/10 * delta_t, where delta_t is the time since the last
observation of the object, in years. If you want to turn off this filtering,
set the dictionary values for each entry to zero.
"""
U2maxRunoff = {0:1.0, 1:4.4, 2:19.6, 3:86.5, 4:382, 5:1692, 6:7488, 7:33121,
               8:146502, 9:np.inf}


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
get_notified:
    Boolean indicating whether a notification should be sent in case of an error
    (both critical and non-critical) occuring in match2SSO.
"""
include_comets = True
keep_tmp = False
time_functions = True
get_notified = True


#==============================================================================
# All parameters below may (allowed, not required) be specified as dictionaries
# where the value depends on the used telescope. If dependent on the telescope,
# use the telescope name/abbreviation you specify with --telescope when you run
# match2SSO from the command line. Alternatively, you may use the alphabetic
# part of the telescope name/abbreviation in the dictionaries below. E.g. when
# specifying "BG2" or "BG" the value will be valid for the BlackGEM-2 telescope
# or all BlackGEM telescopes (BG2+BG3+BG4), respectively.
#==============================================================================

# Directory structure
"""
Beware: in the Google cloud, the tmpFolder should be on a VM, NOT in a bucket
"""

inputFolder={};
for tel in ["ML1"]:
    subdir = ""
    if proc_env == "test":
        subdir = "test-env/"
    elif proc_env == "staging":
        subdir = "staging-env/"
    inputFolder[tel] = ("/idia/projects/meerlicht/data/{}red/{}/"
                        .format(subdir, tel))

for tel in ["BG2", "BG3", "BG4"]:
    superbucket = ""
    if proc_env == "test":
        superbucket = "blackgem-test-env/"
    elif proc_env == "staging":
        superbucket = "blackgem-staging-env/"
    inputFolder[tel] = "gs://{}blackgem-red/{}/".format(superbucket, tel)


runFolder={}; tmpFolder={}; logFolder={}
runFolder["ML"] = "/idia/projects/meerlicht/"
runFolder["BG"] = "/home/sa_105685508700717199458/"
for tel in ["ML", "BG"]:
    tmpFolder[tel] = "{}RunMatch2SSO/tmp/".format(runFolder[tel])
    logFolder[tel] = "{}RunMatch2SSO/log/".format(runFolder[tel])


# If you do not care about specifying software versions in the headers of the
# match2SSO output catalogues, just use an empty string for versionsFile below.
# If you do want to include the lunar and jpl_eph versions in the headers, make
# a text file listing lunar and jpl_eph on different lines, where after the
# software name and an extra space, the version number is listed.
# E.g. "lunar a693", where a693 is the name of the last lunar commit.
versionsFile = "/Software/versions.txt"


# File listing the Minor Planet Center Observatory Codes
obsCodesFile = "/Software/match2SSO/ObsCodes.html"


# JPL DE ephemeris file (describing planetary and lunar ephemerides) needed to
# integrate MPCORB to the observation epoch
JPL_ephemerisFile = {}
for tel in ["ML", "BG"]:
    JPL_ephemerisFile[tel] = ("{}CalFiles/linux_m13000p17000.441"
                              .format(runFolder[tel]))


# Relevant data columns from detection catalogue
colNumber = "NUMBER" # Source number, unique within the catalogue
colRA = "RA_PSF_D"   # [deg]
colDec = "DEC_PSF_D" # [deg]
colMag = "MAG_ZOGY"  # Used in MPC reports, NOT used in the matching
colFlux = "FNU_ZOGY" # [micro Jy], only needed if magnitude column doesn't exist
colSNR = "SNR_ZOGY"  # Used to reject negative-S/N sources. Set to None to not use


# Relevant header keywords from detection catalogue
keyDate = "DATE-OBS"      # isot format
keyRACentre = "RA-CNTR"   # RA of the pointing center, [deg]
keyDecCentre = "DEC-CNTR" # Dec of the pointing center, [deg]


# Optional header keywords from detection catalogue (set to None to not use).
# keyDummy can speed up the processing of empty catalogues. keyLimmag is used
# to determine header keyword N-SSO for the output: the number of reasonably
# bright SSOs in the telescope's FOV (where V-magnitude <= lim. mag). Set to
# None to include all objects in the FOV.
keyDummy = "TDUMCAT" # boolean indicating if catalogue is empty
keyLimmag = "T-LMAG" # transient limiting magnitude


# Astcheck parameters
matchingRadius = 20    # matching radius in arcsec
limitingMagnitude = 25 # limiting V-magnitude


# Maximal number of asteroids that are returned as a match by astcheck. Set it
# high enough so that, if/when match2SSO makes prediction catalogues, the total
# number of solar system objects in the FOV is lower than this number.
# Otherwise the list of asteroids used for the prediction catalogues will be
# truncated.
maximalNumberOfAsteroids = 1000


# Telescope parameters
"""
Run "import pytz; print(pytz.all_timezones)" within Python for possible time
zones. If FOV_is_circle is True, we assume a circular field-of-view. Otherwise,
a square FOV is assumed, where the sides are aligned with RA & Dec. For a
circular FOV, the FOV_width is the diameter of the circle. For a square FOV,
the FOV_width corresponds to the width and height of the FOV.
"""
timeZoneTelescope = {"ML": "Africa/Johannesburg", "BG": "America/Santiago"}
FOV_is_circle = False
FOV_width = 1.6544 # Size of the FOV in degrees
mpc_code = {"ML": "L66", "BG": "X17"} # only used in day mode


# Header for the MPC report (listed per MPC observatory code). The report is
# not automatically sent out to the Minor Planet Center, so no need to worry.
MPCreportHeader = "".join([
        "CON Radboud University, Houtlaan 4, 6525XZ, Nijmegen, The Netherlands\n",
        "CON [d.pieterse@astro.ru.nl]\n",
        "OBS P. J. Groot, S. L. D. Bloemen\n",
        "MEA D. L. A. Pieterse, P. M. Vreeswijk\n",
        "TEL 0.65-m reflector + CCD\n",
        "NET Gaia-DR3\n",
        "AC2 mpc-response@blackgem.org\n"])


# Send notification to gcloud logging for errors instead of sending email
notify_in_gcloud = {"ML": False, "BG": False}

# Email settings
sender = "<danielle@blackgem.org>"
# comma-separated BlackGEM email addresses of recipients:
recipients = "danielle@blackgem.org"
reply_to = "danielle@blackgem.org"
smtp_server = "smtp-relay.gmail.com"
port = 465
use_SSL = True
