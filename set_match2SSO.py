#===============================================================================
#Directory structure
#===============================================================================
inputDataFolder = "/media/danielle/LaCie/Danielle/bbversion_1_0_0/"
softwareFolder = ("/home/danielle/Documents/ESA_project/Observing_analysis/" + 
	          "Asteroid_linking/")
databaseFolder = ("/home/danielle/Documents/ESA_project/Observing_analysis/" + 
	          "Asteroid_linking/match2SSO/tmp/")
logFolder = ("/home/danielle/Documents/ESA_project/Observing_analysis/" + 
	     "Asteroid_linking/match2SSO/log/")
submissionFolder = ("/home/danielle/Documents/ESA_project/Observing_analysis/" +
                    "Asteroid_linking/match2SSO/MPCsubmissions/")


#===============================================================================
#Astcheck parameters
#===============================================================================
matchingRadius = 20 #matching radius in arcsec
limitingMagnitude = 25 #limiting V-magnitude

#Maximal number of asteroids that are returned as a match by astcheck. Set to a
#value large enough so that it won't restrict the output.
maximalNumberOfAsteroids = 1000


#===============================================================================
#Telescope parameter (see pytz.all_timezones within Python for possible time
#zones).
#===============================================================================
timeZoneTelescope = {"ML": "Africa/Johannesburg", "BG": "America/Santiago"}


#===============================================================================
#Links to databases of known asteroids (MPC) and known comets (JPL)
#===============================================================================
URL_asteroidDatabase = "https://www.minorplanetcenter.net/iau/MPCORB/MPCORB.DAT"
URL_cometDatabase = "https://ssd.jpl.nasa.gov/dat/ELEMENTS.COMET"

#Maximal uncertainty parameter allowed for the asteroids that are used for the
#matching (see https://www.minorplanetcenter.net/iau/info/UValue.html). Allowed
#values for the maximum uncertainty parameter are 0 to 9 or None. If None, all 
#objects will be taken into account - including those with letter uncertainties.
maxUncertainty = 2


#===============================================================================
#JPL DE ephemeris file (describing planetary and lunar ephemerides) needed to
#integrate MPCORB to the observation epoch
#===============================================================================
JPL_ephemerisFile = "{}linux_m13000p17000.441".format(softwareFolder)


#===============================================================================
#Libraries
#===============================================================================
#CPP library to convert from __cplusplus macro to the language version string
CPPmacro2version = {"202002L":"C++20", "201703L":"C++17", "201402L":"C++14",
                    "201103L":"C++11", "199711L":"C++98"}

