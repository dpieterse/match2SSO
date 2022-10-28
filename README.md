# match2SSO
Python code to match ***single*** telescope detections to ***known solar system objects*** based on their ***position***. To keep the matching accurate, only ***known objects with small orbital uncertainties*** are used in the matching. The maximum allowed uncertainty is set in _set_match2SSO.py_.

Asteroid matching is done by default. Comet matching is only done when ```include_comets = True``` (in the settings file). Currently, comet matching is not performed by default, as a selection based on orbital uncertainty is less straight forward for comets and Bill Gray's _integrat.cpp_ does not work on JPL's comet database.

_match2SSO_ was created for the [MeerLICHT](http://www.meerlicht.uct.ac.za/) & [BlackGEM](https://astro.ru.nl/blackgem/) telescopes. Originally written in Jupyter Notebook, it's also available as an executable Python script. It's compatible with MeerLICHT's image processing software (BlackBOX & ZOGY) version 1.0.0 and up.

_match2SSO_ makes grateful use of the [_lunar_](https://github.com/Bill-Gray/lunar) and [_jpl_eph_](https://github.com/Bill-Gray/jpl_eph) repositories that were written by Bill Gray under Project Pluto. The core of _match2SSO_ is [_astcheck_](https://www.projectpluto.com/astcheck.htm): a C++ script in the _lunar_ repository that matches detections to known solar system objects.


## Installation
- Install Python (code was tested on Python 3.8.10) and C++ 
- Clone Bill Gray's _lunar_ repository to a software folder (https://github.com/Bill-Gray/lunar) and build the library & executables as described in this repository's README file. Also build _integrat_ by running ```make integrat```.
- Add the lunar folder to PATH in your bash file
- Clone Bill Gray's _jpl_eph_ repository to the software folder (https://github.com/Bill-Gray/jpl_eph) and build the library & executables as described in this repository's README file.
- Clone the match2SSO repository.
- Download JPL's DE ephemeris file to the software folder. (ftp://ssd.jpl.nasa.gov/pub/eph/planets/Linux/)
- Download MPC's Observatory codes list to the software folder. (https://www.minorplanetcenter.net/iau/lists/ObsCodes.html)


## Input & output files
#### Input
_match2SSO_ runs on a FITS catalogue of detections. For MeerLICHT & BlackGEM, it is best run on the transient catalogues that were produced after difference imaging.

Other input files that the software needs are mentioned above under "Installation".
In addition, _match2SSO_ uses MPCORB.DAT (MPC's asteroid database) and COMET.ELEMENTS (JPL's comet database), but these are downloaded when running the script and hence do not need to be pre-downloaded.

#### Output
- Solar Sytem Object (SSO) catalogue containing the matches between detections and known solar system objects (__sso.fits_). SSO catalogue columns and header keywords are listed here: https://www.overleaf.com/read/zrhqwcbkfqns
- Predictions catalogue containing the asteroids that were predicted to be in the FOV during the observation (__sso_predict.fits_). Prediction catalogue columns and header keywords are listed in the overleaf document mentioned above.
- MPC submission file (__submit.txt_), to allow easy submission of the known object detections to the Minor Planet Center. Submission files are created, but not submitted automatically. 


## Running match2SSO
Run the _match2SSO.py_ script from the command line. Alternatively, the code is available in a Jupyter Notebook. It can be run in three modes:
- Historic mode: run on existing data. 
- Day mode: needs to be executed once before the start of an observing night, to allow the night mode to be run in real-time during that night. Allows speedy and parallelized processing in night mode.
- Night mode: run during the observing run (in real time) on a single detection catalogue. Needs the products made during the day mode.

More details on these modes are listed under "Code description" below.


### Main command line parameters:
- ```--mode``` Historic, night or day mode
- ```--catalog``` Name of detection catalogue to run the matching on
- ```--date``` Run _match2SSO_ run on all detection catalogues corresponding to this observing night.
- ```--catlist``` Run _match2SSO_ on all detection catalogues listed in this file.

Allowed combinations of the above-mentioned parameters are:
- Day mode
- Day mode + date
- Night mode + catalog
- Historic mode + catalog
- Historic mode + date
- Historic mode + catlist


### Other command line parameters:
- ```--telescope``` (Abbreviated) telescope name. Allows dictionaries in the settings file _set_match2SSO.py_ with different parameter values for different telescopes.
- ```--logname``` Name of the log file


### Multi-processing
Although multi-processing was not implemented within the code, efforts have been made to allow calling the code multiple times in parallel.
- **Night mode:** runs independently on a single catalogue. The steps in the code that make parallelization impossible have been moved to the day mode. Prepare for the night mode by running the day mode once before the start of an observing night. The night mode can be run on multiple catalogues of the same night in parallel without issues. 
- **Historic mode:** parallelization is only possible for data that was not taken on the same night. Multiple nights can be processed in parallel easily by running the historic mode on those nights individually.


## Code description
![Click here for a flow chart of match2SSO.](https://github.com/dpieterse/match2SSO/blob/master/match2SSO_flow.png?raw=true)

The steps match2SSO performs in the different modes are:
- **Day mode**
    1. Creates a run directory in preparation of the nightly processing
    2. Downloads asteroid and comet databases to the database folder.
    3. Integrates the asteroid database to midnight of the observation night
    4. Combines the comet and integrated asteroid databases into a SOF-formatted known objects database.
    5. Creates symbolic links in the run directory to the used databases and the observatory codes list.
    6. Runs astcheck on a fake detection in order to create the CHK files that astcheck will need for faster / parallel processing when running on observations. 
    7. Removes the fake detection in- & output (but not the CHK files!).
    (Products of steps 1-2 are saved to the databaseFolder, those of steps 3-5
    to the run directory.)
- **Night mode**
    1. Converts the detection catalogue into an MPC-formatted text file.
    2. Runs astcheck on the central coordinates of the observation, to make predictions on the known solar system objects in the FOV.
    3. Makes a prediction catalogue of the known solar system objects in the FOV during the observation.
    4. Runs astcheck on the MPC-formatted text file, to find matches between the detections and known solar system objects. 
    5. Makes an SSO catalogue containing the matches.
    6. Makes an MPC submission file of the matches.
- **Historic mode**<br/> Runs on a single detection catalogue (observation), a night of observations or a list of observations. The observations are grouped and processed per observing night. The historic mode:
    1. Creates a run directory per observation night
    2. [_Optional_] Downloads asteroid and comet databases to the database folder. (If this step is skipped, the most recently downloaded version of the database will be used.)
    3. Integrates the asteroid database to midnight of the observation night.
    4. Combines the comet and integrated asteroid databases into a SOF-formatted known objects database.
    5. Creates symbolic links in the run directory to the used databases and the observatory codes list.
    6. Run the matching per detection catalogue:<br/>
          a. Converts the detection catalogue into an MPC-formatted text file.
          <br/>b. Runs astcheck on the central coordinates of the observation, to make predictions on the known solar system objects in the FOV.
          <br/>c. Makes a prediction catalogue of the known solar system objects in the FOV during the observation.
          <br/>d. Runs astcheck on the MPC-formatted text file, to find matches between the detections and known solar system objects. 
          <br/>e. Makes an SSO catalogue containing the matches.
          <br/>f. Makes an MPC submission file of the matches.<br/>
    7. [_Optional_] Removes the run directory, including files in it (SOF-formatted known objects database, symbolic links, MPC-formatted detection file, astcheck output text file). Also removes the integrated asteroid database.
    
## License
Copyright 2022 Dani&euml;lle Pieterse

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
