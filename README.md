# match2SSO
Python code to match telescope detections to known solar system bodies (in MPCORB).

Created for the MeerLICHT & BlackGEM telescopes. Originally written in Jupyter Notebook, but also available as an executable Python script. Compatible with MeerLICHT's image processing software (BlackBOX & ZOGY) version 1.0.0 and up.

_match2SSO_ makes grateful use of the _lunar_ and _jpl_eph_ repositories that were written by Bill Gray under Project Pluto. The core of _match2SSO_ is _astcheck_: a C++ script in the _lunar_ repository that matches detections to known solar system objects. More information on _astcheck_ can be found at: https://www.projectpluto.com/astcheck.htm

## Installation
- Install Python (code was tested on Python 3.8.10) and C++ 
- Clone Bill Gray's _lunar_ repository to a software folder (https://github.com/Bill-Gray/lunar)
- Clone Bill Gray's _jpl_eph_ repository to the software folder (https://github.com/Bill-Gray/jpl_eph)
- Clone the match2SSO repository.
- Download JPL's DE ephemeris file to the software folder. (ftp://ssd.jpl.nasa.gov/pub/eph/planets/Linux/)
- Download MPC's Observatory codes list to the software folder. (https://www.minorplanetcenter.net/iau/lists/ObsCodes.html)

## Running match2SSO
Run the executable _match2SSO.py_ script from the command line. Alternatively, the code is available in a Jupyter Notebook. It can be run in three modes:
- Historic mode: run on existing data. 
- Night mode: runs on a single catalogue. Can be run in real time. 
- Day mode: needs to be executed once before the start of an observing night, to allow the night mode to be run in real-time during that night. Allows speedy and parallelized processing in night mode.

### Main command line parameters:
- **--mode** Historic, night or day mode
- **--catalog** Name of detection catalogue to run the matching on
- **--date** Run _match2SSO_ run on all detection catalogues corresponding to this observing night.
- **--catlist** Run _match2SSO_ on all detection catalogues listed in this file.

Allowed combinations of the above-mentioned parameters are:
- Day mode
- Day mode + date
- Night mode + catalog
- Historic mode + catalog
- Historic mode + date
- Historic mode + catlist

### Other command line parameters:
- **--telescope** (Abbreviated) telescope name. Allows dictionaries in the settings file _set_match2SSO.py_ with different parameter values for different telescopes.
- **--newdatabases** Boolean to indicate if new known object databases need to be downloaded. If False, the last downloaded version of the database will be used (if it exists).
- **--includecomets** Boolean to indicate if comets should be included in the matching. Including comet matching currently leads to issues.
- **--logname** Name of the log file
- **--keep_tmp** Boolean to indicate if temporary files need to be kept.
- **--overwrite** Boolean to indicate whether existing files are allowed to be overwritten.
- **--timing** Boolean to indicate whether functions need to be wall-timed. Timing information is saved to the log.

### Multi-processing
Although multi-processing was not implemented within the code, efforts have been made to allow calling the code multiple times in parallel.
- Night mode: runs independently on a single catalogue. The steps in the code that make parallelization impossible have been moved to the day mode. Prepare for the night mode by running the day mode once before the start of an observing night. The night mode can be run on multiple catalogues of the same night in parallel without issues. 
- Day mode: parallelization is only possible for data that was not taken on the same night. Multiple nights can be processed in parallel easily by running the day mode on those nights individually.

## Code description
![Click here for a flow chart of match2SSO.](https://github.com/dpieterse/match2SSO/blob/master/match2SSO_flow.png?raw=true)
A description of what the code does in the different modes is given in the function _run_match2SSO_ within _match2SSO_.

### Input files
_match2SSO_ runs on a FITS catalogue of detections. For MeerLICHT & BlackGEM, it is best run on the transient catalogues that were produced after difference imaging.

Other input files that the software needs are mentioned above under "Installation".
In addition, _match2SSO_ uses MPCORB.DAT (MPC's asteroid database) and COMET.ELEMENTS (JPL's comet database), but these are downloaded when running the script and hence do not need to be pre-downloaded.

### Output files
- Solar Sytem Object (SSO) catalogue containing the matches (__sso.fits_). SSO catalogue columns and header keywords are listed here: https://www.overleaf.com/read/zrhqwcbkfqns
- MPC submission file (__submit.txt_), to allow easy submission of the known object detections to the Minor Planet Center. Submission files are created, but not submitted automatically. 
