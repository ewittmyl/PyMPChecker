# PyMPChecker

This repository allows offline and online searching of the Minor Planet
database for GOTO images.

The `Checker` class provides the functionality, principally
via a cone search. Example usage is below:

```
import PyMPChecker as pympc

mpc = pympc.Checker()
mpc.cone_search(4.98112, 7.13445,'2018-08-06T03:27:35',10.)
Sorting catalog list of 782053 minor planets
into 156 jobs across 4 processors.
Got 1 matches

mpc.table
            name      dist          RA        Dec    mag
0  (10) Hygiea 0  4.177362  0:19:55.19  7:08:04.0  10.87
```

The offline version sorts a catalog of solar system bodies by 
ascending separation (in arcseconds) from a given sky coordinate (ra dec). 
Python 2 code by Christopher Klein 
(cklein@astro.berkeley.edu), updated (including modifications to work with Python 3)
by Duncan Galloway from 2018 Sep.

The online version queries the Minor Planet webpage, and parses
the results from that search. Use the flag `online=True` to 
switch to the online version.

You can also search for all bodies within a FITS image. In this case, the observation epoch and search radius are derived from the file:
```
mpc.image_search('image_file.fits')

mpc.table
```

### How do I get set up? ###

* Clone the repository and install using pip:
```
python3 -m pip install .
```
* Once installation is complete (and about once a month thereafter) run
`scripts/MP_Database_Update.csh` to download the latest orbital
parameters file from the MPC and parse it into the PyMPC format. The parsing is
done by `MPCORB2MonthlyCatalog.py` and the resulting datafile is stored in 
`Monthly_Orbit_Catalogsi`. `PyMPChecker_Brute.py` will look in 
Monthly_Orbit_Catalogs for the month corresponding to the queried observation 
time, and an error will result if the necessary datafile is not available.
* For querying on GOTO images, you'll also need to install the [`gotocat`](https://github.com/GOTO-OBS/gotocat)
repository

Note on uncertainty parameter (provided as the last character
of the name string, in the offline search):
* 0-2 is good (on scale with coordinate uncertainty)
* 3-4 is mediocre (3-sigma semi-major axis is dozens of arcseconds)
* 5+ is bad (3-sigma semi-major axis is hundreds of arcseconds)
* ? indicates that none was given
* E indicates that the orbital eccentricity was assumed
* For one-opposition orbits D indicates that a double (or multiple)
    designation is involved and F indicates that an e-assumed
    double (or multiple) designation is involved.

I use `A` if the orbital parameters were generated from archival data.
This indicates that no uncertainty parameter was provided by the MPC.

Exact uncertainty information can be obtained (if available) by running the 
object and date through the MPES: 
  http://www.cfa.harvard.edu/iau/MPEph/MPEph.html

Software dependencies:
*  numpy
*  astropy
*  PyEphem (http://rhodesmill.org/pyephem/)

Deprecated dependencies:
*  Python 2.4
*  pyfits
*  astLib (http://astlib.sourceforge.net/?q=install0.3)
    
Many thanks to Brandon Craig Rhodes (http://rhodesmill.org/brandon) for writing 
and distributing PyEphem and for assistance in its implementation.

### Contribution guidelines ###

* Please contribute!

### Who do I talk to? ###

* Duncan.Galloway@monash.edu

