# These functions have been pulled out of PyMPChecker_Brute.py

import ephem # PyEphem module
# pyfits now deprecated; use astropy instead -- dkg
# import pyfits # PyFITS
from astropy.io import fits
# astLib also a bit dusty, use astropy.wcs instead -- dkg
from astropy.wcs import WCS
import numpy as np
import math
from astropy.io.fits import getheader

# Function sep_sort 
def sep_sort(input_list):
    s_catalog = input_list[0]
    s_ra_radians = input_list[1]
    s_dec_radians = input_list[2]
    s_full_date = input_list[3]
    list_radius = input_list[4]
# Read in catalog list and return top of catalog list sorted in ascending 
# separation from given coordinates.
    s_date = ephem.date(s_full_date)
    # Define the date in PyEphem. 
    s_unsorted_bodies = []
    # Then, for each object in the short catalog provided, the orbital 
    # parameters are read into PyEphem and then the orbit computed for the date 
    # of observation. The resulting RA/DEC of the object is compared to the 
    # RA/DEC of the query coordinates to measure the angular separation 
    # (reported here in arcseconds).
    for body in s_catalog:
        test_obj = ephem.readdb(body)
        test_obj.compute(s_date)
        test_obj_ra = test_obj.a_ra
        test_obj_dec = test_obj.a_dec
        test_obj_separation = 206264.806*(float(ephem.separation((test_obj_ra, 
            test_obj_dec), (s_ra_radians, s_dec_radians))))
        # Test if the object is within the search radius. If it is, add it to 
        # the array of unsorted bodies.
        if (test_obj_separation < list_radius):
            s_unsorted_bodies.append((test_obj.name, 
                test_obj_separation, str(test_obj_ra), 
                str(test_obj_dec), str(test_obj.mag)))
    # The recombined list will be sorted all together after all the jobs are 
    # done.
    return s_unsorted_bodies

# Function split_seq
def split_seq(seq, size):
# Quick function to divide up a large list into multiple small lists, attempting 
# to keep them all the same size. PyMPChecker uses this to cut up the relatively
# large catalog into multiple smaller lists for assignment by the parallel
# processing job server. This is necessary to load-balance across the jobs and
# processors.
    newseq = []
    splitsize = 1.0/size*len(seq)
    for i in range(size):
        newseq.append(seq[int(round(i*splitsize)):int(round((i+1)*splitsize))])
    return newseq

def extract_fits_params(fits_file, imagetype='IMAGE'):
    '''
    Strip the date and position information from the header of the FITS image.
    returns 
    date, footprint
    '''

    try:
#        hdul = fits.open(sys.argv[1])
        hdr = getheader(fits_file, imagetype)
    except:
        raise SystemExit("File not found")

    hdr_date = hdr['DATE-OBS']

    w = WCS(hdr)
    f = w.calc_footprint()

    return hdr_date, f
