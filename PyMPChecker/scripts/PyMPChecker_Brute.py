"""
Python program PyMPChecker_Brute.py by Christopher Klein 
(cklein@astro.berkeley.edu). Program sorts a catalog of solar system bodies by 
ascending separation (in arcseconds) from a given sky coordinate (ra dec). 
The sky coordinates and observation time is given in the input arguments.

Updated by Duncan Galloway (duncan.galloway@monash.edu) from 2018 Sep

Note on uncertainty parameter:
* 0-2 is good (on scale with coordinate uncertainty)
* 3-4 is mediocre (3-sigma semi-major axis is dozens of arcseconds)
* 5+ is bad (3-sigma semi-major axis is hundreds of arcseconds)
* ? indicates that none was given
* E indicates that the orbital eccentricity was assumed
* For one-opposition orbits D indicates that a double (or multiple)
    designation is involved and F indicates that an e-assumed
    double (or multiple) designation is involved.

I use "A" if the orbital parameters were generated from archival data.
This indicates that no uncertainty parameter was provided by the MPC.

Exact uncertainty information can be obtained (if available) by running the 
object and date through the MPES: 
  http://www.cfa.harvard.edu/iau/MPEph/MPEph.html

Software dependencies:
    Python 2.4
    numpy
    astropy
    PyEphem (http://rhodesmill.org/pyephem/)
Deprecated dependencies:
    pyfits
    astLib (http://astlib.sourceforge.net/?q=install0.3)
Many thanks to Brandon Craig Rhodes (http://rhodesmill.org/brandon) for writing 
and distributing PyEphem and for assistance in its implementation.
"""

import sys
import operator
import numpy as np
import ephem # PyEphem module
# import pp # Parallel Python module
from astLib import astCoords
from multiprocessing import Pool, cpu_count
import pkg_resources

import PyMPChecker as pympc

# some definitions

mag_upper_limit=6.2

from time import time
t1 = time()     # Record time at the beginning of program execution.

# **** BEGIN MAIN PROGRAM HERE ****

if (len(sys.argv) != 14) and (len(sys.argv) != 2):
    raise SystemExit("Input Usage: python PyMPChecker_Brute.py ra_h ra_m " +  
    "ra_s dec_d dec_m dec_s date_year date_month date_day date_hour " +  
    "date_minute date_second sorted_list_radius\nor\npython PyMPChecker.py " +
    "fits_image")

# Example execution call: 
#   python PyMPChecker_Brute.py 6 7 8.30 20 18 42.0 2008 7 27 7 8 9 500
#   python PyMPChecker_Brute.py fake_0_00.fits

neg_dec_flag = False

if len(sys.argv) == 2:
#    hdr = pyfits.getheader(sys.argv[1])
#    hdr_date = hdr["UTSHUT"]

# Wrapped up this part into a routine so we can call it independently of
# this script

    year, month, day, hour, min, sec, ra_degrees, dec_degrees, sorted_list_radius = pympc.extract_fits_params(sys.argv[1])

if len(sys.argv) == 14:
# Parse the 13 input arguements as the observation coordinates, time, and 
# requested search radius.
    ra_h = int(sys.argv[1])
    ra_m = int(sys.argv[2])
    ra_s = float(sys.argv[3])
    dec_d = int(sys.argv[4])
    dec_m = int(sys.argv[5])
    dec_s = float(sys.argv[6])
    year = sys.argv[7]
    month = sys.argv[8]
    day = sys.argv[9]
    hour = sys.argv[10]
    min = sys.argv[11]
    sec = sys.argv[12]
    sorted_list_radius = int(sys.argv[13])
    # Run some checks on the input ra and dec values.
    if ra_h < 0: raise SystemExit("Right Ascension Hours cannot be negative.")
    if ra_h > 23: raise SystemExit("Right Ascension Hours " + 
    "cannot be greater than 23.")
    if ra_m < 0: raise SystemExit("Right Ascension Minutes cannot be negative.")
    if ra_m > 59: raise SystemExit("Right Ascension Minutes " + 
    "cannot be greater than 59.")
    if ra_s < 0: raise SystemExit("Right Ascension Seconds cannot be negative.")
    if ra_s >= 60: raise SystemExit("Right Ascension Seconds " + 
    "cannot be greater than or equal to 60.")
    if dec_d <= -90: raise SystemExit("Declination Degrees cannot be " + 
    "less than or equal to -90.")
    if dec_d >= 90: raise SystemExit("Declination Degrees " + 
    "cannot be greater than or equal to 90.")
    if dec_m < 0: raise SystemExit("Declination Minutes cannot be negative.")
    if dec_m > 59: raise SystemExit("Declination Minutes " + 
    "cannot be greater than 59.")
    if dec_s < 0: raise SystemExit("Declination Seconds cannot be negative.")
    if dec_s >= 60: raise SystemExit("Declination Seconds " + 
    "cannot be greater than or equal to 60.")
    if sorted_list_radius <= 0: raise SystemExit("Sorted List Radius cannot be "
    + "less than or equal to 0.")
    # Deal with negative declination values. This makes the later converstion to 
    # decimal degrees and radians correct.
    if sys.argv[4][:1] == "-":
        dec_m = -dec_m
        dec_s = -dec_s
        neg_dec_flag = True

# Recreate the full date from the parsed information. When PyEphem reads in the
# date, it correctly interprets negative values, month values outside [1:12], 
# hour values outside [0:23], and min/sec values outside [0:59].
full_date = year + "/" + month + "/" + day + " " + hour + ":" + min + ":" + sec
print full_date, ra_degrees, dec_degrees, sorted_list_radius

# Figure out which monthly orbit catalog to read in based on the month of the 
# query observation time.
if month == "1": month2 = "01"
if month == "2": month2 = "02"
if month == "3": month2 = "03"
if month == "4": month2 = "04"
if month == "5": month2 = "05"
if month == "6": month2 = "06"
if month == "7": month2 = "07"
if month == "8": month2 = "08"
if month == "9": month2 = "09"
if month == "10": month2 = "10"
if month == "11": month2 = "11"
if month == "12": month2 = "12"
observation_month = year + "_" + month2

start_read_t = time()
# Turn appropriate orbit datafile into catalog_list.
catalog_list = []
try:
#    f_catalog = file("Monthly_Orbit_Catalogs/" + observation_month + 
    f_catalog = file(pkg_resources.resource_filename('PyMPChecker','') 
        +'/Monthly_Orbit_Catalogs/'+ observation_month + "_ORB.DAT", "r")
except (IOError):
    raise SystemExit("You have requested a date for which no orbital " + 
    "parameters exist.")
for line in f_catalog:
    catalog_list.append(line.rstrip())
f_catalog.close()
end_read_t = time()

if len(sys.argv) == 14:
    # Print out summary of the program's input parameters.
    if (neg_dec_flag):
        print "Begin search of Minor Planet Catalog for known bodies near the sky", 
        print "coordinates (ra, dec) %d:%d:%.2f, %s%d:%d:%.2f at time %s." % (ra_h, 
        ra_m, ra_s, "-", abs(dec_d), -dec_m, -dec_s, ephem.date(full_date))
    else:
        print "Begin search of Minor Planet Catalog for known bodies near the sky", 
        print "coordinates (ra, dec) %d:%d:%.2f, %d:%d:%.2f at time %s." % (ra_h, 
        ra_m, ra_s, dec_d, dec_m, dec_s, ephem.date(full_date))

if len(sys.argv) == 2:
    print ("Begin search of Minor Planet Catalog for known bodies within " + 
    "the image %s." % sys.argv[1])

# Calculate the observation coordinates in degrees and radians.
if len(sys.argv) == 14:
    ra_degrees = 15*(ra_h + (ra_m + ra_s/60)/60)
    dec_degrees = (dec_d + (dec_m + dec_s/60)/60)
ra_radians = 0.0174533*ra_degrees
dec_radians = 0.0174533*dec_degrees

start_calculate_t = time()

num_jobs = (len(catalog_list)/20000)*(cpu_count())
print "Sorting catalog list of", len(catalog_list), "minor planets",
print "into", num_jobs, "jobs across", cpu_count(), "processors."

sub_catalog_list = pympc.split_seq(catalog_list, num_jobs)

processing_inputs = []
for sublist in sub_catalog_list:
    processing_inputs.append([sublist, ra_radians, dec_radians, full_date, sorted_list_radius])

# This is the bit that actually does the work

p = Pool(cpu_count())
result = p.map_async(pympc.sep_sort, processing_inputs)
close_bodies = result.get()

recomb_close_bodies = []
for job in range(num_jobs):
    for list in range(len(close_bodies[job])):
        recomb_close_bodies.append(close_bodies[job][list])

final_sorted_bodies = sorted(recomb_close_bodies, 
key=operator.itemgetter(1))

end_calculate_t = time()

# Print the list of close minor planets.
if len(final_sorted_bodies) > 0:
# if len(sys.argv) == 14:
    print ("All minor planets found within %i arcseconds of query point:" % 
    sorted_list_radius)
    print ("Object\t\tSeparation (\")\tRA\t\tDEC\t\tMag (V)" + 
    "\t\tUncertainty Parameter")
    output_file = file("nearby_minor_planets.txt", "w")
    output_file.write("Object\t\tSeparation (\")\tRA\t\tDEC\t\tMag (V)" + 
    "\t\tUncertainty Parameter\n")
    for body in final_sorted_bodies:
        output_sep = body[1]
        output_name = str(body[0])[:-2]
        output_ra = body[2]
        output_dec = body[3]
        output_mag = body[4]
        if float(output_mag) < 6.2: output_mag = "ERR"
        output_unc = (body[0])[-1]
        print ("%s\t\t%.2f\t\t%s\t%s\t%s\t\t%s" % (output_name, 
        output_sep, output_ra, output_dec, output_mag, output_unc))
        output_file.write("'%s'\t\t%.2f\t\t%s\t%s\t%s\t\t%s\n" % (output_name, 
        output_sep, output_ra, output_dec, output_mag, output_unc))
    output_file.close()

# If the FITS image was provided, then the program will write a region
# file to show the positions of the identified minor planets. 
# Furthermore, if a minor planet is more than 2200 arc seconds away, it
# is ignored because it is definitely outside the CCD image area.
if len(sys.argv) == 2 & len(final_sorted_bodies) > 0:
    # Make an SAOImage DS9 region file to write to.
    region_file_output = file(sys.argv[1][:-5] + ".reg", "w")
    region_file_output.write("global color=green font=\"helvetica 10 " + 
    "normal\" select=1 highlite=0 edit=0 move=0 delete=1 include=1 fixed=0 " +  
    "source\nfk5\n")
    # Print out the final results of the program and write the region file.
    print ("Object\t\tSeparation (\")\tRA\t\tDEC\t\tMag (V)" + 
    "\t\tUncertainty Parameter")
# don't know how to get the WCS rotation anymore. According to M. Dyer,
# "there's no easy way to get the rotation angle through astropy,
# annoyingly. The best method I could find was based on
# https://stackoverflow.com/questions/17332853/how-to-find-the-angle-between-north-and-horizon-for-given-altaz-coords-using-pye
# . You can get the matrix by doing `w.wcs.cd` "
#    rotation_degrees = WCS.getRotationDeg()
    rotation_degrees = None
    if (rotation_degrees < 5):
    # If the image has a position angle less than 5 degrees, it's close enough
    # to permit defining a box to exclude minor planets which fall outside the 
    # CCD area but which are within 2200 arcseconds of the center.
# This has already been read in
#        half_size = WCS.getHalfSizeDeg()
        half_size = half_size_array
        half_ra = half_size[0]
        half_dec = half_size[1]
        # Define bounding box of image. I assume the image is aligned with the 
        # RA/DEC coordinate system, i.e., "left" is increasing RA and "up" is 
        # increasing DEC. This will break down if the image is rotated on the sky,
        # that is, if the image has a non-zero position angle. I pad the region 
        # by 20% to be safe.
        left_ra_degrees = ra_degrees + half_ra*1.05
        right_ra_degrees = ra_degrees - half_ra*1.05
        top_dec_degrees = dec_degrees + half_dec*1.05
        bottom_dec_degrees = dec_degrees - half_dec*1.05
    
    for body in final_sorted_bodies:
        output_sep = body[1]
        output_name = str(body[0])[:-2]
        output_ra = body[2]
        output_ra_degrees = astCoords.hms2decimal(output_ra, ":")
        output_dec = body[3]
        output_dec_degrees = astCoords.dms2decimal(output_dec, ":")
        output_mag = body[4]
        if float(output_mag) < mag_upper_limit: output_mag = "ERR"
        output_unc = (body[0])[-1]
        if (((rotation_degrees < 5) 
            and (output_ra_degrees > right_ra_degrees)
            and (output_ra_degrees < left_ra_degrees)
            and (output_dec_degrees > bottom_dec_degrees)
            and (output_dec_degrees < top_dec_degrees))
            or (rotation_degrees >= 5)):
            print ("%s\t\t%.2f\t\t%s\t%s\t%s\t\t%s" % (output_name, 
            output_sep, output_ra, output_dec, output_mag, output_unc))
            
            # The region file is loosely color-and-radius coded to reflect the 
            # uncertainty parameter associated with the orbit of the minor planet.
            # See the note on uncertainty parameter below for further explanation.
            # In a nutshell: Small, green circle => minor planet ought to be in 
            #                    circle.
            #                Medium, yellow circle => minor planet might be in 
            #                    circle, but circle represents middle of uncertainty
            #                    region.
            #                Large, red circle => minor planet middle is at circle,
            #                    but MPC gave very large uncertainty.
            #                Large, purple circle => minor planet middle is at
            #                    circle, but MPC gave no uncertainty, an E, D, or F.
            if output_unc == "0" or output_unc == "1" or output_unc == "2":
                region_file_output.write("circle(" + output_ra + "," + 
                output_dec + 
                ",10\")	# color=green font=\"helvetica 18 normal\" text={" + 
                output_name + ", mag: " + output_mag + ", U=" + 
                output_unc + "}\n")
            elif output_unc == "3" or output_unc == "4":
                region_file_output.write("circle(" + output_ra + "," + 
                output_dec + 
                ",20\")	# color=yellow font=\"helvetica 18 normal\" text={" + 
                output_name + ", mag: " + output_mag + ", U=" + 
                output_unc + "}\n")
            elif (output_unc == "5" or output_unc == "6" or output_unc == "7" or 
            output_unc == "8" or output_unc == "9"):
                region_file_output.write("circle(" + output_ra + "," + 
                output_dec + 
                ",30\")	# color=red font=\"helvetica 18 normal\" text={" + 
                output_name + ", mag: " + output_mag + ", U=" + 
                output_unc + "}\n")
            elif (output_unc == "A"):
                region_file_output.write("circle(" + output_ra + "," + 
                output_dec + 
                ",12.5\")	# color=blue font=\"helvetica 18 normal\" text={" + 
                output_name + ", mag: " + output_mag + ", U=" + 
                output_unc + "}\n")
            else:
                region_file_output.write("circle(" + output_ra + "," + 
                output_dec + 
                ",30\")	# color=purple font=\"helvetica 18 normal\" text={" + 
                output_name + ", mag: " + output_mag + ", U=" + 
                output_unc + "}\n")
    # Close the region file.
    region_file_output.close()














total_time = time() - t1
read_time = end_read_t - start_read_t
calculate_time = end_calculate_t - start_calculate_t
# Print the total program execution time.
print """Total time elapsed: %.2f seconds. 
Catalog read time: %.2f seconds (%.2f). 
Orbital calculation time: %.2f seconds (%.2f).""" % (total_time, 
    read_time, 100*read_time/total_time, calculate_time, 
    100*calculate_time/total_time)


