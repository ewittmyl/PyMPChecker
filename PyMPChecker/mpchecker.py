# This is the class created to replace the functionality of PyMPChecker_Brute.py

from __future__ import print_function

import sys
import os
import re
import operator
import ephem # PyEphem module
from multiprocessing import Pool, cpu_count
import pkg_resources
import time
import pandas as pd
import psycopg2
from datetime import datetime
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np

from .utils import *
from .online_search import *

# Following package is to enable database queries

try:
    import gotocat as gc
except:
    print ('PyMPChecker: ** WARNING ** gotocat is unavailable, database queries will fail')

class Checker:

    def __init__(self, interactive=True, verbose=False):
        '''
        Create the object and determine the list of available data files
        '''

        self.interactive = interactive
        self.verbose = verbose

        self.dir = pkg_resources.resource_filename('PyMPChecker', '/Monthly_Orbit_Catalogs/')
        all_files = os.listdir(self.dir)
        patt = '[0-9]{4}_[0-9]{2}'
        self.epoch = [re.search(patt,x).group(0) for x in all_files if re.search(patt,x)]

        # Check here if the epoch is too different from the current one

        self.month = datetime.now().month
        self.year = datetime.now().year
        self.epoch_num = np.fromiter(map(lambda x: float(x[0:4])+(float(x[5:7])-1)/12.,
                                         self.epoch), dtype=np.float)
        self.yelapsed = min(abs(self.epoch_num-(self.year+(self.month-1)/12.)))
        if self.yelapsed > 0.5/12.:
            print ('PyMPChecker: ** WARNING ** your ephemeris database is >1 month old, you should update')

        # Database parameters (don't connect until you need to)
        # These are also only for the phase 1 information; need to update the database
        # interaction

        self.conn_string = "host='goto-observatory.org' dbname='photometry' user='goto' password='H4iidNHr'"
        self.conn = None

        self.observation_month_read = None

        self.search_run = False
        self.offline_search_run = False
        self.online_search_run = False

    def cone_search(self, ra_degrees, dec_degrees, obs_time, radius=1, online=False,
                    _debug=False):
        '''
        Fundamental search for the package. Will run a search centred around the provided
        position, at the provided epoch, and return the results as a pandas table
        :param ra_degrees:
        :param dec_degrees:
        :param obs_time:
        :return:
        Example shown below; this search should return (325395) 2009 CQ5, although
        the offline searched failed on the day
        mpc.cone_search(227.68176, 33.36493, '2019-09-17T20:44:00',radius=10.)
        '''

        self.online_search_run = False
        self.offline_search_run = False
        self.search_run = False

        if online:
            # This option uses a lightly edited version of Krzysztof's online
            # checker

            start_calculate_t = time.time()

            # date string decimals on seconds is limited to 6
            # search radius is in arcminutes, and smallest allowed value is 1

            result = find_asteroid(None, ra_degrees, dec_degrees, datestr=obs_time[:26],
                                   r=max([radius/60,1]),debug=_debug)
            if result is not None:

                res_table = pd.DataFrame(result,
                    columns=['name', 'dist', 'RA_deg', 'Dec_deg', 'mag', 'RA_off', 'Dec_off', 'comment'])
                if radius < 60.:
                    res_table = res_table.drop(res_table[res_table['dist'] > radius].index)
                self.n_results = len(res_table)
                if self.n_results > 0:
                    self.table = res_table

            else:
                self.n_results = 0

            end_calculate_t = time.time()
            self.online_search_run = True

        else:

            # Here we do the offline search
            # Position assumed to be in degrees

            ra_radians = 0.0174533 * ra_degrees
            dec_radians = 0.0174533 * dec_degrees

            # Time should be a stamp like '2018-08-06T03:27:05', as read in from the file header; not
            # sure the best way to check that it matches that format

            year = obs_time[0:4]
            month = obs_time[5:7]
            day = obs_time[8:10]
            hour = obs_time[11:13]
            min = obs_time[14:16]
            sec = obs_time[17:19]
            full_date = year + "/" + month + "/" + day + " " + hour + ":" + min + ":" + sec

            observation_month = year + "_" + month

    # You don't need to read in the file for this specific month; the filename just tells you when
    # the file was created

            yelapsed = abs(self.epoch_num-(float(year)+(float(month)-1)/12.))
            inearest = np.argmin(yelapsed)
            if ((not self.observation_month_read) |
                (self.epoch[inearest] != self.observation_month_read)):

    # Only read in the catalog list if we haven't already
    # pick the latest one available

                # self.observation_month_read = sorted(self.epoch)[-1]
                self.observation_month_read = self.epoch[inearest]
                if yelapsed[inearest] > 0.5/12.:
                    print ('''
PyMPChecker: ** WARNING ** nearest ephemeris file is >= 1 month away from target date
                           accuracy may be compromised''')

                start_read_t = time.time()
    # Turn appropriate orbit datafile into catalog_list.
                self.catalog_list = []
    #            f_catalog = open(self.dir+observation_month + "_ORB.DAT", "r")
                f_catalog = open(self.dir+self.observation_month_read + "_ORB.DAT", "r")
                for line in f_catalog:
                    self.catalog_list.append(line.rstrip())
                f_catalog.close()
                end_read_t = time.time()

                self.read_elapsed = end_read_t-start_read_t

            start_calculate_t = time.time()

            num_jobs = (len(self.catalog_list) // 20000) * (cpu_count())
            if self.verbose:
                print ("Sorting catalog list of", len(self.catalog_list), "minor planets",)
                print ("into", num_jobs, "jobs across", cpu_count(), "processors.")

            sub_catalog_list = split_seq(self.catalog_list, num_jobs)

            processing_inputs = []
            for sublist in sub_catalog_list:
                processing_inputs.append([sublist, ra_radians, dec_radians, full_date, radius])

            # This is the bit that actually does the work

            p = Pool(cpu_count())
            result = p.map_async(sep_sort, processing_inputs)
            close_bodies = result.get()

            # Added to fix issue
            # https://github.com/GOTO-OBS/PyMPChecker/issues/5

            p.close()
            p.join()

            recomb_close_bodies = []
            for job in range(num_jobs):
                for list in range(len(close_bodies[job])):
                    recomb_close_bodies.append(close_bodies[job][list])

            # The result here is a list of tuples, as the example below:
            # ('2014 WN13 0', 793.2264337773519, '0:19:25.96', '7:31:27.2', '20.9')
            # giving the name (with uncertainty character appended), distance from the provided position
            # (arcsec), RA, Dec, and magnitude

            self.final_sorted_bodies = sorted(recomb_close_bodies,
                                              key=operator.itemgetter(1))

            self.n_results = len(self.final_sorted_bodies)

            if self.n_results > 0:
                self.table = pd.DataFrame(self.final_sorted_bodies,
                                     columns=['name','dist','RA_str','Dec_str','mag'])

    # Convert the table into more useful data formats

                self.table['accuracy'] = [x[-1] for x in self.table.name]
                self.table.name = [x[:-2] for x in self.table.name]
                self.table.dist = pd.to_numeric(self.table.dist)
                self.table.mag = pd.to_numeric(self.table.mag)
                c = SkyCoord(ra=self.table.RA_str, dec=self.table.Dec_str,
                             unit=(u.hourangle, u.deg))
                self.table['RA_deg'] = c.ra.deg
                self.table['Dec_deg'] = c.dec.deg

            end_calculate_t = time.time()

            self.offline_search_run = True

        self.search_run = self.offline_search_run or self.online_search_run
        self.search_elapsed = end_calculate_t - start_calculate_t

        if self.verbose:
            if self.n_results > 0:
                print('Got {} matches'.format(self.n_results))
            else:
                print ('Got no matches in given region')


    def image_db_search(self, image_id, phase=None, online=False):
        '''
        This version will search over an image in the database, based on the footprint
        :param image_id:
        :return:
        '''

        if phase is None:
            print ("PyMPChecker: ** WARNING ** no phase supplied, assuming most recent data")
            self.phase = 4
        else:
            self.phase = phase

        g = gc.GOTOdb(phase=self.phase)

        # fov = g.query('SELECT fov, date FROM images WHERE id = %s', image_id)
        fov = g.query('SELECT fov, '+g.image_date+' FROM '+g.imagedb+' WHERE id = %s', image_id)

        self.footprint = np.array([list(i) for i in eval(fov['fov'].values[0])])
        # self.obs_date = str(fov['date'].values[0])
        self.obs_date = str(fov[g.image_date].values[0])

    #     if not self.conn:
    #         self.conn = psycopg2.connect(conn_string)
    #         self.cursor = self.conn.cursor()
    #
    #     self.comm="""SELECT id,filename,jd,date,target FROM images
    # WHERE instrument='UT4' AND filter='L' AND filename LIKE '%median.fits'
    # AND target='"""+sel_field+"'"

        self.footprint_search(online=online)


    def image_search(self, filename, imagetype='IMAGE', online=False):
        '''
        Check for minor planets within a FITS image
        '''

        try:
#            year, month, day, hour, min, sec, ra_degrees, dec_degrees, sorted_list_radius = extract_fits_params(filename)
             obs_date, footprint = extract_fits_params(filename, imagetype=imagetype)
        except:
            print ("** ERROR ** can't read information from file {}".format(filename))
            return

        self.obs_date = obs_date
        # print (obs_date)
        self.footprint = footprint

        self.footprint_search(online=online)


    def footprint_search(self, online=False):
        '''
        This routine is used by both the image search routines, since they both rely
        on searching within a GOTO footprint (but get the footprint differently)
        :return:
        '''

        if self.verbose:
            print ("Running footprint search over footprint\n{}".format(self.footprint))

        if ((self.footprint[1,0] > self.footprint[0,0]) |
            (self.footprint[2,0] > self.footprint[3,0])):

            # The FoV spans the RA = 0 line here, so we need to be a bit more careful
            # about how we calculate the centre
            # Example is for image #470457 or 470498 (use this for image_db_search)

            print ("Special here for crossing RA=0.0")
            center = [np.mean(self.footprint[:,0]+[360.0,0.,0.,360.0]),
                      np.mean(self.footprint[:,1])]
            half_size_array = [max(abs(self.footprint[:,0]+[360.0,0.0,0.0,360.0]
                                       -center[0]))*np.cos(np.deg2rad(center[1])),
                               max(abs(self.footprint[:,1]-center[1]))]#*3600.
            center[0] = center[0] % 360.
        else:
            center = [np.mean(self.footprint[:,0]),np.mean(self.footprint[:,1])]
            half_size_array = [max(abs(self.footprint[:,0]-center[0]))*np.cos(np.deg2rad(center[1])),
                               max(abs(self.footprint[:,1]-center[1]))]#*3600.
    #    print "half_size_array = ",half_size_array

        if self.verbose:
            print ("Translates to center {}\n and half_size_array {}".format(center, half_size_array))

        ra_degrees = center[0]
        dec_degrees = center[1]
    #    half_size_array = WCS.getHalfSizeDeg()

        sorted_list_radius = math.sqrt(half_size_array[0]*half_size_array[0] + half_size_array[1]*half_size_array[1])
    # provided value needs to be in arcsec
        sorted_list_radius *= 3600.

        self.cone_search(ra_degrees, dec_degrees, self.obs_date, sorted_list_radius,
                         online=online)

# Need to "mask" here to avoid the bits of the search cone that fall outside the field

        if self.n_results > 0:
            mask = ((self.table.RA_deg >= min(self.footprint[:,0])) &
                   (self.table.RA_deg <= max(self.footprint[:,0])) &
                   (self.table.Dec_deg >= min(self.footprint[:,1])) &
                   (self.table.Dec_deg <= max(self.footprint[:,1])))

            self.table = self.table[mask]
            self.n_results = len(self.table)
            if self.verbose:
                print ("... of which {} fall within the image footprint".format(self.n_results))

            if self.n_results == 0:
                # Don't return an empty table! This could cause confusion
                del self.table

    def update_database(self, year_month=None, update_db_file=False):
        '''
        Reread the source database, if necessary, and generate this month's file
        This method replaces MP_Database_Update.csh and MPCORB2MonthlyCatalog.py
        '''

        self.db_file = 'MPCORB.DAT'
        self.db_file_and_path = self.dir+self.db_file
        db_file_present = os.path.isfile(self.db_file_and_path)

        if (not db_file_present) | update_db_file:

# Would use the wget package but it doesn't allow the --no-check-certificate option

            if self.verbose:
                print ("Downloading a new copy of the database file, please wait...")
            try:
                os.system('wget --no-check-certificate https://www.minorplanetcenter.org/iau/MPCORB/'
                      +self.db_file+' -O '+self.db_file_and_path)
            except:
                print ("** ERROR ** can't download database file, aborting")
                return

            if self.verbose:
                print ("...done")

        if not year_month:

# Epoch is not specified, so carry out the calculation for the current year & month
# This command returns a different format between Python 2 and 3

#            current_date = ephem.now()
            year_month = ephem.now().datetime().strftime('%Y_%m')

        if (observation_month in self.epoch):
            print ('** ERROR ** specified year/month is already present, aborting')
            return

        '''
        Code below copied from PMCORB2MonthlyCatalog.py
        This basically does the same thing for every month, so you don't need to use a separate
        file depending on the month. (Instead, you want to update the file periodically to make
        sure the local version of the database is up to date with the remote version)
        
        Has not yet been tested!
        '''

        raw_catalog_list = []

        f = open(self.db_file_and_path, "r")

        for line in f:
            if line == "\n":
                continue       # in case there's a blank
            else:              # line in the original data file
                raw_catalog_list.append(line.rstrip())
        f.close()

        for n in range(100):
            if raw_catalog_list[n] == "----------------------------------------------------------------------------------------------------------------------------------------------------------------":
                start_line = n + 1
# to define the start of the actual data table,
# which comes after ~30 lines of header text

        cropped_catalog_list = []
# crop off the header
        for n in range(len(raw_catalog_list) - start_line):
            cropped_catalog_list.append(raw_catalog_list[n + start_line])

        full_catalog = []

        for obj_mpc in cropped_catalog_list:
            # designation = obj_mpc[0:8].strip()
            abs_m_H = obj_mpc[8:14].strip()
            slope_G = obj_mpc[14:20].strip()
            epoch = obj_mpc[20:26].strip()
            mean_anomaly_M = obj_mpc[26:36].strip()
            peri = obj_mpc[37:47].strip()
            node = obj_mpc[48:58].strip()
            inclin = obj_mpc[59:69].strip()
            eccen = obj_mpc[70:80].strip()
            motion_n = obj_mpc[80:92].strip()
            a = obj_mpc[92:104].strip()
            unc_U = obj_mpc[105:107].strip()
            # reference = obj_mpc[107:117].strip()
            # num_obs = obj_mpc[117:123].strip()
            # num_ops = obj_mpc[123:127].strip()
            # opposition_observ_span = obj_mpc[127:137].strip()
            # rms = obj_mpc[137:142].strip()
            # coarse_pert = obj_mpc[142:146].strip()
            # precise_pert = obj_mpc[146:150].strip()
            # computer = obj_mpc[150:161].strip()
            readable_designation = obj_mpc[166:194].strip()
            # date_last_obs = obj_mpc[194:203].strip()
            """
            MPC format has a "packed" date, allowing the epoch to be stored in 
            fewer digits. However, this must be converted to mm/dd/yyyy format 
            for XEphem.
            """
            # month
            if epoch[3:4] == "1": epoch_x = "01/"
            elif epoch[3:4] == "2": epoch_x = "02/"
            elif epoch[3:4] == "3": epoch_x = "03/"
            elif epoch[3:4] == "4": epoch_x = "04/"
            elif epoch[3:4] == "5": epoch_x = "05/"
            elif epoch[3:4] == "6": epoch_x = "06/"
            elif epoch[3:4] == "7": epoch_x = "07/"
            elif epoch[3:4] == "8": epoch_x = "08/"
            elif epoch[3:4] == "9": epoch_x = "09/"
            elif epoch[3:4] == "A": epoch_x = "10/"
            elif epoch[3:4] == "B": epoch_x = "11/"
            elif epoch[3:4] == "C": epoch_x = "12/"

            # day
            if epoch[4:5] == "1": epoch_x = epoch_x + "01.0/"
            elif epoch[4:5] == "2": epoch_x = epoch_x + "02.0/"
            elif epoch[4:5] == "3": epoch_x = epoch_x + "03.0/"
            elif epoch[4:5] == "4": epoch_x = epoch_x + "04.0/"
            elif epoch[4:5] == "5": epoch_x = epoch_x + "05.0/"
            elif epoch[4:5] == "6": epoch_x = epoch_x + "06.0/"
            elif epoch[4:5] == "7": epoch_x = epoch_x + "07.0/"
            elif epoch[4:5] == "8": epoch_x = epoch_x + "08.0/"
            elif epoch[4:5] == "9": epoch_x = epoch_x + "09.0/"
            elif epoch[4:5] == "A": epoch_x = epoch_x + "10.0/"
            elif epoch[4:5] == "B": epoch_x = epoch_x + "11.0/"
            elif epoch[4:5] == "C": epoch_x = epoch_x + "12.0/"
            elif epoch[4:5] == "D": epoch_x = epoch_x + "12.0/"
            elif epoch[4:5] == "E": epoch_x = epoch_x + "14.0/"
            elif epoch[4:5] == "F": epoch_x = epoch_x + "15.0/"
            elif epoch[4:5] == "G": epoch_x = epoch_x + "16.0/"
            elif epoch[4:5] == "H": epoch_x = epoch_x + "17.0/"
            elif epoch[4:5] == "I": epoch_x = epoch_x + "18.0/"
            elif epoch[4:5] == "J": epoch_x = epoch_x + "19.0/"
            elif epoch[4:5] == "K": epoch_x = epoch_x + "20.0/"
            elif epoch[4:5] == "L": epoch_x = epoch_x + "21.0/"
            elif epoch[4:5] == "M": epoch_x = epoch_x + "22.0/"
            elif epoch[4:5] == "N": epoch_x = epoch_x + "23.0/"
            elif epoch[4:5] == "O": epoch_x = epoch_x + "24.0/"
            elif epoch[4:5] == "P": epoch_x = epoch_x + "25.0/"
            elif epoch[4:5] == "Q": epoch_x = epoch_x + "26.0/"
            elif epoch[4:5] == "R": epoch_x = epoch_x + "27.0/"
            elif epoch[4:5] == "S": epoch_x = epoch_x + "28.0/"
            elif epoch[4:5] == "T": epoch_x = epoch_x + "29.0/"
            elif epoch[4:5] == "U": epoch_x = epoch_x + "30.0/"
            elif epoch[4:5] == "V": epoch_x = epoch_x + "31.0/"

            # year
            if epoch[0:1] == "I": epoch_x = epoch_x + "18"
            elif epoch[0:1] == "J": epoch_x = epoch_x + "19"
            elif epoch[0:1] == "K": epoch_x = epoch_x + "20"
            epoch_x = epoch_x + epoch[1:3]

            if unc_U == "": unc_U = "?"
            expanded_designation = readable_designation + " " + unc_U

            # Write XEphem format orbit to the full_catalog list.
            full_catalog.append(expanded_designation + ",e," + inclin + ","
                                + node + "," + peri + "," + a + "," + motion_n + "," + eccen + "," +
                                mean_anomaly_M + "," + epoch_x + "," + "2000" + ",H " + abs_m_H +
                                "," + slope_G + "\n")

            # f2 = file("Monthly_Orbit_Catalogs/" + current_year + "_" + current_month +
            f2 = open(self.dir + year_month + "_ORB.DAT", "w")
            for obj in full_catalog:
                f2.write(obj)
            f2.close()

