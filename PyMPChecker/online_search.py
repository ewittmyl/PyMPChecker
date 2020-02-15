# Copied from Krzysztof's code at
# https://github.com/GOTO-OBS/goto/blob/master/scripts/MPchecker/mpchecker.py
#
# Duncan.Galloway@monash.edu 2019 Sep 19
#

import sys
import glob
from astropy.wcs import WCS
from astropy.io import fits
import time
from astropy.coordinates import Angle
from astropy import units
import requests
import numpy as np



def find_asteroid(sciim,ra0=None,dec0=None,r=110,datestr=None,debug=False):
            # print(sciim,ra0,dec0)
            
            url_base='http://cgi.minorplanetcenter.net/cgi-bin/mpcheck.cgi'
            
            if not datestr:
                hdulist=fits.open(sciim)
                n=0
                for n in range(len(hdulist)):
                    head=hdulist[n].header
                    if 'DATE-OBS' in head: 
                        break
                hdulist.close()
                if not 'DATE-OBS' in head: 
                    print('No date in header!')
                    return 'Date header problem!'
                datestr=head['DATE-OBS']

                if ra0==None or dec==None:
                    w=WCS(head)
                    xsize,ysize=head['NAXIS1'],head['NAXIS2']
                    
                    ra0,dec0=w.all_pix2world(xsize/2.,ysize/2.,1)
                    print('Field centre:',ra0,dec0)
                    print('Search radius:',r)
                
                
            
            #print(head['DATE-OBS'])date=time.strptime
            for fmt in ("%Y-%m-%dT%H:%M:%S.%f","%Y-%m-%dT%H:%M:%S"):
                    try:
                        date=time.strptime(datestr,fmt)
                    except ValueError:
                        #print(fmt,'bad')
                        continue
                    break
            
            if debug:
                    print(datestr,date)
            dayfr=date.tm_mday+(date.tm_hour+date.tm_min/60+date.tm_sec/3600)/24
            fulldate='%d %d %.2f'%(date.tm_year,date.tm_mon,dayfr)
            ra=Angle(float(ra0)/180*12, unit=units.hour).hms#.to_string(sep=' ',pad=True,precision=2)
            dec=Angle(float(dec0), unit=units.deg).dms#.to_string(sep=' ',pad=True,precision=1)
            
            #print(ra.h,ra.m,ra.s,dec.d,dec.m,dec.s)
            url=url_base+'?year='+str(date.tm_year)+'&month='+str(date.tm_mon)+'&day='+'%.2f'%(dayfr)+\
            '&which=pos&ra='+'%02d'%(ra.h)+'%20'+'%02d'%(ra.m)+'%20'+'%05.2f'%(ra.s)+\
            '&decl='+'%+03.0f'%(dec.d)+'%20'+'%02d'%(abs(dec.m))+'%20'+'%04.1f'%(abs(dec.s))+\
            '&radius='+str(r)+'&limit=24&TextArea=&oc=950&sort=d&mot=h&tmot=s&pdes=u&needed=f&ps=n&type=p' 

            if debug:
                    print(url)
            
            #print(url)
            #exit(0)
            while True:
                try:
                    contents=requests.get(url).text
                except:
                    print('Repeating mpc query for %.4f %.4f'%(ra0,dec0))
                    time.sleep(5)
                    continue
                break
            # print (contents)
# Results of the query are a webpage, as follows; the rest of this code is 
# concerned with parsing the result rows
# <html>
# <head>
# <title>MPChecker/CMTChecker/NEOChecker/NEOCMTChecker</title>
# </head>
# <body>
# <h1>MPChecker/CMTChecker/NEOChecker/NEOCMTChecker</h1>
# Here are the results of your search(es) in the requested field(s)
#  (positions are determined from elements integrated to a nearby epoch)
# :<p><hr><p>
# The following objects, brighter than <i>V</i> = 24, were found in the  10.0-arcminute region around R.A. = 15 10 43.62, Decl. = +33 21 53.7 (J2000.0) on 2019 09 17.86 UT:
#  
#  <pre>
#  Object designation         R.A.      Decl.     V       Offsets     Motion/hr   Orbit  <a href="https://www.minorplanetcenter.net/iau/info/FurtherObs.html">Further observations?</a>
#                            h  m  s     &#176;  '  "        R.A.   Decl.  R.A.  Decl.        Comment (Elong/Decl/V at date 1)
# 
# (325395) 2009 CQ5        15 10 43.5 +33 21 45  18.6   0.0W   0.1S    16+    69+    6o  NEO : Desirable between 2019 Sept. 18-Oct. 18.  At the first date, object will be within 60 deg of the sun.
#  </pre>
#  <p><hr><p>
#  Number of objects checked =  819915
#  <p><hr><p>
#  <h2>Explanatory Notes</h2>
#  <ul>
#  <li>The positions are J2000.0 and are "quick look" positions designed for
#      identification, not the rigorous comparison of observations with theory.
#  <li>For recent dates (within 40 days of now), positions of NEOs are generated f
#  rom elements for the
#      nearest 0h UT epoch.  For older dates, positions are generated from element
#  s
#      for the nearest 40-day epoch date.
#  <li>Offsets, intended for use by supernova hunters, are given in arc-minutes
#      as the coordinates of the parent galaxies are rarely
#      given to arc-second precision.
#  <li>The motions are in arcseconds per stated time unit (if minutes or hours)
#      or degrees/day.
#  <li>Right-ascension motions include the cos(Decl.) term.
#  <li>The brief orbit descriptor is either:
#    <ul>
#    <li> the number of oppositions (if marked with `o'),
#    <li> the arc-length in days (if marked with `d') or
#    <li> `V' if it is a Generalized V&#228;is&#228;l&#228; solution.
#    </ul>
#  <li>Comets are listed regardless of how faint they are.  No magnitude
#      estimates are supplied for comets.  The heliocentric distance, r,
#      is displayed in the Comments column.  Most comets more than 5 AU
#      from the sun will be beyond the reach of most observers, but the
#      information is displayed to enable identification should an outburst
#      occur.
#  <li>If you are requesting objects other than just planets and natural
#      satellites, a count is displayed of the number of objects that were
#      checked.  If the count is less than 700000 (16000 for the NEO and NEOCMT
#      checker, and 800 for the CMT checker) and the date is after 1946,
#      then the files of elements
#      used by this service may have been truncated and this fact should
#      be <a href="http://cgi.minorplanetcenter.net/cgi-bin/feedback.cgi?U=cgi/Che
#  ckMP&amp;S=Truncated%20minor-planet%20elements">reported</a>.
#      Note that only the first 500 numbered minor planets are available for
#      checking pre-1900 dates.
#  <li>The information that is displayed for each object for pre-2009 dates is
#      much more limited than it is for current dates.
#  </ul>
# <p><hr><p>
# <p><center>These calculations have been performed on the <a href="https://www.minorplanetcenter.net/iau/Ack/TamkinFoundation.html">Tamkin Foundation Computing Network</a>.</center>
# <p><hr><p>
# </body>
# </html>
                
            indx=contents.find('<pre>')
            #print(contents)
            #exit(0)
            if indx>0 and indx<len(contents):
                indx=contents.find('Object',indx)
                #indx=contents.find('\n',indx)
                #indx=contents.find('\n',indx+1)
                #indx=contents.find('\n',indx+1)
                indx2=contents.find('</pre>',indx)
                #print(contents[indx:indx2].split('\n'))
                result=[]
                
                for line in contents[indx:indx2].split('\n')[3:-1]:
                    c=line.split()
                    #print(c)
                    try:
                        decindx=[i for i, s in enumerate(c[4:10]) if '+'  in s or '-' in s][0]+4
                    except:
                        print('Cannot find declination sign!')
                        decindx=4
                    #print(decindx)
                    ra=(float(c[decindx-3])+float(c[decindx-2])/60.+float(c[decindx-1])/3600.)*180/12.
                    
                    if '-' in c[decindx]:
                        sign=-1
                    else: 
                        sign=1
                    dec=(float(c[decindx])+sign*float(c[decindx+1])/60.+sign*float(c[decindx+2])/3600.)
                    #dec=float(c[decindx])
                    #print(dec)
                    imag = 0
                    try:
                        mag=float(c[decindx+3])
                        imag = 1
                    except:
                        mag=99.999
                    ra_off=c[decindx+3+imag]
                    dec_off=c[decindx+4+imag]
                    # convert offset to arcsec
                    dist = 60.*np.sqrt(float(ra_off[:-1])**2+float(dec_off[:-1])**2)
                    comment=' '.join(c[decindx+8+imag:])

                    result.append([' '.join(c[:decindx-3]),dist,ra,dec,mag,ra_off,dec_off,comment])

                return result
                #return(contents[indx:indx2].split('\n')[3])
            else:
                return None
            #exit(0)
            
