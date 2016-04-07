# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from netCDF4 import Dataset
import glob
import numpy as np
from SoundingFilter import sndfilter
from datetime import datetime

# <codecell>

'''
    This code reads in ARM netCDF data for a radiosonde and then passes
    the data to the sounding filter program.  The sounding filter program
    pulls out the standard and significant levels and returns the new arrays
    that are then stored as netCDF files with some different metadata.

    This program is meant to find all of the ARM netCDF files in a directory
    and then run the filter on them.

'''

def makeFile(fn, bt, lat, lon, src, **snd):
    # Function to save the filtered data to a new netCDF file
    output = Dataset(fn, 'w', format='NETCDF3_CLASSIC')
    output.createDimension('time', None)
    dt_str = datetime.strftime(datetime.now(), "%Y-%m-%d %H:%M:%S")
    output.Date_netCDF_file_created = dt_str
    output.source = src # original netCDF file
    #output.filter_version = sndfilter.githash
    output.filter_version = '187343787dbb3a73c8fef376021dd753a0d78568'
    output.author = 'Greg Blumberg'
    output.author_affiliation = 'OU/CIMMS'
    output.contact = 'wblumberg@ou.edu'
    output.comment = "These are files containing sounding that that was filtered using the standard and significant level filter developed by Tim Supinie and Greg Blumberg." 
    output.filter_link = 'https://github.com/tsupinie/SoundingFilter'

    var = output.createVariable('base_time', 'i4')
    var[:] = bt
    var.units = "seconds since 1970-01-01 00:00:00+00:00"

    var = output.createVariable('time_offset', 'i4', ('time',))
    var[:] = np.arange(0,len(snd['wspd'])*2, 2)
    var.units = 'seconds since the value in base_time (a dummy var)'

    var = output.createVariable('pres', 'f4', ('time',))
    var[:] = snd['pres']
    var.units = 'mb'
    var.comment = 'air pressure'

    var = output.createVariable('dp', 'f4', ('time',))
    var[:] = snd['dewp']
    var.units = 'Celsius'
    var.comment = 'Dewpoint temperature'

    var = output.createVariable('tdry', 'f4', ('time',))
    var[:] = snd['temp']
    var.units = 'Celsius'
    var.comment = 'dry bulb temperature'

    var = output.createVariable('wspd', 'f4', ('time',))
    var[:] = snd['wspd']
    var.units = 'kts'
    var.comment = 'wind speed'

    var = output.createVariable('deg', 'f4', ('time',))
    var[:] = snd['wdir']
    var.units = 'deg'
    var.comment = 'wind direction from North'
    print lat
    var = output.createVariable('lat', 'f4')
    var[:] = lat[0]
    var.units = 'decimal degrees'
    var.comment = 'latitude'

    var = output.createVariable('lon', 'f4')
    var[:] = lon[0]
    var.units = 'decimal degrees'
    var.comment = 'longitude'

    var = output.createVariable('alt', 'f4', ('time',))
    var[:] = snd['hght']
    var.units = 'm above MSL'
    var.comment = 'altitude'

    output.close()

# A dictionary that converts the variable names from the netCDF file
# into something that the sounding filter can understand.
name_translate = {
    'tdry':'temp',
    'dp':'dewp',
    'pres':'pres',
    'wspd':'wspd',
    'deg':'wdir',
    'alt':'hght'}

sgp_files = np.sort(glob.glob('/raid/FRDD/Dave.Turner/data/sgp/sonde/*.cdf'))
for i in sgp_files: # loop over all the ARM netCDF files found
    snd_file = Dataset(i)

    bt = snd_file.variables['base_time'][:]
    lat = snd_file.variables['lat'][:]
    lon = snd_file.variables['lon'][:]
    src = i

    # Pull out the variables from the netCDF file that correspond to different profiles
    snd = dict( (name_translate[k], snd_file.variables[k][:]) for k in ['tdry', 'dp', 'pres', 'wspd', 'deg', 'alt'])
    snd_file.close()
    try:
        # Run the filter
        snd_filtered = sndfilter.soundingFilter(**snd)
    except Exception,e:
        # if there's an issue, print out the exception and just continue onto the next file
        print i
        print e
        continue

    # Make the filename and actually save the data to a new file.
    dates = '.'.join(i.split('/')[-1].split('.')[2:4])
    print dates
    dt = datetime.strptime(dates, '%Y%m%d.%H%M%S')
    print dt
    fn = 'ounsondessl1blumC1.c1.' + dates + '.cdf'
    
    print i, len(snd_filtered['temp'][:])
    makeFile(fn, bt, lat, lon, src, **snd_filtered) 

