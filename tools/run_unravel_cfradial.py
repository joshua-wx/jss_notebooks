import os #used for system commands
import argparse
import warnings
# Just ignoring warning messages.
warnings.simplefilter("ignore")
import tempfile #used to create temporary folders to store data
import zipfile #used to extract tar files
import urllib #used to download data via http
from datetime import datetime #used to manipulate time :)
from glob import glob #used for manipulating pathnames

import pyart
from unravel import dealias


def main():
    
    #read file list
    in_ffn_list = sorted(glob(IN_FOLDER + '/*.nc'))
    
    for i, ffn in enumerate(in_ffn_list):
        #read into pyart
        radar = pyart.io.read(ffn)
        #processing using unravel
        ultimate_vel = dealias.process_3D(radar, nyquist_velocity=NI, velname="VRADH", dbzname="DBZH")
        #add field
        radar.add_field_like("VRADH", "corrected_velocity_3D", ultimate_vel, replace_existing=True)
        #save to file
        output_ffn = OUT_FOLDER + '/' + os.path.basename(ffn)[:-3] + '_unravel.nc'
        pyart.io.write_cfradial(output_ffn ,radar)
        #update user
        print('finished file', str(i), 'of', str(len(in_ffn_list)))
    print('completed unravel processing')
        
if __name__ == '__main__':
    """
    Global vars
    """    
    
    # Parse arguments
    parser_description = "unravel cfradial"
    parser = argparse.ArgumentParser(description=parser_description)
    parser.add_argument(
        '-i',
        '--input',
        dest='input',
        default=None,
        type=str,
        help='input folder of cfradial',
        required=True)
    parser.add_argument(
        '-o',
        '--output',
        dest='output',
        default=None,
        type=str,
        help='output folder of unravelled cfradial',
        required=True)
    parser.add_argument(
        '-n',
        '--nyquist',
        dest='nyquist',
        default=None,
        type=float,
        help='nyquist velocity of cfradial',
        required=True)
    args = parser.parse_args()
    IN_FOLDER  = args.input
    OUT_FOLDER = args.output
    NI = args.nyquist
    
    if not os.path.exists(OUT_FOLDER):
        os.mkdir(OUT_FOLDER)
    
    with warnings.catch_warnings():
        # Just ignoring warning messages.
        warnings.simplefilter("ignore")
        main()