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

def main():
    
    #parse inputs
    radar_id_str = str(RID).zfill(2) #convert radar id to a string and fill with a leading 0 if only one digit
    date_dt      = datetime.strptime(DT_STR, '%Y%m%d')

    #build request filename url
    tar_fn       = radar_id_str + '_' + date_dt.strftime('%Y%m%d') + '.pvol.zip'
    request_url  = '/'.join([BASE_URL, 'odim_pvol', radar_id_str, date_dt.strftime('%Y'), 'vol', tar_fn])
    print('my request is ',request_url)

    #download the zip file
    if not os.path.isfile(tar_fn):
        urllib.request.urlretrieve(request_url, tar_fn)

    #extract the zip file to a temporary directory
    temp_dir = tempfile.mkdtemp()
    zip_fh = zipfile.ZipFile(tar_fn)
    zip_fh.extractall(path = temp_dir)
    zip_fh.close()

    #list all the volumes extracted from the zip file
    file_list = sorted(glob(temp_dir + '/*.h5'))

    #read using pyart and save to cfradial
    for file in file_list:
        #read file
        out_ffn = OUT_FOLDER + '/' + os.path.basename(file)[:-3] + '.nc'
        #output original
        radar = pyart.aux_io.read_odim_h5(file, file_field_names=True)
        pyart.io.write_cfradial(out_ffn, radar)

    print('finished processing:', tar_fn)
    
    os.system('rm -rf ' + tar_fn)
    os.system('rm -rf ' + temp_dir)
    
    
if __name__ == '__main__':
    """
    Global vars
    """    
    #path config
    BASE_URL     = 'http://dapds00.nci.org.au/thredds/fileServer/rq0' #base url for NCI dataset
    
    # Parse arguments
    parser_description = "request and convert odhimh5 to cfradial"
    parser = argparse.ArgumentParser(description=parser_description)
    parser.add_argument(
        '-d',
        '--date',
        dest='date',
        default=None,
        type=str,
        help='date to extract from archive in YYYYMMDD',
        required=True)
    parser.add_argument(
        '-r',
        '--rid',
        dest='rid',
        default=None,
        type=int,
        help='Radar ID',
        required=True)
    parser.add_argument(
        '-o',
        '--output',
        dest='output',
        default=None,
        type=str,
        help='output folder',
        required=True)
    args = parser.parse_args()
    RID     = args.rid
    DT_STR  = args.date
    OUT_FOLDER = args.output
    
    if not os.path.exists(OUT_FOLDER):
        os.mkdir(OUT_FOLDER)
    
    with warnings.catch_warnings():
        # Just ignoring warning messages.
        warnings.simplefilter("ignore")
        main()