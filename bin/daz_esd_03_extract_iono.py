#!/usr/bin/env python3
"""
v1.1 2024 Muhammet Nergizci, Leeds
v1.0 2022-01-03 Milan Lazecky, Leeds Uni

This script will extract ionosphere shifts using either IRI2016 model or the CODE dataset.
===============
Input & output files
===============
Inputs :
- esds.csv - contains data with heading:

Outputs :
 - esds_with_iono.csv - added iono columns
=====
Usage
=====
daz_esd_03_extract_iono.py [--indaz esds.csv] [--use_code] [--outdaz esds_with_iono.csv]
python daz_esd_03_extract_iono.py 027A_04887_262625.tide.csv 027A_04887_262625.ion.csv

Notes:
    --use_code  Will apply CODE to get TEC values rather than the default IRI2016 estimates. Note IRI2016 is still used to estimate iono peak altitude. Tested only in LiCSAR environment.
"""
import getopt, os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'lib')))
from daz_esd_lib import *
from daz_iono_lib import *


# keeping assumption of hei=450 km --- best to test first (yet at the moment we are anyway quite coarse due to 1-value-per-frame)
use_iri_hei = False

class Usage(Exception):
    """Usage context manager"""
    def __init__(self, msg):
        self.msg = msg


#%% Main
def main(argv=None):
    
    #%% Check argv
    if argv == None:
        argv = sys.argv
    
    #%% Set default
    ionosource = 'code'

    #%% Read options
    try:
        opts, args = getopt.getopt(argv[1:], "h", ["help"])
        if len(args) < 2:
            raise Usage("Frame identifier is required as a positional argument.")

        indazfile = args[0]
        outdazfile= args[1]
        
        for o, a in opts:
            if o == '-h' or o == '--help':
                print(__doc__)
                return 0
          
    except getopt.error as msg:
        raise Usage(msg)

    ####Define the folderpaths:
    LiCS_proc = '/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/proc/current'
    batchdir = os.environ.get('BATCH_CACHE_DIR', './')  # Default to current directory if not set
    workdir = os.path.join(batchdir, 'daz_esd')
    homedir = os.getcwd()
    frame=indazfile.split('.')[0]
    
    try:
        if os.path.exists(os.path.join(workdir, outdazfile)):
            raise Usage('Output iono_esds csv file already exists. Please remove it to rerun the code again.')
        if not os.path.exists(os.path.join(workdir, indazfile)):
            raise Usage('Input daz_esd values do not exist. Please be sure to run daz_esd_01_prepare_input.py and daz_esd_02_extract_SET.py.')
    except Usage as err:
        print("\nERROR:")
        print("  " + str(err.msg))
        print("\nFor help, use -h or --help.\n")
        return 2

    ###open the $frame.tide.csv file
    esds=pd.read_csv(os.path.join(workdir, indazfile))
    framespd= esds.drop_duplicates(subset='swath_overlap')
    framespd.reset_index(drop=True, inplace=True)
    framespd=framespd[['swath_overlap', 'master', 'center_lat', 'center_lon', 'center_time','center_range','az_pix', 'heading', 'incidence', 'dfDC', 'ka']]
    esds, framespd = extract_iono_full(esds, framespd, ionosource = ionosource, use_iri_hei=use_iri_hei)
            
    if 'daz_bovl_notide' in esds.columns:
        for col in ['daz_bovl_notide', 'daz_esd2_bovl_notide']:
            esds[col + '_noiono'] = esds[col] - esds['daz_iono_mm']
    else:
        col = 'daz_bovl'
        try:
            esds[col + '_noiono'] = esds[col] - esds['daz_iono_mm']  # 2023/08: changed sign to keep consistent with the GRL article
        except KeyError:
            print('probably a bug, please check column names - in any case, the correction is stored as daz_iono_mm column')


    print('saving files')
    esds.to_csv(os.path.join(workdir, outdazfile))
    framespd.to_csv(os.path.join(workdir, frame+'.frame.step3.csv'))
    print('done')

#%% main
if __name__ == "__main__":
    sys.exit(main())