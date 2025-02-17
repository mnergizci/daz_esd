#!/usr/bin/env python3
"""
v1.1 Muhammet Nergizci, University of Leeds
v1.0 2022-01-03 Milan Lazecky, Leeds Uni
https://github.com/comet-licsar/daz

This script will calculate solid Earth tides (currently using 'solid' implemented in GMT earthtides) and merge with esds.txt.
If the SET file exists, it will just merge it.

===============
Input & output files
===============
Inputs :
 - $frame.csv - it should be output from daz_esd_01_prepa
Outputs :
 - esds.csv - contains data with heading:
frame,orbits_precision(!),daz_tide_mm,epochdate,daz_mm,years_since_beginning,daz_mm_notide
 - tides.csv - contains data with heading:
frame, epoch, dEtide, dNtide, dUtide

=====
Usage
=====
daz_02_extract_SET.py [--indaz esds.txt] [--tidescsv tides.csv] [--outdaz esds.csv]
daz_esd_02_extract_SET.py 027A_04887_262625.csv  earthtides.csv  027A_04887_262625.tide.csv

 --tidescsv - input or output (if does not exist) file containing SET.

"""

import getopt, os, sys
import pandas as pd
# Add the parent directory to the sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'lib')))
from daz_esd_lib import *

class Usage(Exception):
    """Usage context manager"""
    def __init__(self, msg):
        self.msg = msg

#%% Main
def main(argv=None):
    
    #%% Check argv
    if argv is None:
        argv = sys.argv

    #%% Read options
    try:
        opts, args = getopt.getopt(argv[1:], "h", ["help"])
        if len(args) < 3:
            raise Usage("Frame identifier is required as a positional argument.")
        
        indazfile = args[0]
        tidescsv = args[1]
        outdazfile= args[2]

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
    print(indazfile, tidescsv, outdazfile)
    print(os.path.join(workdir, indazfile))
    try:
        # if os.path.exists(os.path.join(workdir, outdazfile)):
        #     raise Usage('Output esds csv file already exists. Please remove it to rerun the code again.')
        if not os.path.exists(os.path.join(workdir, indazfile)):
            raise Usage('Input daz_esd values do not exist. Please be sure to run daz_esd_01_prepare_input.py.')
    except Usage as err:
        print("\nERROR:")
        print("  " + str(err.msg))
        print("\nFor help, use -h or --help.\n")
        return 2


    # Read input data
    df_temp = pd.read_csv(os.path.join(workdir, indazfile))
    # df_temp_sorted = df_temp[df_temp['swath_overlap'].str.startswith('2-8')]
    # df_temp_sorted.reset_index(drop=True, inplace=True)
    #df_inp=df_temp_sorted[['epoch', 'master', 'heading', 'center_time', 'center_lon', 'center_lat']]
    df_inp=df_temp[['epoch', 'master', 'heading', 'center_time', 'center_lon', 'center_lat']]
    output_path = os.path.join(workdir, frame + '.SET.csv')
    df_inp.to_csv(output_path, index=False)

    #processing itself
    if not os.path.exists(os.path.join(workdir,tidescsv)):
        # print(args[0], args[1], args[2])
        cmd= f"bash get_set_esd.sh {os.path.join(frame + '.SET.csv')} {tidescsv}"
        os.system(cmd)
    else:
        print('SET file already exists. Will use it for merging')

    if not os.path.exists(os.path.join(workdir,tidescsv)):
        print('ERROR - the SET file was not generated, exiting')
        return 2


    # now (finally) load esds and tides to python, merge etc.
    earthtides = pd.read_csv(os.path.join(workdir,tidescsv))
    esds = merge_tides(df_temp,earthtides)
    esds.to_csv(os.path.join(workdir, outdazfile), index=False)

#%% main
if __name__ == "__main__":
    sys.exit(main())