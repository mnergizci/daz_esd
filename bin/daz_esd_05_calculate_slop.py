#!/usr/bin/env python3
"""
v1.0 2022-01-03 Milan Lazecky, Leeds Uni

This script will estimate velocities from daz values.

===============
Input & output files
===============
Inputs :
 - frames.csv - should have ITRF and iono corr.
 - esds.csv - should have iono corr.

Outputs :
 - esds_final.csv
 - frames_final.csv

=====
Usage
=====
daz_05_calculate_slopes.py [--s1ab] [--indaz esds_with_iono.csv] [--infra frames_with_itrf.csv] [--outfra frames_final.csv] [--outdaz esds_final.csv]

Parameters:
    --s1ab ...... also estimate (and store to outfra) the s1ab offset prior to velocity estimation. Now done only for the noiono+notide (final) daz
    --nosubset .. by default, we limit the dataset to start since March 2016 as it appeared too noisy before. This can be adjusted/cancelled using this switch.
"""
#%% Change log
'''
v1.1 2024-04-06 ML
 - added S1AB offset estimation
v1.0 2022-01-03 Milan Lazecky, Uni of Leeds
 - Original implementation - based on codes from 2021-06-24
'''

import getopt, os, sys
import argparse
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'lib')))
from daz_esd_lib import *
from daz_timeseries_lib import *


class Usage(Exception):
    """Usage context manager"""
    def __init__(self, msg):
        self.msg = msg


#%% Main
def main(argv=None):
    #%% Check argv
    if argv is None:
        argv = sys.argv

    #%% Define argument parser
    parser = argparse.ArgumentParser(description='Process files and options.')
    parser.add_argument('indazfile', help='Input DAZ file')
    parser.add_argument('inframefile', help='Input frame file')
    parser.add_argument('outdazfile', help='Output DAZ file')
    parser.add_argument('outframefile', help='Output frame file')
    parser.add_argument('--s1ab', action='store_true', help='Flag to enable S1AB')
    parser.add_argument('--nosubset', action='store_true', help='Flag to disable subset')
    parser.add_argument('--roll_assist', action='store_true', help='Flag to rolling assist with ITRF14')

    #%% Parse arguments
    try:
        args = parser.parse_args(argv[1:])
        
        indazfile = args.indazfile
        inframefile = args.inframefile
        outdazfile = args.outdazfile
        outframefile = args.outframefile
        s1ab = args.s1ab
        subset = not args.nosubset
        roll_assist = args.roll_assist
    except argparse.ArgumentError as msg:
        raise Usage(msg)

    ####Define the folderpaths:
    #jasmin_path:
    
    LiCS_proc = '/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/proc/current'
    batchdir = os.environ.get('BATCH_CACHE_DIR', './')  
    
    ###local_path:
    #batchdir ='data'

    ##here same with both env
    workdir = os.path.join(batchdir, 'daz_esd')
    homedir = os.getcwd()
    print(workdir)
    frame=indazfile.split('.')[0]
    try:
        # if os.path.exists(os.path.join(workdir, outdazfile)) or os.path.exists(os.path.join(workdir, outframefile)):
            # raise Usage('Output files already exists. Please remove/rename them to rerun the code again.')
        if not os.path.exists(os.path.join(workdir, indazfile)) and os.path.exists(os.path.join(workdir, inframefile)):
            raise Usage('Input daz and frame files do not exist. Please be sure to run the previous codes properly?')
    except Usage as err:
        print("\nERROR:")
        print("  " + str(err.msg))
        print("\nFor help, use -h or --help.\n")
        return 2


    
    ### processing itself:
    indaz=os.path.join(workdir, indazfile)
    inframe=os.path.join(workdir, inframefile)
    esds, framespd = load_csvs(esdscsv = indaz, framescsv = inframe, core_init=True)
    # esds= pd.read_csv(os.path.join(workdir, indazfile))
    # framespd= pd.read_csv(os.path.join(workdir, inframefile))

    # setting 'subset' - means, only data > 2016-03-01 as before it is too noisy
    if subset:
        print('Subsetting dataset to include only data after 2016-03-01')
        esds['epoch'] = pd.to_datetime(esds['epoch'])
        esds = esds[esds['epoch'] > pd.Timestamp('2017-01-01')]
    if s1ab:
        print('Estimating S1AB offset per frame')
        # estimate the offset first, then apply correction, and then use Huber as usual
        framespd = estimate_s1ab_allframes(esds, framespd, col = 'daz_bovl_notide_noiono', rmsiter = 50)
        print('Applying S1AB corrections (only to daz_mm_notide_noiono and stored as daz_mm_final)')
        esds['daz_mm_final'] = esds['daz_bovl_notide_noiono'].copy()
        esds, framespd = correct_s1ab(esds, framespd, cols=['daz_mm_final'])
    # 2021-10-12: the original way:
    for col in ['daz_bovl_notide_noiono', 'daz_bovl', 'daz_iono_mm', 'daz_tide_mm', 'daz_esd2_bovl_notide_noiono']:
        if col in esds:
            print('estimating velocities of '+col)
            esds, framespd = df_calculate_slopes(esds, framespd, alpha = 1, eps = 1.35, bycol = col, subset = subset, roll_assist = roll_assist)
    # to back up before continuing:
    print('saving datasets')
    framespd.to_csv(os.path.join(workdir,outframefile))
    esds.to_csv(os.path.join(workdir,outdazfile))
    print('done')

#%% main
if __name__ == "__main__":
    sys.exit(main())