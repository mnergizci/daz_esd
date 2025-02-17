#!/usr/bin/env python3
"""
v1.0 2022-01-03 Milan Lazecky, Leeds Uni

This script will extract plate motion model values from ITRF2014 PMM (extracted from UNAVCO website).
It will extract values in a buffer around given frame centre coordinate, and then average them.

===============
Input & output files
===============
Inputs :
 - frames.csv - either orig frames.csv or the one including iono columns
[- vel_gps.nc - GPS velocities, must contain VEL_E, VEL_N variables, expected in NNR]

Outputs :
 - frames_with_itrf.csv - new columns with the ITRF2014 PMM in azimuth direction

=====
Usage
=====
daz_04_extract_PMM.py [--add_eu] [--infra frames.csv] [--outfra frames_with_itrf.csv] [--velnc vel_gps_kreemer.nc] 

Note: param velnc is optional, but if provided as nc file with VEL_E, VEL_N variables, it will be used as GPS velocities.

"""
#%% Change log
'''
v1.0 2022-01-03 Milan Lazecky, Uni of Leeds
 - Original implementation - based on codes from 2021-06-24
'''

import getopt, os, sys
import argparse
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

    #%% Define argument parser
    parser = argparse.ArgumentParser(description='Process files and optionally add EU flag.')
    parser.add_argument('inionfile', help='Input ION file')
    parser.add_argument('outpmmfile', help='Output PMM file')
    parser.add_argument('velncfile', help='Velocity NC file')
    parser.add_argument('--add_eu', action='store_true', help='Flag to add EU')
    # parser.add_argument('-h', '--help', action='help', help='Show this help message and exit')

    #%% Parse arguments
    try:
        args = parser.parse_args(argv[1:])
        
        inionfile = args.inionfile
        outpmmfile = args.outpmmfile
        velncfile = args.velncfile
        add_eu = args.add_eu
        
    except argparse.ArgumentError as msg:
        raise Usage(msg)

    # Your processing function should be called here

    ####Define the folderpaths:
    LiCS_proc = '/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/proc/current'
    batchdir = os.environ.get('BATCH_CACHE_DIR', './')  # Default to current directory if not set
    workdir = os.path.join(batchdir, 'daz_esd')
    homedir = os.getcwd()
    frame=inionfile.split('.')[0]    
    #print(inionfile, outpmmfile, velncfile)    
    try:
        # if os.path.exists(os.path.join(workdir, outpmmfile)):
        #     raise Usage('Output PMM file already exists. Please remove/rename it to rerun the code again.')
        if not os.path.exists(os.path.join(workdir, inionfile)):
            raise Usage(f'Input {inionfile} values do not exist. Please recheck again?')
    except Usage as err:
        print("\nERROR:")
        print("  " + str(err.msg))
        print("\nFor help, use -h or --help.\n")
        return 2


    # get plate motion model
    ################### ITRF2014
    #takes again long - some 2-3 hours
    print('getting plate motion model values using the average for the 222x222 km around the frame centre')
    print('(using ITRF2014 and external nc file for GPS, if available)')
    framespd=pd.read_csv(os.path.join(workdir, inionfile))
    velnc=os.path.join(workdir, velncfile)
    # print(velnc)
    framespd = df_get_itrf_gps_slopes(framespd, velnc=velnc, add_eu= add_eu)
    if add_eu:
        print('Finished extracting both NNR and EU PMM values, note:')
        print("framespd['vel_eur'] = framespd['slope_from_daz'] - framespd['slope_vel_itrf_nnr'] + framespd['slope_vel_itrf_eu']")
    framespd.to_csv(os.path.join(workdir, outpmmfile), index=False)
    print('done')

#%% main
if __name__ == "__main__":
    sys.exit(main())


'''
some notes that might be relevant:
gridagg['GPS_N_2014_EU'] = gridagg['GPS_N'] - gridagg['ITRF_N_2008'] + gridagg['ITRF_2014_EU_N']
gridagg['GPS_E_2014_EU'] = gridagg['GPS_E'] - gridagg['ITRF_E_2008'] + gridagg['ITRF_2014_EU_E']
'''