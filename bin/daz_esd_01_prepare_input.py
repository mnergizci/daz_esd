#!/usr/bin/env python3
"""
v1.0 Muhammet Nergizci, University of Leeds

This script processes LiCSAR data, calculates along-track ICC and ESD offsets, and save the data to csv file in your $BATCHdir.


Update:
The daz_total is from 
===============
Input & output files
===============
Inputs:
 - frame_id which you interest

Outputs:
 - CSV file containing processed data

=====
Usage
=====
daz_esd.py $frame [--icc=0] [--esd_it2_flag=0]

 --icc - flag to include icc_offset in calculations (default 0)
 --esd_it2_flag - flag to include only esd_it2 in calculations (default 0)
"""

import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import warnings
import getopt
import sys
import argparse
from daz_lib_licsar import *
# Add the parent directory to the sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'lib')))
from daz_esd_lib import *
import framecare as fc
warnings.filterwarnings("ignore")


class Usage(Exception):
    """Usage context manager"""
    def __init__(self, msg):
        self.msg = msg


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser(description='Process some data.')
    parser.add_argument('frame', type=str, help='Frame identifier')
    parser.add_argument('--orbdiff_fix', action='store_true', help='Enable orbit difference fix')
    parser.add_argument('--daz_wrt', action='store_true', help='Enable DAZ write')

    args = parser.parse_args(argv[1:])

    #print(f"args: {args}")

    daz_wrt = args.daz_wrt
    orbdiff_fix = args.orbdiff_fix
    frame = args.frame

    # Print the values of the flags for debugging
    #print(f"daz_wrt: {daz_wrt}, orbdiff_fix: {orbdiff_fix}")

    # Call the processing function
    process_data(frame, daz_wrt, orbdiff_fix)

###colors for printing out!

def process_data(frame, daz_wrt, orbdiff_fix):
    # Locate the LiCSAR_proc directory where the log and results file database.
    LiCS_proc = '/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/proc/current'
    batchdir = os.environ.get('BATCH_CACHE_DIR', './')  # Default to current directory if not set
    track = frame[:4]
    track_id = str(int(track[:-1]))
    folder_path = os.path.join(LiCS_proc, track_id, frame)
    
    # Define DataFrame in pandas
    df = pd.DataFrame(columns=['epoch','master', 'esd_master', 'swath', 'overlap','swath_overlap', 'ph_mean','icc_offset','esd_it1','esd_it2','scale_factor','az_pix','ka','dfDC','daz_bovl', 'S1AorB'])
    
    
    for subfolder in ['log', 'SLC']:
        subfolder_path = os.path.join(folder_path, subfolder)
        if not os.path.isdir(subfolder_path):
            continue
        if subfolder == 'log':
            log_path = subfolder_path
            frame = os.path.basename(os.path.dirname(subfolder_path))
            df = results_files_process(log_path, df, frame)
            log_files_process(log_path, df)
        if subfolder == 'SLC':
            slc_path = subfolder_path
            prime_epoch = os.listdir(slc_path)[0]
            prime_path = os.path.join(slc_path, prime_epoch)
            scale_factors, az_pix, ka, dfDC = get_dfDC_slc(prime_path, f0=5405000500, burst_interval=2.758277)
            for i in range(3):
                scale_factor_temp = float(scale_factors[i])
                az_pix_temp = float(az_pix[i])
                kas_temp=float(ka[i])
                dfDC_temp=float(dfDC[i])
                dataframe_to_update = df[(df['epoch'].str.startswith(prime_epoch)) & (df['swath'] == str(i + 1))]
                for index, row in dataframe_to_update.iterrows():
                    df.at[index, 'az_pix'] = az_pix_temp
                    df.at[index, 'scale_factor'] = scale_factor_temp
                    df.at[index, 'ka']=kas_temp
                    df.at[index, 'dfDC']=dfDC_temp

            # Extract center range values from SLC files
            center_ranges = []
            for sw in [1, 2, 3]:
                slc_file = os.path.join(prime_path, f'{prime_epoch}.IW{sw}.slc.par')
                if os.path.exists(slc_file):
                    near_range = slc_process(slc_file, "near_range")
                    far_range = slc_process(slc_file, "far_range")
                    if near_range is not None and far_range is not None:
                        center_range = (near_range + far_range) / 2
                        center_ranges.append(center_range)
                    else:
                        center_ranges.append(None)
            # Ensure that center_ranges has the correct number of elements
            if len(center_ranges) != 3:
                raise ValueError("Expected 3 center ranges, got {}".format(len(center_ranges)))
            # Update dataframe with center ranges
            for sw in range(1, 4):
                dataframe_to_update = df[(df['epoch'].str.startswith(prime_epoch)) & (df['swath'] == str(sw))]
                for index, row in dataframe_to_update.iterrows():
                    df.at[index, 'center_range'] = center_ranges[sw-1] if len(center_ranges) > sw-1 else None

            #Obtain center_time(start_time) for each individual burst from TOPS file.
            #it is not exactly center time for bovl center, so it can be improved
            for sw in [1, 2,3]:
                TOPS_file=os.path.join(prime_path, f'{prime_epoch}.IW{sw}.slc.TOPS_par')
                if os.path.exists(TOPS_file):
                    for overlap in df['overlap'].unique():
                        sensing_start_time=slc_process(TOPS_file, f"sensing_start_time_{overlap}:")
                        if sensing_start_time:
                            sensing_start_time_str = str(timedelta(seconds=sensing_start_time))
                            df.loc[(df['epoch'].str.startswith(prime_epoch))  & (df['swath'] == str(sw)) & (df['overlap'] == overlap), 'center_time'] = sensing_start_time_str
                    
    # Ensuring all necessary columns are numeric
    df['ph_mean'] = pd.to_numeric(df['ph_mean'], errors='coerce')
    df['scale_factor'] = pd.to_numeric(df['scale_factor'], errors='coerce')
    df['az_pix'] = pd.to_numeric(df['az_pix'], errors='coerce')
    # Ensure all float columns are rounded to 8 decimal places
    df = df.apply(lambda x: x.round(8) if x.dtype.kind in 'fc' else x)
    # Set pandas display option to show floats with 8 decimal places
    pd.options.display.float_format = '{:.8f}'.format


    # Calculating esd_offset
    df['esd_it2'] = df['ph_mean'] * (-1) * df['scale_factor'] / df['az_pix']
    
    # Calculating daz_bovl based on daz_wrt
    df['daz_bovl'] = (df['icc_offset'] + df['esd_it1'] + df['esd_it2']) * df['az_pix'] * 1000
    df['daz_esd2_bovl']=df['esd_it2']* df['az_pix'] * 1000
    # Create 'swath_overlap' column and remove 'swath' and 'overlap' columns
    df['swath'] = df['swath'].astype(int)
    df['overlap'] = df['overlap'].astype(int)
    df = df.sort_values(by=['epoch', 'swath', 'overlap'])
    df.reset_index(drop=True, inplace=True)
    df['swath_overlap'] = df['swath'].astype(str) + '-' + df['overlap'].astype(str)
    df.drop(columns=['swath', 'overlap'], inplace=True)
    
    # Split 'epoch' into 'epoch' and 'master', and convert to datetime
    df['master'] = df['epoch'].apply(lambda x: x.split('_')[0])
    df['epoch'] = df['epoch'].apply(lambda x: x.split('_')[1])

    if daz_wrt:
        print('It will read from Lazeckys dataframe rather than calculating')
        dazML=pd.DataFrame()
        try:
            a=extract2txt_esds_frame(frame)
            dazML=dazML.append(a)
        except:
            print('frame '+frame+' is empty')
        df=df.merge(dazML[['epoch', 'daz_total_wrt_orbits']], on='epoch', how='left')
        df['daz_bovl']=(df['daz_total_wrt_orbits']+df['esd_it2'])* df['az_pix'] * 1000
        
    df['esd_master']=df.apply(lambda row: row['master'] if row['esd_master'] == '' else row['esd_master'], axis=1)
    df['epoch'] = pd.to_datetime(df['epoch'], format='%Y%m%d')
    df['master'] = pd.to_datetime(df['master'], format='%Y%m%d')
    df['esd_master'] = pd.to_datetime(df['esd_master'], format='%Y%m%d')
    df['daz_bovl'] = pd.to_numeric(df['daz_bovl'], errors='coerce')
    ##to see just proper decimals

    ####calculating the ctr_lat anc ctr_long for individual burst!
    ##check the bovls_rad2mm.py line 65-77 in licsar_proc. Kudos to Dr. Milan Lazecky
    bovljson = os.path.join(folder_path, frame+'.bovls.geojson')
    if not os.path.exists(bovljson):
        print('extracting burst polygons from LiCSInfo database')
        gpd_bursts = fc.frame2geopandas(frame, use_s1burst=True)
        try:
            gpd_bursts.to_file(bovljson, driver='GeoJSON')
        except:
            print('cannot write to LiCSAR_proc, saving locally')
            bovljson = frame+'.bovls.geojson'
            gpd_bursts.to_file(bovljson, driver='GeoJSON')
    #gpd_overlaps, overlap_gdf1, overlap_gdf2, overlap_gdf3 = extract_burst_overlaps(frame)
    gpd_overlaps, sw_overlaps_dict = fc.extract_burst_overlaps(frame, os.path.dirname(os.path.realpath(bovljson)))

    ###This is to try to move the geojson to procdir.
    homedir=os.getcwd()
    bovjson_file=os.path.join(homedir, frame + '.bovls.geojson')
    if os.path.exists(bovjson_file):
        try:
            # Try to move the file
            shutil.move(bovjson_file, folder_path)
            print(f"Moved {bovjson_file} to {folder_path}")
        except PermissionError:
            print(f"Permission denied: could not move {bovjson_file} to {folder_path}. geojson saved {os.getcwd()}.")
        except Exception as e:
            print(f"An error occurred while moving {bovjson_file} to {folder_path}: {e}")
    else:
        print(f"{bovjson_file} does not exist")

#    print('Coordinates of bovl centers saving to dataframe!')
    ######adding the ctr_lat and ctr_lon!
    gpd_swath_overlap(gpd_overlaps)
    print('1')
    add_centroid_columns(gpd_overlaps)
    print('2')
    df_cord=pd.merge(df, gpd_overlaps[['swath_overlap', 'center_lat', 'center_lon']], on='swath_overlap', how='left')
    # Ensure all float columns are rounded to 8 decimal places
    print('3')
    df_cord = df_cord.apply(lambda x: x.round(8) if x.dtype.kind in 'fc' else x)
    # Set pandas display option to show floats with 8 decimal places
    pd.options.display.float_format = '{:.8f}'.format
    df_cord['center_lat'] = df_cord['center_lat'].round(6)
    df_cord['center_lon'] = df_cord['center_lon'].round(6)
    print('4')
#######create the Heading and Incidence tif file, then obtain values in ctr_lat, lot for indv. burst
#########
    LiCS_public= os.environ.get('LiCSAR_public', './') 
    metadata_folder=os.path.join(LiCS_public,track_id,frame,'metadata')
    rng_azi_NEU_from_NEU(metadata_folder)
    # print('heading and incidence tif created, lets add them into dataframe!')

###
    df_centers = df_cord[['swath_overlap', 'center_lon', 'center_lat']].drop_duplicates(subset='swath_overlap')
    
    # Create a temporary directory
    temp_dir = os.path.join(batchdir, 'temp_for_daz_esd')
    os.makedirs(temp_dir, exist_ok=True)
    
    # Save the coordinates to a text file
    outdir = os.path.join(temp_dir, frame + '_sw_coords.txt')
    df_centers.to_csv(outdir, index=False, header=None, sep='\t')
    
    # Prepare the awk command
    coords_txt_path = os.path.join(temp_dir, frame + '_coords.txt')
    awk_command = f"awk '{{print $2, $3}}' {frame}_sw_coords.txt > {frame}_coords.txt"
    try:
        os.chdir(temp_dir)
        subprocess.run(awk_command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running awk command: {e}")
    
    # Run gmt grdtrack commands for the tif files
    for tifs in os.listdir(temp_dir):
        if tifs.endswith('angle.tif') and tifs.startswith(frame):
            tifs_path = os.path.join(temp_dir, tifs)
            if tifs.endswith('heading_angle.tif') and tifs.startswith(frame):
                command = f'gmt grdtrack -G{tifs_path} {coords_txt_path} > {frame}_heading.txt'
            else:
                command = f'gmt grdtrack -G{tifs_path} {coords_txt_path} > {frame}_incidence.txt'
            try:
                subprocess.run(command, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(f"An error occurred while running gmt grdtrack command: {e}")
    
    # Combine heading and inc for txt.
    final_command = f"paste {frame}_sw_coords.txt {frame}_heading.txt {frame}_incidence.txt | awk '{{print $1, $2, $3, $6, $9}}' > {frame}_lon_lan_head_inc.txt"
    try:
        subprocess.run(final_command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running final paste command: {e}")
    
    lot_lan_head_inc_path = os.path.join(temp_dir, f'{frame}_lon_lan_head_inc.txt')
    # Open the head_inc_txt as df to merge in the main dataframe.
    col_names = ['swath_overlap', 'center_lon', 'center_lat', 'heading', 'incidence']
    df_lot_lan_head_inc = pd.read_csv(lot_lan_head_inc_path, delim_whitespace=True, names=col_names)
    df_angle=pd.merge(df_cord, df_lot_lan_head_inc, on='swath_overlap')
    df_angle.drop(columns=['center_lat_x', 'center_lon_x'], inplace=True)
    df_angle.rename(columns={'center_lat_y': 'center_lat', 'center_lon_y': 'center_lon'}, inplace=True)
    
    ##flag S1AorB
    df_angle=flag_s1b_esds(df_angle)


    ##POD correction
    if orbdiff_fix:
        print('fixing the orb diff values by -39mm shifting')
        df_angle = fix_pod_offset(df_angle, using_orbits=False)        
    
    # Ensure all float columns are rounded to 8 decimal places
    df_angle = df_angle.apply(lambda x: x.round(8) if x.dtype.kind in 'fc' else x)
    # Set pandas display option to show floats with 8 decimal places
    pd.options.display.float_format = '{:.8f}'.format
    
    df_angle.drop(columns=['daz_total_wrt_orbits'], inplace=True)
    # Saving the DataFrame to CSV
    os.chdir(homedir)
    output_dir = os.path.join(batchdir, 'daz_esd')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    output_path = os.path.join(output_dir, frame + '.step1.csv')
    df_angle.to_csv(output_path, index=False)
    print(f"dataframe is ready for SET and Iono correction! it is saved to {output_path}")


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Usage as err:
        print("\nERROR:")
        print("  " + str(err.msg))
        print("\nFor help, use -h or --help.\n")
        sys.exit(2)