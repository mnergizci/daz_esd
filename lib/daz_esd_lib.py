import pygmt
import os,sys
import pandas as pd
import sys
import math
import numpy as np
import pickle
import matplotlib.pyplot as plt
from datetime import datetime
from datetime import timedelta
from sklearn.linear_model import LinearRegression
import warnings
from scipy.constants import speed_of_light
from lics_unwrap import *
import requests
from lxml import html
warnings.filterwarnings("ignore")
from osgeo import gdal
from astropy.timeseries import LombScargle
from scipy.optimize import curve_fit
import geopandas as gpd
gdal.UseExceptions()



def results_files_process(filename, df, orbit):
    """
    This function is to process the result files output to pandas dataframe. 
    In this code we split the files names concerning underscore.
    If the folder starts with 20's, it is the result file. 
    """
    try:
        for file in os.listdir(filename):
            result_find = file.split("_")
            if result_find[0].startswith("20"):
                with open(os.path.join(filename, file), "r") as a_file:
                    data = [line.strip().split() for line in a_file if len(line.strip().split())==11]
                new_rows = [{'epoch': file.replace('.results', ''), 'swath': line[0], 'overlap': line[1], 'ph_mean': line[2]} for line in data[1:]]
                df = df.append(pd.DataFrame(new_rows), ignore_index=True)
    except Exception as e:
        print(f"There is an Error: {e}")
    return df

def log_files_process(filename, df):
    for i, row in df.iterrows():
        coreg_file = "coreg_quality_" + row['epoch'] + '.log'
        try:
            with open(os.path.join(filename, coreg_file), "r") as a_file:
                lines = a_file.readlines()
                icc = sum(float(line.strip().split()[2]) for line in lines if "daz =" in line)
                icc = round(icc, 7)  # Round to 7 decimal places
                esd_it = sum(float(line.strip().split()[3]) for line in lines if line.strip().startswith("az_ovr iteration 1:"))
                esd_it = round(esd_it, 7)  # Round to 7 decimal places
                esd_mast_list= [line.strip().split()[-1] for line in lines if  'Spectral diversity estimation' in line]
                esd_mast = esd_mast_list[0] if esd_mast_list else ''
            df.at[i, 'icc_offset'] = icc
            df.at[i, 'esd_it1'] = esd_it
            df.at[i, 'esd_master'] = esd_mast
        except Exception as e:
            print(f"There is an Error: {e}")


"""
Until here, We exploit the log folder and use .results and .log files for our dataframe. 
Following code is relative to SLC folder and calculating f_dc.
"""

def process_file(filename, keyword):
    with open(filename, "r") as a_file:
        for line in a_file:
            stripped_line = line.strip()
            line_list = stripped_line.split()
            if keyword in line_list:
                return float(line_list[1])
    return None

def slc_process(filename, parameter_type):
    keywords = {
        "steering_angle": "az_steering_rate:",
        "heading_angle": "heading:",
        "center_lat": "center_latitude:",
        "center_lon": "center_longitude:",
        "center_range": "center_range_slc:",
        "center_time": "center_time:",
        "pixel_spacing": "azimuth_pixel_spacing:",
        "incidence_angle": "incidence_angle:",
        "near_range": "near_range_slc:",
        "far_range": "far_range_slc:"
    }
    
    if parameter_type.startswith("sensing_start_time_"):
        keyword = parameter_type  # Use the dynamic keyword directly
    else:
        keyword = keywords.get(parameter_type, "")
    
    return process_file(filename, keyword) 

def gpd_swath_overlap(df):
    '''
    first run, gpd_overlaps, sw_overlaps_dict = fc.extract_burst_overlaps(frame)
    afterwards this function create_swath_overlap columns to match the ctr_lat,lon with the main df.
    '''
    df['swath_overlap'] = None  # Initialize the new column

    # Group by 'swath_1' and calculate cumulative count for each group
    df['overlap_count'] = df.groupby('swath_1').cumcount() + 1

    # Create 'swath_overlap' column based on 'swath_1' and the cumulative count
    df['swath_overlap'] = df.apply(lambda row: f"{row['swath_1']}-{row['overlap_count']}", axis=1)

    # Drop the helper column
    df.drop(columns=['overlap_count'], inplace=True)

    return df
    
def add_centroid_columns(gpd):
    '''
    This function helps to create the ctr_lat, and ctr_lon for individual burst.
    
    '''
    # Calculate the centroid of each polygon
    gpd['centroid'] = gpd['geometry'].centroid
    
    # Extract the latitude and longitude of the centroid
    gpd['center_lat'] = gpd['centroid'].y
    gpd['center_lon'] = gpd['centroid'].x
    
    # Drop the temporary 'centroid' column if not needed
    gpd.drop(columns=['centroid'], inplace=True)
    
    return gpd
#######################
#####flag the dataset. 

# Define the flag_s1b function
def flag_s1b(epochdates, masterdate, mastersat='A', returnstr=False):
    """
    Args:
        epochdates (list of dt.datetime.date)
        masterdate (dt.datetime.date)
        mastersat (str): 'A' or 'B'
        returnstr (bool): if True, returns 'A' or 'B', otherwise returns 1 for 'B'
    """
    if mastersat == 'B':
        masterdate = masterdate + pd.Timedelta('6 days')
    isB = []
    for epoch in epochdates:
        # Ensure both dates are Timestamps for subtraction
        epoch = pd.Timestamp(epoch)
        masterdate = pd.Timestamp(masterdate)
        # ok, give +- 1 day tolerance due to midnight issue
        if np.abs(np.mod((epoch - masterdate).days, 12)) <= 1:
            if returnstr:
                val = 'A'
            else:
                val = 0
        else:
            if returnstr:
                val = 'B'
            else:
                val = 1
        isB.append(val)
    isB = np.array(isB)
    return isB

#The main function also take the frame.csv. God bless me I played the orig codes a lot.
def flag_s1b_esds(esds):
    esds['S1AorB'] = 'X'
    grouped = esds.groupby(['epoch', 'master'])
    for name, group in grouped:
        masterdate = pd.Timestamp(group.master.values[0])
        epochdates = group.epoch
        mastersat = 'A'  # Assuming 'A' for now, modify as needed
        group['S1AorB'] = flag_s1b(epochdates, masterdate, mastersat, returnstr=True)
        esds.update(group)
    return esds



########################
###open the esds and frame csv properly as dataframe

def df_preprepare_esds(esdsin, framespdin, firstdate = '', countlimit = 25):
    #basic fixes
    esds = esdsin.copy(deep=True)
    framespd = framespdin.copy(deep=True)
        # this helps for nans in 'master' (causing it float)
    framespd = framespd.dropna()
    esds['years_since_beginning'] = 0.0
    framespd['count_all'] = 0
    # framespd['daz_bovl_std_all'] = 0.0
    firstdatei = firstdate
    for swath_overlap, group in esds.groupby('swath_overlap'):
        if not firstdate:
            firstdatei = group['epoch'].min()
        swath_overlaps = framespd[framespd['swath_overlap'] == swath_overlap]
        if swath_overlaps.empty:
            #print('Warning, frame {} not found in framespd, using defaults'.format(frame))
            #azimuth_resolution = 14.0
            print('Warning, swath_overlap {} not found in framespd, skipping'.format(swath_overlap))
            esds = esds.drop(esds.loc[esds['swath_overlap']==swath_overlap].index)
            continue
        else:
            azimuth_resolution = float(swath_overlaps['az_pix'])
        count = group.epoch.count()
        if count < countlimit:
            print('small number of {} samples in swath '.format(str(count))+swath_overlap+' - removing')
    
            esds = esds.drop(esds.loc[esds['swath_overlap']==swath_overlap].index)
            framespd = framespd.drop(framespd.loc[framespd['swath_overlap']==swath_overlap].index)
            continue
        #remove median from daz_total
        # medianvalue = group['daz_total_wrt_orbits'].median()
        # group['daz_total_wrt_orbits'] = group['daz_total_wrt_orbits'] - medianvalue
        # #save the median correction values to framespd
        # framespd.at[frameta.index[0], 'daz_median_shift_mm'] = medianvalue*azimuth_resolution*1000
        framespd.at[swath_overlaps.index[0], 'count_all'] = int(count)
        # group['daz_mm'] = group['daz_total_wrt_orbits']*azimuth_resolution*1000
        # group['daz_cc_mm'] = group['daz_cc_wrt_orbits']*azimuth_resolution*1000
        group['years_since_beginning'] = group['epoch'] - firstdatei
        
        group['years_since_beginning'] = group['years_since_beginning'].apply(lambda x: float(x.days)/365.25)
        #get std, after detrending - but no need to save daz_detrended_mm now....
        # group['daz_detrended_mm'] = signal.detrend(group['daz_mm'])
        # framespd.at[frameta.index[0], 'daz_mm_std_all'] = np.std(group['daz_detrended_mm'])
        #update esds
        # # esds.update(group['daz_total_wrt_orbits'])
        # esds.update(group['daz_mm'])
        # esds.update(group['daz_cc_mm'])
        esds.update(group['years_since_beginning'])
    # extra check/fix?
    framespd=framespd.dropna()
    for subswat_overlap in framespd['swath_overlap']:
        if subswat_overlap not in esds['swath_overlap'].values:
            framespd = framespd.drop(framespd.loc[framespd['swath_overlap']==swath_overlap].index)
    # perhaps not necessary, but just in case..
    for swath_overlap in esds.swath_overlap.unique():
        if swath_overlap not in framespd['swath_overlap'].values:
            esds = esds.drop(esds.loc[esds['swath_overlap']==swath_overlap].index)
    #got those in mm, so we can remove the px values now
    # esds = esds.drop('daz_total_wrt_orbits', axis=1)
    # esds = esds.drop('daz_cc_wrt_orbits', axis=1)
    #we also do not need the RSLC3 information
    if 'esd_master' in esds.columns:
        esds = esds.drop('esd_master', axis=1)
    # some more to add - daz_tide:
    if 'daz_tide_mm' in esds.columns:
        esds['daz_bovl_notide'] = esds['daz_bovl'] - esds['daz_tide_mm']
        esds['daz_bovl_notide_noiono']=esds['daz_bovl'] - esds['daz_tide_mm'] - esds['daz_iono_mm']
    if 'daz_esd2_bovl' in esds.columns:
        print('daz inlcues esd2_bovl??!!2')
        esds['daz_bovl_esd2_notide_noino']=esds['daz_esd2_bovl']-esds['daz_tide_mm']-esds['daz_iono_mm']
    return esds, framespd



def load_csvs(esdscsv = 'esds.csv', framescsv = 'frames.csv', core_init = False):
    framespd = pd.read_csv(framescsv)
    esds = pd.read_csv(esdscsv)
    # Remove Unnamed columns
    esds = esds.loc[:, ~esds.columns.str.contains('^Unnamed')]
    framespd = framespd.loc[:, ~framespd.columns.str.contains('^Unnamed')]
    if 'version' in esds.columns:
        esds = esds.drop('version', axis=1)
    
    # Check and convert 'epoch' to datetime if it exists, else raise an error
    if 'epoch' in esds.columns:
        esds['epoch'] = esds['epoch'].apply(lambda x: pd.to_datetime(str(x)).date())
    else:
        raise ValueError("The column 'epoch' does not exist in the esds DataFrame.")

    if core_init:
        mindate = esds['epoch'].min()
        #maxdate = esds['epochdate'].max()
        esds, framespd = df_preprepare_esds(esds, framespd, mindate)
        #esds = esds.reset_index(drop=True)
        framespd = framespd.reset_index(drop=True)
    return esds, framespd
    
######################
####Fix orbit by 39 mm:
def fix_pod_offset(esds, using_orbits = False):
    """Function to fix the shift after new orbits in 2020-07-29/30, either using real POD diff if possible (if in LiCSAR), or applying 39 mm constant value.
    Args:
        esds (pd.Dataframe):   as loaded (i.e. with the relevant daz columns)
        using_orbits (bool):   if True, it will try use directly PODs to find diff (only with daz_lib_licsar)
    Returns:
        pd.DataFrame :  original esds with applied correction
    """
    col='daz_bovl'
    if not using_orbits:
        print('subtracting towards 2020-07-30')
        ep = esds[esds['epoch'] <= pd.Timestamp('2020-07-30')][col]
        offset_px = np.float64(39)  # just a mean value
        esds.loc[esds['epoch'] <= pd.Timestamp('2020-07-30'), col] = ep - offset_px
    return esds




###Codes from Dr. Qi Ou
#######################
class OpenTif:
    """ a Class that stores the band array and metadata of a Gtiff file.
    Can be initiated using a Tiff file alone, or with additional bands from
    its associated uncertainties, incidence and heading files. """
    def __init__(self, filename, sigfile=None, incidence=None, heading=None, N=None, E=None, U=None):
        self.ds = gdal.Open(filename)
        self.basename = os.path.splitext(os.path.basename(filename))[0]
        self.band = self.ds.GetRasterBand(1)
        self.data = self.band.ReadAsArray()
        self.xsize = self.ds.RasterXSize
        self.ysize = self.ds.RasterYSize
        self.left = self.ds.GetGeoTransform()[0]
        self.top = self.ds.GetGeoTransform()[3]
        self.xres = self.ds.GetGeoTransform()[1]
        self.yres = self.ds.GetGeoTransform()[5]
        self.right = self.left + self.xsize * self.xres
        self.bottom = self.top + self.ysize * self.yres
        self.projection = self.ds.GetProjection()
        pix_lin, pix_col = np.indices((self.ds.RasterYSize, self.ds.RasterXSize))
        self.lat, self.lon = self.top + self.yres*pix_lin, self.left+self.xres*pix_col

        # convert 0 and 255 to NaN
        #self.data[self.data == 0.00000000000000] = np.nan
        self.data[self.data == 255] = np.nan

        if sigfile is not None:
            self.dst = gdal.Open(sigfile)
            self.bandt = self.dst.GetRasterBand(1)
            self.sigma = self.bandt.ReadAsArray()
            self.sigma[self.sigma == 0] = np.nan
            if self.dst.RasterXSize != self.xsize or self.dst.RasterYSize != self.ysize:
                try:
                    self.sigma = self.sigma[:self.ysize, :self.xsize]
                except Warning:
                    print('Error: Sigma and Velocity file not the same size!')
                    print('sig has size = ' + str(self.dst.RasterXSize) + ', ' + str(self.dst.RasterYSize))
                    print('vel has size = ' + str(self.ds.RasterXSize) + ', ' + str(self.ds.RasterYSize))
                # self.sigma = np.ones((self.ysize, self.xsize))
        else:
            self.sigma = np.ones((self.ysize, self.xsize))

        if incidence is not None:
            self.ds_inc = gdal.Open(incidence)
            self.band_inc = self.ds_inc.GetRasterBand(1)
            self.inc = np.deg2rad(self.band_inc.ReadAsArray())
            self.inc[self.inc == 0] = np.nan
            if self.ds_inc.RasterXSize != self.xsize or self.ds_inc.RasterYSize != self.ysize:
                try:
                    self.inc = self.inc[:self.ysize, :self.xsize]
                except Warning:
                    print('Error: Inc and Velocity file not the same size!')
                    print('inc has size = ' + str(self.ds_inc.RasterXSize) + ', ' + str(self.ds_inc.RasterYSize))
                    print('vel has size = ' + str(self.ds.RasterXSize) + ', ' + str(self.ds.RasterYSize))

        if heading is not None:
            self.ds_head = gdal.Open(heading)
            self.band_head = self.ds_head.GetRasterBand(1)
            self.head = np.deg2rad(self.band_head.ReadAsArray())
            self.head[self.head == 0] = np.nan
            if self.ds_head.RasterXSize != self.xsize or self.ds_head.RasterYSize != self.ysize:
                try:
                    self.head = self.head[:self.ysize, :self.xsize]
                except Warning:
                    print('Error: Heading and Velocity file not the same size!')
                    print('head has size = ' + str(self.ds_head.RasterXSize) + ', ' + str(self.ds_head.RasterYSize))
                    print('vel has size = ' + str(self.ds.RasterXSize) + ', ' + str(self.ds.RasterYSize))

        if N is not None:
            self.ds_N = gdal.Open(N)
            self.band_N = self.ds_N.GetRasterBand(1)
            self.N = self.band_N.ReadAsArray()
            # self.N[self.N == 0] = np.nan
            if self.ds_N.RasterXSize != self.xsize or self.ds_N.RasterYSize != self.ysize:
                try:
                    self.N = self.N[:self.ysize, :self.xsize]
                except Warning:
                    print('Error: Heading and Velocity file not the same size!')
                    print('head has size = ' + str(self.ds_N.RasterXSize) + ', ' + str(self.ds_N.RasterYSize))
                    print('vel has size = ' + str(self.ds.RasterXSize) + ', ' + str(self.ds.RasterYSize))

        if E is not None:
            self.ds_E = gdal.Open(E)
            self.band_E = self.ds_E.GetRasterBand(1)
            self.E = self.band_E.ReadAsArray()
            # self.E[self.E == 0] = np.nan
            if self.ds_E.RasterXSize != self.xsize or self.ds_E.RasterYSize != self.ysize:
                try:
                    self.E = self.E[:self.ysize, :self.xsize]
                except Warning:
                    print('Error: Heading and Velocity file not the same size!')
                    print('head has size = ' + str(self.ds_E.RasterXSize) + ', ' + str(self.ds_E.RasterYSize))
                    print('vel has size = ' + str(self.ds.RasterXSize) + ', ' + str(self.ds.RasterYSize))

        if U is not None:
            self.ds_U = gdal.Open(U)
            self.band_U = self.ds_U.GetRasterBand(1)
            self.U = self.band_U.ReadAsArray()
            # self.U[self.U == 0] = np.nan
            if self.ds_U.RasterXSize != self.xsize or self.ds_U.RasterYSize != self.ysize:
                try:
                    self.U = self.U[:self.ysize, :self.xsize]
                except Warning:
                    print('Error: Heading and Velocity file not the same size!')
                    print('head has size = ' + str(self.ds_U.RasterXSize) + ', ' + str(self.ds_U.RasterYSize))
                    print('vel has size = ' + str(self.ds.RasterXSize) + ', ' + str(self.ds.RasterYSize))

def export_tif(data, df, filename):
    # Export data to tif format.
    driver = gdal.GetDriverByName("GTiff")
    outdata = driver.Create(filename, df.xsize, df.ysize, 1, gdal.GDT_Float32)
    outdata.SetGeoTransform([df.left, df.xres, 0, df.top, 0, df.yres])  ##sets same geotransform as input
    # outdata.SetProjection(df.projection)  ##sets same projection as input
    outdata.GetRasterBand(1).WriteArray(data)
    outdata.FlushCache()
    outdata.FlushCache()  # need to flush twice to export the last tif properly, otherwise it stops halfway.



def rng_azi_NEU_from_NEU(path):
    '''
    This code is adapted from Dr. Qi Ou.
    '''
    frame_dir=os.path.dirname(path)
    frame = os.path.basename(frame_dir)
    orientation = frame[3]
    track = frame[:4]

    n_file = os.path.join(path, "{}.geo.N.tif".format(frame))
    e_file = os.path.join(path, "{}.geo.E.tif".format(frame))
    u_file = os.path.join(path, "{}.geo.U.tif".format(frame))
    n = OpenTif(n_file)
    e = OpenTif(e_file)
    u = OpenTif(u_file)

    if orientation == "A":
        inc_rad = np.arccos(u.data)
        head_rad = np.arcsin(n.data / np.sin(inc_rad))
        heading=np.degrees(head_rad)
        incidence=np.degrees(inc_rad)
    elif orientation == "D":
        inc_rad = np.arccos(u.data)
        head_rad = np.arcsin(- n.data / np.sin(inc_rad)) - np.pi
        heading=np.degrees(head_rad)
        incidence=np.degrees(inc_rad)
    else:
        raise ValueError("The 4th character of frameID is neither A nor D, please check your frame name.")

    # Create the output directory if it doesn't exist
    batchdir = os.environ.get('BATCH_CACHE_DIR', './')
    temp_dir=os.path.join(batchdir, 'temp_for_daz_esd')
    
    os.makedirs(temp_dir, exist_ok=True)

    export_tif(heading, n, os.path.join(temp_dir, "{}_heading_angle.tif".format(frame)))
    export_tif(incidence, n, os.path.join(temp_dir, "{}_incidence_angle.tif".format(frame)))
#################
#################dfDC calculation from LiCS_daz library and ESD study on TR Eq. I can make here much better. Please improve the azimuth pixel size variable regarding in individual burst.



def get_param_gamma(param, parfile, floatt = True, pos = 0):
    a = grep1line(param,parfile).split()[1+pos]
    if floatt:
        a = float(a)
    return a

def s1_azfm(r, t0, azp):
  """azfr = s1_azfm(r, t0, azp)
  
  Calculate azimuth FM rate given slant range, reference slant-range delay and the azimuth FM rate polynomial for ScanSAR data
  
  **Arguments:**
  
  * r:    slant range (meters)
  * t0:   reference slant range time for the polynomial (center swath delay in s)
  * azp:  polynomial coefficients
  
  **Output:**
  
  * the function returns the azimuth FM rate"""

  tsr = 2.0 * r / speed_of_light;
  dt = tsr - t0;
  azfr = azp[0] + dt * (azp[1] + dt*(azp[2] + dt*(azp[3] + dt*azp[4])));
  return azfr;
    
def get_dfDC_slc(path_to_slcdir, f0=5405000500, burst_interval = 2.758277, returnka = True, returnperswath = False, returnscalefactor=True):
    #f0 = get_param_gamma('radar_frequency', parfile)
    #burst_interval = get_param_gamma('burst_interval', topsparfile)
    epoch = os.path.basename(path_to_slcdir)
    print(epoch)
    frame = path_to_slcdir.split('/')[-3]
    pi=np.pi
    
    if len(frame)!=17:
        frame=path_to_slcdir.split('/')[-4]
    parfile = os.path.join(path_to_slcdir, epoch+'.slc.par')
    #parfile = glob.glob(path_to_slcdir+'/????????.slc.par')[0]
    #topsparfiles = glob.glob(path_to_slcdir+'/????????.IW?.slc.TOPS_par')
    #iwparfiles = glob.glob(path_to_slcdir+'/????????.IW?.slc.par')
    #
    
    lam = speed_of_light / f0
    dfDC = []
    kas = []
    ctr_range = []
    far_range = []
    near_range = []
    scalefactor= []
    afmrate_srdly = []
    afmrate_ply= []
    kr_list = []
    numbursts = [ int(frame.split('_')[2][:2]), int(frame.split('_')[2][2:4]), int(frame.split('_')[2][4:6])]
    azps_list = []
    az_line_time_list = []
    #krs = []
    #print('This is a proper solution but applied to primary SLC image. originally it is applied by GAMMA on the RSLC...')
    #for n in range(len(topsparfiles)):
    for n in [1,2,3]:
        topsparfile = os.path.join(path_to_slcdir, epoch+'.IW'+str(n)+'.slc.TOPS_par')
        iwparfile = os.path.join(path_to_slcdir, epoch+'.IW'+str(n)+'.slc.par')
        if (not os.path.exists(iwparfile)) or (not os.path.exists(topsparfile)):
            dfDC.append(np.nan)
            kas.append(np.nan)
            numbursts[n-1] = np.nan
        else:
            #topsparfile = topsparfiles[n]
            #iwparfile = iwparfiles[n]
            az_steering_rate = get_param_gamma('az_steering_rate', topsparfile) # az_steering_rate is the antenna beam steering rate
            #az_ster_rate.append(az_steering_rate)
            r1 = get_param_gamma('center_range_slc', iwparfile)
            #get the satellite velocity
            midNstate = int(get_param_gamma('number_of_state_vectors', iwparfile)/2)+1
            sv = 'state_vector_velocity_' + str(midNstate)
            velc1 = get_param_gamma(sv, iwparfile, pos=0)
            velc2 = get_param_gamma(sv, iwparfile, pos=1)
            velc3 = get_param_gamma(sv, iwparfile, pos=2)
            vsat = np.sqrt(velc1**2 + velc2**2 + velc3**2)
            midNstate=1
            # now some calculations
            afmrate_srdelay = get_param_gamma('az_fmrate_srdelay_'+ str(midNstate), topsparfile)
            afmrate_srdly.append(afmrate_srdelay) 
            afmrate_poly = []
            afmrate_poly.append(get_param_gamma('az_fmrate_polynomial_' + str(midNstate), topsparfile, pos = 0))
            afmrate_poly.append(get_param_gamma('az_fmrate_polynomial_' + str(midNstate), topsparfile, pos = 1))
            afmrate_poly.append(get_param_gamma('az_fmrate_polynomial_' + str(midNstate), topsparfile, pos = 2))
            try:
                afmrate_poly.append(get_param_gamma('az_fmrate_polynomial_' + str(midNstate), topsparfile, pos = 3))
            except:
                afmrate_poly.append(0)
            try:
                afmrate_poly.append(get_param_gamma('az_fmrate_polynomial_' + str(midNstate), topsparfile, pos = 4))
            except:
                afmrate_poly.append(0)
            afmrate_ply.append(afmrate_poly)
            ka = s1_azfm(r1, afmrate_srdelay, afmrate_poly) #unit: Hz/s == 1/s^2
            kr = -2.0 * vsat * az_steering_rate*(pi / 180.0) / lam
            kr_list.append(kr)
            if (kr != 0.0):
                #kt = ka * kr/(kr - ka)
                # but maybe should be kt = (kr*ka)/(ka-kr) # see https://iopscience.iop.org/article/10.1088/1755-1315/57/1/012019/pdf  --- and remotesensing-12-01189-v2, and Fattahi et al...
                # ok, gamma reads kr to be ka... corrected
                kt = kr * ka/(ka - kr)
            else:
                kt = -ka
            #finally calculate dfDC:
            #burst_interval = get_param_gamma('burst_interval', topsparfile)
            kas.append(ka)
            #krs.append(kr)
            dfDC.append(kt*burst_interval) #burst_interval is time within the burst... we can also just calculate.. see Grandin: eq 15: hal.archives-ouvertes.fr/hal-01621519/document
            #ok, that's the thing - burst_interval is actually t(n+1) - t(n) - see remotesensing-12-01189-v2
            #so it should be kt * -burst_interval, that is why GAMMA has the -kt J ... ok, good to realise this

            ####calculating scalefactor (Nergizci)
            azps=np.float64(get_param_gamma('azimuth_pixel_spacing', iwparfile))
            azps_list.append(azps)
            az_line_time=np.float64(get_param_gamma('azimuth_line_time', iwparfile))
            az_line_time_list.append(az_line_time)
        
            dfdc=kt*burst_interval
            sf=(azps)/(dfdc*az_line_time*2*np.pi)
            scalefactor.append(sf)
    return scalefactor, azps_list, kas, dfDC

##########for_solid_earth_tide from Lazecky############
def EN2azi(N, E, heading = -169):
    #ascending: around -13
    #descending: around -169 (or around 13, doesn't matter)
    alpha = np.deg2rad(heading)
    return E*np.sin(alpha)+N*np.cos(alpha)


def merge_tides(esds, earthtides):
    esds['daz_tide_mm'] = 0.0

    for i, row in esds.iterrows():
        epoch = row['epoch']
        lat = row['center_lat']
        lon = row['center_lon']
        
        # Find matching rows in earthtides
        tiderow = earthtides[
            (earthtides['epoch'] == epoch) & 
            (earthtides['lat'].round(6) == round(lat, 6)) & 
            (earthtides['lon'].round(6) == round(lon, 6))
        ]
        
        if tiderow.empty:
            print(f'Error - no matching entry in tides for epoch {epoch}, lat {lat}, lon {lon}')
            continue

        E = tiderow['E'].values[0]
        N = tiderow['N'].values[0]
        heading = row['heading']
        daz_tide_mm = round(EN2azi(N, E, heading) * 1000, 6)
        esds.at[i, 'daz_tide_mm'] = daz_tide_mm
    esds['daz_bovl_notide']=esds['daz_bovl']-esds['daz_tide_mm']
    esds['daz_esd2_bovl_notide']=esds['daz_esd2_bovl']-esds['daz_tide_mm']
    
    print('\nDone')
    return esds

#############FOR PPM:
#################################For PPM!
def EN2azi(N, E, heading = -169):
    #ascending: around -13
    #descending: around -169 (or around 13, doesn't matter)
    alpha = np.deg2rad(heading)
    #thanks Chris Rollins!!!!
    return E*np.sin(alpha)+N*np.cos(alpha)

def get_ITRF_ENU(lat, lon, model='itrf2014', refto='NNR'):
    '''Gets plate motion model values from UNAVCO website
    
    Args
        lat (float): latitude
        lon (float): longitude
        model (string): choose model code, e.g. 'gsrm_2014', 'itrf2014',.. see lines after l. 571 of `view-source:https://www.unavco.org/software/geodetic-utilities/plate-motion-calculator/plate-motion-calculator.html`
        refto (string): choose reference plate, e.g. 'NNR', 'EU',... see lines on the link above after line 683
    
    Returns
        float, float: east and north component of plate motion in given coordinates
    '''
    url = "https://www.unavco.org/software/geodetic-utilities/plate-motion-calculator/plate-motion/model"
    # data to be sent to api
    data = {'name':'modelform',
        'id':'modelform',
        'lat':str(lat),
        'lon':str(lon),
        'model':model,
        'reference':refto,
        'format':'ascii'}
    try:
        # Sending post request and saving response as response object
        r = requests.post(url, data=data)
        r.raise_for_status()  # Raise an error for bad status codes
        #print(r.text)  # Debug: Print the response content to see what's returned
        # Parse the response content
        cont = html.fromstring(r.content)
        pre = cont.xpath('//pre')[0].text_content()
        #print(pre)  # Debug: Print the pre tag content to see the extracted text
        values = pre.split()
        E = float(values[2])
        N = float(values[3])
        return E, N
    except requests.exceptions.Timeout:
        print("Connection error: Request timed out")
        raise
    except requests.exceptions.TooManyRedirects:
        print("Connection error: Too many redirects")
        raise
    except requests.exceptions.RequestException as e:
        print(f"Connection error: {e}")
        raise
    except (ValueError, IndexError) as e:
        print(f"Parsing error: {e}")
        raise



# get ITRF N, E values
def get_itrf_gps_EN(df, samplepoints=3, velnc='vel_gps_kreemer.nc', refto='NNR', rowname = 'centroid'):
    '''Gets EN velocities from ITRF2014 plate motion model (auto-extract from UNAVCO website)
    In case velnc exists, it will be used as well, to generate GPS_N/E.. 
    I prepared the vel_gps_kreemer.nc file from data available in supplementary files of article DOI:10.1002/2014GC005407
    i.e. from ggge20572-sup-0015-suppinfofig14.Z : vel_1deg_NNR.gmt
    Basically, I used the existing velocities and correctly georeferenced them to WGS-84 frame, stored as 1 deg grid in NetCDF
    '''
    usevel = False
    if velnc:
        if os.path.exists(velnc):
            print('found velocities nc file - using it instead of ITRF2014')
            vels=xr.open_dataset(velnc)
            usevel = True
    itrfs_N = []
    itrfs_E = []
    itrfs_rms_N = []
    itrfs_rms_E = []
    if usevel:
        GPS_N = []
        GPS_E = []
        GPS_rms_N = []
        GPS_rms_E = []
    iii = 0
    fullcount = len(df)
    for ind, row in df.iterrows():
        iii = iii+1
        print('getting ITRF for {0}/{1} cells'.format(iii, fullcount))
        clon = row[rowname+'_lon']
        clat = row[rowname+'_lat']
        print('extracting PMM values from ITRF2014')
        # use a median over 'whole' frame:
        itrfEs = []
        itrfNs = []
        leng=round(clon*10+23.4/2)+1-round(clon*10-23.4/2)
        for i in range(round(clon*10-23.4/2),round(clon*10+23.4/2)+1,int(leng/samplepoints)):
            lon = i/10
            for j in range(round(clat*10-23.4/2),round(clat*10+23.4/2)+1,int(leng/samplepoints)):
                lat = j/10
                
                try:
                    # print(clon, clat)
                    # print(lon, lat)
                    itrfE, itrfN = get_ITRF_ENU(lat, lon, refto=refto)
                    #print(itrfE, itrfN)
                    #itrfE, itrfN = 0,0 #debug
                    itrfEs.append(itrfE)
                    itrfNs.append(itrfN)
                    #itrfs.append(EN2azi(N, E, heading))
                except:
                    print('connection error')
        if usevel:
            print('extracting PMM values from GPS stations')
            gps1_E = vels.sel(lon=slice(clon-125/111, clon+125/111), lat=slice(clat-125/111, clat+125/111)).VEL_E #.values
            gps1_N = vels.sel(lon=slice(clon-125/111, clon+125/111), lat=slice(clat-125/111, clat+125/111)).VEL_N #.values
            GPS_E.append(float(gps1_E.mean()))
            GPS_N.append(float(gps1_N.mean()))
            GPS_rms_E.append(float(gps1_E.std(ddof=1)))
            GPS_rms_N.append(float(gps1_N.std(ddof=1)))
        itrfs_E.append(np.mean(itrfEs))
        itrfs_N.append(np.mean(itrfNs))
        itrfs_rms_E.append(np.std(itrfEs,ddof=1))
        itrfs_rms_N.append(np.std(itrfNs,ddof=1))
    df['ITRF_N'] = itrfs_N
    df['ITRF_E'] = itrfs_E
    df['ITRF_RMSE_E'] = itrfs_rms_E
    df['ITRF_RMSE_N'] = itrfs_rms_N
    if usevel:
        df['GPS_N'] = GPS_N
        df['GPS_E'] = GPS_E
        df['GPS_RMSE_E'] = GPS_rms_E
        df['GPS_RMSE_N'] = GPS_rms_N
    return df



def df_get_itrf_gps_slopes(framespd, velnc='vel_gps_kreemer.nc', add_eu = True):
    '''
    will get both itrf 2014 pmm and gps stations averages, converted to LOS of frames in framespd
    '''
    framespd = get_itrf_gps_EN(framespd, samplepoints=3, velnc=velnc, refto='NNR', rowname = 'center')
    framespd['slope_plates_vel_azi_itrf2014'] = 0.0
    #framespd['slope_plates_vel_azi_itrf2014_point'] = 0.0
    if 'GPS_N' in framespd.columns:
        framespd['slope_plates_vel_azi_gps'] = 0.0
    print('converting to LOS')
    for ind,frameta in framespd.iterrows():
        #frame = frameta['frame']
        heading = frameta['heading']
        N = frameta['ITRF_N']
        E = frameta['ITRF_E']
        framespd.at[ind, 'slope_plates_vel_azi_itrf2014'] = EN2azi(N, E, heading)
        if 'GPS_N' in framespd.columns:
            N = frameta['GPS_N']
            E = frameta['GPS_E']
            framespd.at[ind, 'slope_plates_vel_azi_gps'] = EN2azi(N, E, heading)
    if add_eu:
        print('repeating for Eurasia-referenced framework')
        framespd = framespd.rename({'ITRF_N': 'ITRF_N_NNR', 'ITRF_E': 'ITRF_E_NNR'})
        framespd = get_itrf_gps_EN(framespd, samplepoints=3, velnc=None, refto='EU', rowname='center')
        for ind, frameta in framespd.iterrows():
            # frame = frameta['frame']
            heading = frameta['heading']
            N = frameta['ITRF_N']
            E = frameta['ITRF_E']
            framespd.at[ind, 'slope_plates_vel_azi_itrf2014_eur'] = EN2azi(N, E, heading)
        framespd = framespd.rename({'ITRF_N':'ITRF_N_EUR', 'ITRF_E':'ITRF_E_EUR'})
    return framespd

def seasonal_model(x, A, omega, phi, C):
    """Sinusoidal seasonal model: A*sin(omega*x + phi) + C"""
    return A * np.sin(omega * x + phi) + C

######Here is for illustration functions!
def compute_trends(x, y, use_curve_fit=False):
    """Compute trends using either linear regression or sinusoidal curve fitting."""
    x = pd.to_numeric(x, errors='coerce').astype(float)
    y = pd.to_numeric(y, errors='coerce').astype(float)

    if y.isna().sum() > 0.5 * len(y):  # Ignore if more than 50% NaN values
        return None, None

    x = np.array(x).reshape(-1, 1)  # Reshape for sklearn
    y = np.array(y)

    if use_curve_fit:
        # Initial guess for seasonal parameters: A (amplitude), omega (2pi/period), phi (phase shift), C (offset)
        guess = [np.std(y), 2 * np.pi / 365, 0, np.mean(y)]
        try:
            params, _ = curve_fit(seasonal_model, x.flatten(), y, p0=guess, maxfev=10000)
            y_fit = seasonal_model(x.flatten(), *params)
        except RuntimeError:
            y_fit = None  # In case curve fitting fails
    else:
        # Linear Trend (Linear Regression)
        linreg = LinearRegression().fit(x, y)
        y_fit = linreg.predict(x)

    # Seasonal Trend (Lomb-Scargle)
    frequency, power = LombScargle(x.flatten(), y).autopower(minimum_frequency=1/365, maximum_frequency=1/10)
    seasonal_period = 1 / frequency[np.argmax(power)] if len(frequency) > 0 else None  # Period in days

    return y_fit, seasonal_period



def daz_plot_esd(df, frame='', epoch_col='epoch', daz_orig=None, daz_ion=None, daz_set=None, output_dir=None):
    """
    Plots the DAZ values over time with both seasonal and linear trends.

    - Uses linear regression for daz_orig.
    - Uses sinusoidal fitting (curve_fit) for daz_iono and daz_tide.
    - Displays seasonal signal period (in days) in the legend.

    Parameters:
    df (pd.DataFrame): DataFrame containing the data.
    frame (str): Frame identifier for saving output images.
    epoch_col (str): Name of the column containing epoch dates.
    daz_orig (str): Name of the column containing original DAZ values.
    daz_ion (str): Name of the column containing ionospheric corrected DAZ values.
    daz_set (str): Name of the column containing SET corrected DAZ values.
    output_dir (str): Directory to save the plot. If None, the plot will not be saved.
    """

    # Convert 'epoch' to datetime if it's not already
    df[epoch_col] = pd.to_datetime(df[epoch_col])
    df["days_since_start"] = (df[epoch_col] - df[epoch_col].min()).dt.days  # Convert time to numerical format

    # Ensure the columns exist
    valid_cols = [col for col in [daz_orig, daz_ion, daz_set] if col in df.columns]
    if not valid_cols:
        print("Error: None of the specified DAZ columns exist in the dataset!")
        print("Available columns:", df.columns)
        return

    # Create the plot
    fig, ax = plt.subplots(figsize=(10, 6))

    colors = {"daz_bovl": "blue", "daz_iono_mm": "orange", "daz_tide_mm": "green"}
    labels = {"daz_bovl": "Original DAZ", "daz_iono_mm": "CODE-GIM", "daz_tide_mm": "SET"}

    for col, color in zip([daz_orig, daz_ion, daz_set], colors.values()):
        if col and col in df:
            x = df[epoch_col]
            y = df[col]

            # Use curve fitting for daz_iono and daz_tide, but linear fit for daz_orig
            use_curve_fit = col in [daz_ion, daz_set]
            y_fit, seasonal_period = compute_trends(df["days_since_start"], df[col], use_curve_fit=use_curve_fit)

            # Scatter plot for raw data
            ax.scatter(x, y, label=f"{labels[col]} (Seasonal: {seasonal_period:.1f} days)" if seasonal_period else labels[col], color=color, alpha=0.5)

            # Plot fitted trend
            if y_fit is not None:
                ax.plot(x, y_fit, linestyle="--", color=color, label=f"{labels[col]} Trend (Curve Fit)" if use_curve_fit else f"{labels[col]} Trend (Linear)")

    # Formatting the plot
    ax.set_xlabel("Epoch")
    ax.set_ylabel("DAZ [mm]")
    ax.set_title("DAZ Trends with Seasonal and Linear Components")
    ax.legend()
    ax.grid(True)

    # Get the minimum and maximum dates from x
    min_date = df[epoch_col].min().date()
    max_date = df[epoch_col].max().date()

    # Set the x-axis tick locations and labels
    xticks = pd.date_range(start=min_date, end=max_date, freq='YS')  # Year start frequency
    xticks_labels = [date.strftime('%Y') for date in xticks]
    plt.xticks(xticks, xticks_labels, rotation='vertical')

    # Optionally, save the plot
    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)
        filename = os.path.join(output_dir, f'{frame}.jpg')
        plt.savefig(filename, dpi=300)
        print(f"Plot saved to {filename}")

    plt.show()
    

def daz_plotting_point_pygmt(asc_gnss_file, dem_file, esd_step5_file, frame_step5_file, data_dir='daz_esd', region=[65, 95, 30, 50], dem_resolution='15s'):
    '''
    # Usage example
    data_dir = 'daz_esd'
    asc_gnss_file = '1_tien_shan_az_asc.txt'
    dem_file = 'earth_relief_tien_shan_15s.nc'
    esd_step5_file = '027A_04790_141413.step5.csv'
    frame_step5_file = '027A_04790_141413.frame.step5.csv'
    region = [65, 95, 30, 50]
    plot_deformation_rate(asc_gnss_file, dem_file, esd_step5_file, frame_step5_file, data_dir, region)
    '''
    
    
    # Initialize the PyGMT figure
    fig = pygmt.Figure()

    # Folder path
    batchdir = os.environ.get('BATCH_CACHE_DIR', './')
    batchdir = '/work/scratch-pw3/licsar/mnergiz/batchdir'
    datadir = os.path.join(batchdir, data_dir)
    
    # Upload data
    asc_GNSS = pd.read_csv(os.path.join(datadir, asc_gnss_file), delimiter=' ', header=None, names=['lon', 'lat', 'vel'])
    dem = os.path.join(batchdir, dem_file)
    
    # DEM downloading
    if not os.path.exists(dem):
        print('DEM is downloading please wait! After downloading, the process will be faster!')
        try:
            # Download the earth relief data and save it to a file
            grid = pygmt.datasets.load_earth_relief(resolution=dem_resolution, region=region)
            # Saving the grid to a NetCDF file
            grid.to_netcdf(dem)
            print(f"Data successfully downloaded and saved to {dem}")
        except Exception as e:
            print(f"An error occurred: {e}")
    else:
        print(f'DEM already exists!')

    # Load ESD and frame data
    esd = pd.read_csv(os.path.join(datadir, esd_step5_file))
    framesd = pd.read_csv(os.path.join(datadir, frame_step5_file))

    # Define the region (bounding box) for the plot
    plot_region = [
        esd['center_lon'].min() - 1, esd['center_lon'].max() + 1,
        esd['center_lat'].min() - 1, esd['center_lat'].max() + 1
    ]

    # Plot the basemap with DEM relief as the background
    fig.basemap(region=plot_region, projection="M13c", frame=["af", '+t"Deformation Rate"'])

    # Add DEM relief background using the locally saved file
    pygmt.makecpt(cmap="gray", series=[-8000, 8000, 1000], continuous=True)
    fig.grdimage(
        grid=dem,
        cmap=True,
        region=plot_region,
        projection="M13c",
        shading=True,
        frame=True
    )

    # Create a colormap for daz_bovl using the "roma" colormap with the specified range
    pygmt.makecpt(cmap="roma", series=[-5, 5])

    # Overlay the data points as circles with colors based on daz_bovl
    fig.plot(
        x=framesd['center_lon'], y=framesd['center_lat'], 
        style="c0.5c",  # circle with diameter of 0.5 cm
        color=framesd['slope_daz_bovl_esd2_notide_noino_mmyear'], 
        cmap=True,
        pen="black"
    )

    # Asceding GNSS data
    fig.plot(
        x=asc_GNSS['lon'], y=asc_GNSS['lat'], 
        style="t0.5c",  # triangle with diameter of 0.5 cm
        color=asc_GNSS['vel'], 
        cmap=True,
        pen="black"
    )

    # Add a color bar for daz_bovl values, with intervals of 5 mm
    fig.colorbar(frame='a5+l"es2"')

    # Show the plot
    fig.show()

