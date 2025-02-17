#!/usr/bin/env python3

from daz_lib import *


#optional:
try:
    import ephem
except:
    print('WARNING: no ephem library, the [optional] dusk/dawn time cannot be calculated')

#for iono correction
import nvector as nv
import iri2016
import pyproj
import numpy as np
import re
import glob
try:
    import wget
    #import zlib
    from LiCSAR_misc import grep1line
except:
    print('error loading libraries to use of CODE (either wget or LiCSAR_misc). Use of CODE will fail')


# get daz iono
################### IONOSPHERE 

def extract_iono_full(esds, framespd, ionosource = 'code', use_iri_hei=False):
    """ Full extraction of ionospheric effect from ionosource.
    Note this will create column with the phase advanced effect recalculated to apparent azimuth offset [mm] that has opposite sign.
    Therefore this conforms the GRL article and you can subtract this correction from the original values, as usual.
    
    Args:
        ionosource (str):   either 'iri' or 'code'
        use_iri_hei (bool): estimating F2 peak altitude using IRI2016. CODE is valid for 450 km, so might be better to set this OFF for using CODE
    Returns:
        esds, framespd
    """
    # estimating the ionosphere - takes long (several hours)
    # also, storing TECS values (i.e. TEC in slant direction, from IRI2016)
    esds['tecs_A'] = 0.0
    esds['tecs_B'] = 0.0
    esds['daz_iono_mm'] = 0.0
    #esds['daz_mm_notide_noiono_grad'] = 0.0
    # framespd['Hiono'] = 0.0
    # framespd['Hiono_std'] = 0.0
    # framespd['Hiono_range'] = 0.0
    framespd['tecs_A'] = 0.0
    framespd['tecs_B'] = 0.0
    for bovl in framespd['swath_overlap']:
        print(bovl)
        resolution = framespd[framespd['swath_overlap'] == bovl]['az_pix'].values[0] # in metres
        try:
            #daz_iono_with_F2 = calculate_daz_iono(frame, esds, framespd)
            #daz_iono_grad, hionos, tecs_A_master, tecs_B_master = calculate_daz_iono(frame, esds, framespd, method = 'gomba', out_hionos = True, out_tec_master = True)
            if use_iri_hei:
                daz_iono_grad, hionos, tecs_A_master, tecs_B_master, tecs_A, tecs_B = calculate_daz_iono(bovl, esds, framespd, method = 'gradient', out_hionos = True, out_tec_all = True, ionosource=ionosource, use_iri_hei=use_iri_hei)
                hiono = np.mean(hionos)
                hiono_std = np.std(hionos)
                
            else:
                
                daz_iono_grad, tecs_A_master, tecs_B_master, tecs_A, tecs_B = calculate_daz_iono(bovl, esds,
                                                                                                         framespd,
                                                                                                         method='gradient',
                                                                                                         out_hionos=False,
                                                                                                         out_tec_all=True,
                                                                                                         ionosource=ionosource,
                                                                                                         use_iri_hei=use_iri_hei)

                hiono = 450
                hiono_std = 0
        except:
            print('some error occurred extracting TEC(s) here')
            continue
        selesds=esds[esds['swath_overlap']==bovl].copy()
        # 2023/08: changing sign to keep consistent with the GRL article
        selesds['daz_iono_mm'] = daz_iono_grad*resolution*1000
        selesds['tecs_A'] = tecs_A
        selesds['tecs_B'] = tecs_B
        # skipping the correction here, since daz_mm_notide might not exist/not needed:
        #selesds['daz_mm_notide_noiono_grad'] = selesds['daz_mm_notide'] + selesds['daz_iono_grad_mm'] #*resolution*1000
        esds.update(selesds)
        framespd.at[framespd[framespd['swath_overlap'] == bovl].index[0], 'Hiono'] = hiono
        if use_iri_hei:
            framespd.at[framespd[framespd['swath_overlap']==bovl].index[0], 'Hiono_std'] = hiono_std
            framespd.at[framespd[framespd['swath_overlap']==bovl].index[0], 'Hiono_range'] = max(hionos)-min(hionos)
        framespd.at[framespd[framespd['swath_overlap']==bovl].index[0], 'tecs_A'] = tecs_A_master
        framespd.at[framespd[framespd['swath_overlap']==bovl].index[0], 'tecs_B'] = tecs_B_master
        #esds.at[esds[esds['frame']==frame].index, 'daz_mm_notide_noiono_F2'] = esds[esds['frame']==frame]['daz_mm_notide'] - esds['daz_iono_with_F2']*resolution*1000
    return esds, framespd


#######################################
# step 3 - get daz iono
################### IONOSPHERE 

def get_tecs(glat, glon, altitude, acq_times, returnhei = False, source='iri', alpha = 0.85):
    '''Gets estimated TEC over given point, up to given altitude
    
    Args:
        glat, glon, altitude: coordinates and max 'iono height' to get the TEC values for. Altitude is in km
        acq_times (list of dt.datetime): time stamps to get the TEC over the given point
        returnhei (boolean):  if True, it would return TEC values but also estimated F2 peak heights (from IRI)
        source (str): source of TEC - either 'iri' for IRI2016 model (must be installed), or 'code' to autodownload from CODE
        alpha (float): for CODE only, estimate of ratio of TEC towards 'to the satellite only'. If 'auto', it will estimate it using iri. 0.85 is good value
    '''
    if glon > 180:
        glon = glon - 180  # fix round-earth
    if returnhei and source == 'code':
        print('WARNING, height is estimated only through IRI model, now setting to it')
        source = 'iri'
    if alpha == 'auto':
        getalpha=True
    else:
        getalpha=False
    altkmrange = [0, altitude, altitude]
    TECs = []
    heis = []
    source='code' ##change this 
    for acqtime in acq_times:
        if source == 'iri':
            iri_acq = iri2016.IRI(acqtime, altkmrange, glat, glon )
            # print('it is done, till here!')
            # print(iri_acq)
            TECs.append(iri_acq.TEC.values[0])
            heis.append(iri_acq.hmF2.values[0])
        elif source == 'code':
            if getalpha:
                iri_acq_gps = iri2016.IRI(acqtime, [0, 20000, 20000], glat, glon )
                iri_acq = iri2016.IRI(acqtime, altkmrange, glat, glon )
                alpha = float(iri_acq.TEC/iri_acq_gps.TEC)
            try:
                tec = get_vtec_from_code(acqtime, glat, glon)
            except:
                # CODE data is not available earlier than with 6 months delay.. or more?
                tec = np.nan
                print('No CODE data for date '+str(acqtime.date())+'. Setting NaN.')
            # decrease the value by some alpha... we expect alpha % of TEC being below the satellite.. should be improved
            #alpha = 0.85
            tec = alpha*tec
            TECs.append(tec)
    if returnhei:
        return TECs, heis
    else:
        return TECs



def download_code_data(acqtime, storedir = '/gws/nopw/j04/nceo_geohazards_vol1/code_iono'):
    instrings = ['CODG', 'CGIM']
    #lastrings = ['I', 'N']
    ffound = False
    for instr in instrings:
        #for lastr in lastrings:
        if not ffound:
            filename = instr + acqtime.strftime('%j') + '0.' + acqtime.strftime('%y') + 'I.Z'
            url = 'http://ftp.aiub.unibe.ch/CODE/' + acqtime.strftime('%Y') + '/' + filename
            fullpath = os.path.join(storedir, filename)
            ionix = fullpath[:-2]
            if not os.path.exists(ionix):
                if not os.path.exists(fullpath):
                    # download this
                    try:
                        wget.download(url, out=storedir)
                    except:
                        ffound = False
                if os.path.exists(fullpath):
                    ffound = True
            else:
                ffound = True
    # since 12/2022 they changed naming convention to e.g. COD0OPSFIN_20230510000_01D_01H_GIM.INX.gz
    # see https://cddis.nasa.gov/Data_and_Derived_Products/GNSS/atmospheric_products.html#iono
    if not ffound:
        filename = 'COD0OPSFIN_'+ acqtime.strftime('%Y') + acqtime.strftime('%j')+'0000_01D_01H_GIM.INX.gz' # TODO: check YMD
        #filename = instr + acqtime.strftime('%j') + '0.' + acqtime.strftime('%y') + 'I.Z'
        url = 'http://ftp.aiub.unibe.ch/CODE/' + acqtime.strftime('%Y') + '/' + filename
        fullpath = os.path.join(storedir, filename)
        ionix = fullpath[:-3]
        if not os.path.exists(ionix):
            if not os.path.exists(fullpath):
                # download this
                try:
                    wgotfile = wget.download(url, out=storedir)
                except:
                    print('error during wget download')
                    ffound = False
            if os.path.exists(fullpath):
                ffound = True
        else:
            ffound = True
    if not ffound:
        print('no CODE layer found for '+filename)
        return False
    if not os.path.exists(ionix):
        rc = os.system('cd ' + storedir + '; 7za x ' + filename + ' >/dev/null 2>/dev/null; rm ' + fullpath)
    if not os.path.exists(ionix):
        print('ERROR: maybe you do not have 7za installed')
        return False
    return ionix


def get_vtec_from_code(acqtime, lat = 0, lon = 0, storedir = '/gws/nopw/j04/nceo_geohazards_vol1/code_iono', return_fullxr = False):
    """ Adapted from Reza Bordbari script, plus using functions from https://notebook.community/daniestevez/jupyter_notebooks/IONEX

    Args:
        acqtime (dt.datetime)
        lat (float)
        lon (float)
        storedir (str)
        return_fullxr (bool): if True, will return full TEC datacube
    """
    #D = acqtime.strftime('%Y%m%d')
    #ipp = np.array([lat,lon])
    # check if exists:
    fna = glob.glob(storedir+'/????'+ acqtime.strftime('%j') + '0.' + acqtime.strftime('%y')+ '?')
    if fna:
        ionix = os.path.join(storedir,fna[0])
    else:
        # download:
        ionix = download_code_data(acqtime, storedir)
    if not ionix:
        return False
    #else:
    #    rc=os.system('rm '+fullpath) # clean the .Z
    # prep 
    #hhmmss=acqtime.strftime('%H%M%S')
    # loading the TEC maps, thanks to https://notebook.community/daniestevez/jupyter_notebooks/IONEX (but improved towards xarray by ML B-)
    tecmaps = get_tecmaps(ionix)
    interval = int(grep1line('INTERVAL',ionix).split()[0])
    timestep = interval/3600
    timecoords = np.arange(0.0,24.0+timestep,timestep)  # we expect start/end time being midnight, should be standard for all CODE files?
    lat_all = np.arange(87.5,-87.5-2.5,-2.5)
    lon_all = np.arange(-180.0,180.0+5,5.0)
    tecxr = xr.DataArray(data=tecmaps, dims=['time','lat','lon'],
                        coords=dict(time=timecoords, lon=lon_all, lat=lat_all) )
    # interpolate through the nan values
    tonan=9999
    tecxr.where(tecxr!=tonan)
    tecxr=tecxr.interpolate_na(dim="lon", method="linear", fill_value="extrapolate")
    tecxr = tecxr*1e+16 # from TECU
    if return_fullxr:
        return tecxr
    else:
        return get_vtec_from_tecxr(tecxr, acqtime, lat, lon)


# get_vtec_from_code(acqtime, lat, lon, storedir = '/gws/nopw/j04/nceo_geohazards_vol1/code_iono', return_fullxr = False):
def get_vtec_from_tecxr(tecxr, acqtime, lat, lon, rotate=True):
    '''Function to be used with tecxr (output from get_vtec_from_code) to get the tec for given coords'''
    h_time = float(acqtime.strftime('%H'))
    m_time = float(acqtime.strftime('%M'))
    s_time = float(acqtime.strftime('%S'))
    # given time in decimal format
    time_dec = h_time + (m_time/60) + (s_time / 3600)
    # ML: 2023/08, based on : https://github.com/insarlab/MintPy/blob/main/src/mintpy/objects/ionex.py
    # that is actually based on
    # Schaer, S., Gurtner, W., & Feltens, J. (1998). IONEX: The ionosphere map exchange format
    #         version 1.1. Paper presented at the Proceedings of the IGS AC workshop, Darmstadt, Germany.
    if rotate:
        # 3D interpolation with rotation as above reference
        htimes = tecxr.time.values
        pretime = int(htimes[htimes <= time_dec][-1])
        postime = int(htimes[htimes >= time_dec][0])
        #
        lon0 = lon + (time_dec - pretime) * 360. / 24.
        lon1 = lon + (time_dec - postime) * 360. / 24.
        #
        tec_val0 = float(tecxr.interp(time=pretime, lon=lon0, lat=lat, method='linear'))
        tec_val1 = float(tecxr.interp(time=postime, lon=lon1, lat=lat, method='linear'))
        #
        tec = ((postime - time_dec) / (postime - pretime) * tec_val0
                   + (time_dec - pretime) / (postime - pretime) * tec_val1)
    else:
        # previous attempt, but still too different from the S1_ETAD CODE outputs (that rotates the Earth towards the Sun..)
        tec = float(tecxr.interp(time=time_dec, lon=lon,lat=lat, method='cubic')) # should be better than linear, but maybe quadratic is more suitable?
    return tec


def parse_map(tecmap, exponent = -1):
    tecmap = re.split('.*END OF TEC MAP', tecmap)[0]
    return np.stack([np.fromstring(l, sep=' ') for l in re.split('.*LAT/LON1/LON2/DLON/H\\n',tecmap)[1:]])*10**exponent


def get_tecmaps(filename):
    exponent = int(grep1line('EXPONENT',filename).split()[0]) # this is exponent of the data
    with open(filename) as f:
        ionex = f.read()
        return [parse_map(t, exponent) for t in ionex.split('START OF TEC MAP')[1:]]


def get_tec(tecmap, lat, lon):
    i = round((87.5 - lat)*(tecmap.shape[0]-1)/(2*87.5))
    j = round((180 + lon)*(tecmap.shape[1]-1)/360)
    return tecmap[i,j]


def get_lat_block(IONEXFile = None,fullDate = None,time_dec = None,time_sort = None,lat_sort = None): 
    fullDate1 = fullDate
    if time_dec > 23:
        time_sort[2] = 0
        fullDate2 = fullDate + np.array([[0],[0],[1]])
    else:
        fullDate2 = fullDate
    
    # number of lines in the IONEX file.
    NoL = CalcNumOfLines(IONEXFile)
    # find the epochs and the two latitude blocks per epoch
    DataFile = open(IONEXFile,'r')
    fseek(DataFile,0,'bof')
    tline = fgetl(DataFile)
    for i in np.arange(1,NoL+1).reshape(-1):
        # epoch one
        e1 = strfind(tline,np.array([num2str(transpose(np.array([[fullDate1],[time_sort(1)],[0],[0]]))),'                        ','EPOCH OF CURRENT MAP']))
        if not (len(e1)==0) :
            condition = 1
            while condition == 1:
                # lat one
                l1 = strfind(tline,np.array([num2str(lat_sort(1),'%.1f'),'-180.0 180.0   5.0 450.0']))
                if not (len(l1)==0) :
                    a = fgetl(DataFile)
                    b = a
                    for j in np.arange(1,4+1).reshape(-1):
                        b = append(b,fgetl(DataFile))
                    latBlock1_t1 = str2num(b)
                tline = fgetl(DataFile)
                # lat two
                l2 = strfind(tline,np.array([num2str(lat_sort(2),'%.1f'),'-180.0 180.0   5.0 450.0']))
                if not (len(l2)==0) :
                    a = fgetl(DataFile)
                    b = a
                    for j in np.arange(1,4+1).reshape(-1):
                        b = append(b,fgetl(DataFile))
                    latBlock2_t1 = str2num(b)
                    condition = 2
        # epoch two
        e2 = strfind(tline,np.array([num2str(transpose(np.array([[fullDate2],[time_sort(2)],[0],[0]]))),'                        ','EPOCH OF CURRENT MAP']))
        if not (len(e2)==0) :
            condition = 1
            while condition == 1:
                # lat one
                l1 = strfind(tline,np.array([num2str(lat_sort(1),'%.1f'),'-180.0 180.0   5.0 450.0']))
                if not (len(l1)==0) :
                    a = fgetl(DataFile)
                    b = a
                    for j in np.arange(1,4+1).reshape(-1):
                        b = append(b,fgetl(DataFile))
                    latBlock1_t2 = str2num(b)
                clear('a','b','j')
                tline = fgetl(DataFile)
                # lat two
                l2 = strfind(tline,np.array([num2str(lat_sort(2),'%.1f'),'-180.0 180.0   5.0 450.0']))
                if not (len(l2)==0) :
                    a = fgetl(DataFile)
                    b = a
                    for j in np.arange(1,4+1).reshape(-1):
                        b = append(b,fgetl(DataFile))
                    latBlock2_t2 = str2num(b)
                    condition = 2
            break
        tline = fgetl(DataFile)
    return latBlock1_t1,latBlock2_t1,latBlock1_t2,latBlock2_t2


#get satellite position - ECEF
def aer2ecef(azimuthDeg, elevationDeg, slantRange, obs_lat, obs_long, obs_alt):
    '''
    obs_alt must be in metres
    solution by https://stackoverflow.com/questions/15954978/ecef-from-azimuth-elevation-range-and-observer-lat-lon-alt
    '''
    sitex, sitey, sitez = latlonhei2ecef(obs_lat,obs_long,obs_alt)
    #some needed calculations
    slat = np.sin(np.radians(obs_lat))
    slon = np.sin(np.radians(obs_long))
    clat = np.cos(np.radians(obs_lat))
    clon = np.cos(np.radians(obs_long))
    #
    azRad = np.radians(azimuthDeg)
    elRad = np.radians(elevationDeg)
    #
    # az,el,range to sez convertion
    south  = -slantRange * np.cos(elRad) * np.cos(azRad)
    east   =  slantRange * np.cos(elRad) * np.sin(azRad)
    zenith =  slantRange * np.sin(elRad)
    #
    #
    x = ( slat * clon * south) + (-slon * east) + (clat * clon * zenith) + sitex
    y = ( slat * slon * south) + ( clon * east) + (clat * slon * zenith) + sitey
    z = (-clat *        south) + ( slat * zenith) + sitez
    #
    return x, y, z


def latlonhei2ecef(lat, lon, alt):
    '''
    altitude should be in metres!!!!!
    '''
    transformer = pyproj.Transformer.from_crs(
        {"proj":'latlong', "ellps":'WGS84', "datum":'WGS84'},
        {"proj":'geocent', "ellps":'WGS84', "datum":'WGS84'},
        )
    x, y, z = transformer.transform(lon, lat, alt, radians=False)
    return x, y, z


def ecef2latlonhei(x, y, z):
    transformer = pyproj.Transformer.from_crs(
        {"proj":'geocent', "ellps":'WGS84', "datum":'WGS84'},
        {"proj":'latlong', "ellps":'WGS84', "datum":'WGS84'},
        )
    lon, lat, alt = transformer.transform(x,y,z,radians=False)
    return lat, lon, alt


def get_altitude(lat, lon):
    '''
    uses USGS site to get elevation data. thanks to:
    https://gis.stackexchange.com/questions/338392/getting-elevation-for-multiple-lat-long-coordinates-in-python
    '''
    #
    # USGS Elevation Point Query Service
    #url = r'https://nationalmap.gov/epqs/pqs.php?'
    #
    # opentopodata
    #url = r'https://api.opentopodata.org/v1/eudem25m?locations=51.875127,-3.341298
    url = r'https://api.opentopodata.org/v1/etopo1?locations={0},{1}'.format(lat, lon)
    #
    result = requests.get(url) # + urllib.parse.urlencode(params)))
    elev = result.json()['results'][0]['elevation']
    #print(elev)
    if float(elev) < -100:
        elev = 0
    return elev



def calculate_daz_iono(bovl, esds, framespd, method = 'gradient', out_hionos = False, out_tec_master = False, out_tec_all = False, ionosource='iri', use_iri_hei=False):
    ''' Function to calculate iono correction for a given frame.

    Args:
        frame (str)                 frame ID/BOVL id for esd!
        esds (pandas.Dataframe)     standard esds table
        framespd (pandas.Dataframe) standard framespd table
        method (str)                gradient or liang (gradient is correct, Liang is kept from first attempts, may not work)
        out_hionos (bool)           whether to output also IRI-estimated peak height of the F2 iono layer (sadly, CODE calculates for 450 km hei which is probably wrong)
        out_tec_master (bool)
        out_tec_all (bool)
        ionosource (str)            iri or code (for IRI2016 or CODE GIM-based ionosphere. the latter improves RMSE!)
        use_iri_hei (bool)          if True, it estimates F2 peak altitude using IRI2016 and uses for the correction (with any ionosource). NOTE, CODE data are for 450 km ALT! so better not use_iri_hei..

    Notes: 'liang' method should include also some extra F2 height correction..
    2023/08: Liang method was first implemented here and a lot happened since that time.. Please consider it obsolete.
    '''
    if method == 'gomba': # renamed it to keep as it is
        method = 'gradient'
    if ionosource == 'iri' and (not use_iri_hei):
        print('using IRI2016, setting the iri to estimate F2 peak altitude')
        use_iri_hei=False  ##I change
    selected_frame_esds = esds[esds['swath_overlap'] == bovl].copy()
    frameta = framespd[framespd['swath_overlap']==bovl]
    # extract some variables
    heading = frameta['heading'].values[0]
    scene_center_lon = frameta['center_lon'].values[0]
    scene_center_lat = frameta['center_lat'].values[0]
    resolution_of_pixel = frameta['az_pix'].values[0]

    try:
        range_avg = frameta['center_range'].values[0] # 2024: GAMMA had a bug wrongly informing on centre_range (may differ by 20 km or so!). fixed most of it
    except:
        range_avg = frameta['center_range'].values[0]
    master = frameta['master'].values[0]
    inc_angle_avg = frameta['incidence'].values[0]
    center_time=frameta['center_time'].values[0]
    dfDC = frameta['dfDC'].values[0]
    ##okey so far
    if dfDC == 0:
        print('warning, bovl '+bovl+' has no info on dfDC, using default')
        dfDC = 4365 #Hz, mean value in the whole dataset
    #a bit recalculate
    theta = np.radians(inc_angle_avg)
    #sathei = int(range_avg * np.cos(theta)/1000) #in km --- will do this better
    master_time = pd.to_datetime(str(master)+'T'+center_time)
    acq_times = pd.to_datetime(selected_frame_esds.epoch.astype(str)+'T'+center_time)
    # include master time!
    acq_times[acq_times.index[-1]+1] = master_time
    #
    # 1. get f2 hei inbetween target center point C and nadir of the satellite satg
    wgs84 = nv.FrameE(name='WGS84')
    Pscene_center = wgs84.GeoPoint(latitude=scene_center_lat, longitude=scene_center_lon, degrees=True)
    burst_len = 7100*2.758277 #approx. satellite velocity on the ground 7100 [m/s] * burst_interval [s]
    ###### do the satg_lat, lon
    azimuthDeg = heading-90 #yes, azimuth is w.r.t. N (positive to E)
    elevationDeg = 90-inc_angle_avg
    slantRange = range_avg
    #try:
    #    scene_alt = get_altitude(scene_center_lat, scene_center_lon)
    #except:
    #    scene_alt = 0
    #to get position of the satellite - UNCLEAR about slantRange - is this w.r.t. DEM? (should ask GAMMA) - if not (only in WGS-84), i should use scene_alt=0!
    # 2023/08: checked using orbits - the slantranfe is wrt ellipsoid! setting scene alt 0
    x, y, z = aer2ecef(azimuthDeg, elevationDeg, slantRange, scene_center_lat, scene_center_lon, 0) #scene_alt)
    satg_lat, satg_lon, sat_alt = ecef2latlonhei(x, y, z)
    Psatg = wgs84.GeoPoint(latitude=satg_lat, longitude=satg_lon, degrees=True)
    # get middle point between scene and sat - and get F2 height for it
    path = nv.GeoPath(Pscene_center.to_nvector(), Psatg.to_nvector())
    # get point in the middle
    Pmid_scene_sat = path.interpolate(0.5).to_geo_point()
    # work in dedicated table
    df = pd.DataFrame(acq_times)
    # 2024: not clear if IRI2016 F2 peak altitude is correct. Allowing the standard 450 km assumption by CODE (I think TECs will get scaled, so the gradient should be still ok. Not tested)
    if use_iri_hei or out_hionos:
        # get hionos in that middle point:
        tecs, hionos = get_tecs(Pmid_scene_sat.latitude_deg, Pmid_scene_sat.longitude_deg, 800, acq_times, returnhei = True)
        hiono_master = hionos[-1]
        selected_frame_esds['hiono'] = hionos[:-1]  ###*1000 # convert to metres, avoid last measure, as this is 'master'
        df['hiono'] = hionos
    else:
        df['hiono'] = 450 # standard altitude used by CODE
        hionos = list(df.hiono.values)
        hiono_master = 450
        selected_frame_esds['hiono'] = 450
        print('No select, hiono as 450')
        #
    ############## now calculate TEC using the SLM knowledge, i.e. different A,B per epoch (!)
    # (note that the last hiono is for the master/reference epoch
    print('iri calculated correctly!')
    tecs_A = []
    tecs_B = []
    #for hiono in hionos:
    for i,a in df.iterrows():
        hiono = a['hiono']*1000 # m
        epochdate = a['epoch']
        # first, get IPP - ionosphere pierce point
        # range to IPP can be calculated using:
        # range_IPP = hiono/np.sin(theta)
        # but maybe not... let's do it simpler:
        range_IPP = slantRange * hiono / sat_alt
        x, y, z = aer2ecef(azimuthDeg, elevationDeg, range_IPP, scene_center_lat, scene_center_lon, 0) #scene_alt)
        ippg_lat, ippg_lon, ipp_alt = ecef2latlonhei(x, y, z)
        Pippg = wgs84.GeoPoint(latitude=ippg_lat, longitude=ippg_lon, degrees=True)
        # then get A', B'
        PsatgA, _azimuth = Psatg.displace(distance=burst_len/2, azimuth=heading-180, method='ellipsoid', degrees=True)
        PsatgB, _azimuth = Psatg.displace(distance=burst_len/2, azimuth= heading, method='ellipsoid', degrees=True)
        # then do intersection ...
        PippAt, _azimuth = Pippg.displace(distance=burst_len, azimuth=heading-180, method='ellipsoid', degrees=True)
        PippBt, _azimuth = Pippg.displace(distance=burst_len, azimuth= heading, method='ellipsoid', degrees=True)
        path_ipp = nv.GeoPath(PippAt, PippBt)
        path_scene_satgA = nv.GeoPath(Pscene_center, PsatgA)
        path_scene_satgB = nv.GeoPath(Pscene_center, PsatgB)
        # these two points are the ones where we should get TEC
        PippA = path_ipp.intersect(path_scene_satgA).to_geo_point()
        PippB = path_ipp.intersect(path_scene_satgB).to_geo_point()
        ######### get TECS for A, B
        TECV_A = get_tecs(PippA.latitude_deg, PippA.longitude_deg, round(sat_alt/1000), [epochdate], False, source=ionosource)[0]
        TECV_B = get_tecs(PippB.latitude_deg, PippB.longitude_deg, round(sat_alt/1000), [epochdate], False, source=ionosource)[0]
        # get inc angle at IPP - see iono. single layer model function
        earth_radius = 6378160 # m
        sin_thetaiono = earth_radius/(earth_radius+hiono) * np.sin(theta)
        TECS_A = TECV_A/np.sqrt(1-sin_thetaiono**2)
        TECS_B = TECV_B/np.sqrt(1-sin_thetaiono**2)
        tecs_A.append(TECS_A)
        tecs_B.append(TECS_B)
    #
    tec_A_master = tecs_A[-1]
    tec_B_master = tecs_B[-1]
    tecs_A = tecs_A[:-1]
    tecs_B = tecs_B[:-1]
    #
    
    selected_frame_esds['TECS_A'] = tecs_A
    selected_frame_esds['TECS_B'] = tecs_B
    
    #
    ##############################
    #if method == 'gomba':
    PRF = 486.486
    k = 40.308193 # m^3 / s^2
    f0 = 5.4050005e9
    c = speed_of_light
    fH = f0 + dfDC*0.5
    fL = f0 - dfDC*0.5
    #tecovl = (TECs_B1 - TEC_master_B1)/(fH*fH) - (TECs_B2 - TEC_master_B2)/(fL*fL)
    #daz_iono = -2*PRF*k*f0/c/dfDC * tecovl
    if method == 'gradient':
        # 08/2021 - empirically checked, correct:
        #tecovl = (selected_frame_esds['TECS_B'] - tec_B_master)/(fL*fL) - (selected_frame_esds['TECS_A'] - tec_A_master)/(fH*fH)
        #daz_iono = 2*PRF*k*f0/c/dfDC * tecovl
        # 04/2022 - actually the squares seem not needed, based directly on iono2phase (see article):
        # note it is 'master' minus 'epoch S', since i start from ifg as M * complex conjugate of S (so basically phase M - phase S, phase is negative for higher delay)
        # 2023 - oh well, it was best to do straightforward equations, see Lazecky et al., 2023, GRL
        tecovl = (tec_A_master - selected_frame_esds['TECS_A'])/fH - (tec_B_master - selected_frame_esds['TECS_B'])/fL
        daz_iono = 2*PRF*k/c/dfDC * tecovl
        
    else:
        # following the Liang 2019:
        # needs ka...
        hsat =  sat_alt
        ka = frameta['ka'].values[0]
        #print('setting Hiono={} km, but it would be close to Gomba method as around 700 km'.format(hiono_m))
        # kion should be rad/s
        # this would get iono effect within the full burst spectrum (from fH to fL)
        phaionL = -4*np.pi*k/c/fL * selected_frame_esds['TECS_B']  # - TEC_master_B1)
        phaionH = -4*np.pi*k/c/fH * selected_frame_esds['TECS_A']  # - TEC_master_B2)
        phaionL_m = -4*np.pi*k/c/fL * tec_B_master
        phaionH_m = -4*np.pi*k/c/fH * tec_A_master
        #
        #phi_ionoramp_burst = fL*fH/(f0*(fH*fH - fL*fL)) * (dphaionL*fH - dphaionH*fL)
        phi_ionoramp_burst = phaionH - phaionL
        phi_ionoramp_burst_m = phaionH_m - phaionL_m
        #i suppose i relate it to burst overlap interval (i.e. 2.7 s /10)
        #burst_ovl_interval = 2.758277/10 #s
        burst_interval = 2.758277
        #kion = phi_ionoramp_burst / burst_ovl_interval # rad/s
        kion = phi_ionoramp_burst / burst_interval # rad/s
        kion_m = phi_ionoramp_burst_m / burst_interval # rad/s
        nesd_iono = -(1/(2*np.pi)) * kion/ka * ((selected_frame_esds['hiono'])/hsat + 0.5) # only for one acq!!!
        nesd_iono_master = -(1/(2*np.pi)) * kion_m/ka * ((hiono_master)/hsat + 0.5) # only for one acq!!!
        # in seconds, need to convert to pixels
        #phiesd = nesd_iono * 2 * np.pi * dfDC
        #daz_iono = PRF * phiesd / (2*np.pi*dfDC)
        daz_iono = PRF * (nesd_iono - nesd_iono_master) # / (2*np.pi*dfDC)
        print('it is okay till here')
        #return daz_iono
    if out_hionos:
        if out_tec_master:
            return daz_iono, hionos, tec_A_master, tec_B_master
        elif out_tec_all:
            # experiment...
            return daz_iono, hionos, tec_A_master, tec_B_master, tecs_A, tecs_B
        else:
            return daz_iono, hionos
    else:
        if out_tec_all:
            return daz_iono, tec_A_master, tec_B_master, tecs_A, tecs_B
        else:
            return daz_iono



# extra step: get hours from dawn, or dusk respectively:
def get_hours_from_dusk_dawn(framespd):
    sun = ephem.Sun()
    location = ephem.Observer()
    framespd['hours_from_dusk_dawn'] = 0.0
    #
    for i, frameta in framespd.iterrows():
        if frameta['centre_time'] == '0':
            continue
        lon = frameta['center_lon']
        lat = frameta['center_lat']
        master_time = pd.to_datetime(str(frameta['master'])+'T'+frameta['centre_time'])
        location.lat = str(lat)
        location.lon = str(lon)
        location.date = master_time
        try:
            a = location.previous_rising(sun).datetime()
            aa = (master_time - a).total_seconds()/60/60
            b = location.next_rising(sun).datetime()
            bb = (master_time - b).total_seconds()/60/60
            c = location.previous_setting(sun).datetime()
            cc = (master_time - c).total_seconds()/60/60
            d = location.next_setting(sun).datetime()
            dd = (master_time - d).total_seconds()/60/60
        except:
            continue
        hours = 9999
        for x in [aa,bb,cc,dd]:
            if abs(x) < abs(hours):
                hours = x
        #if sunminutes <0:
        #    sunrise = location.previous_transit(sun).datetime()
        #bb = master_time.time()
        #timediff = (bb.hour*60 + bb.minute) - (sunrise.hour*60 + sunrise.minute)
        framespd.at[i, 'hours_from_dusk_dawn'] = hours
    return framespd