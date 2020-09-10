#!/usr/bin/env python
# coding: utf-8

## Wind Erosion Model
## Written by Charlie Weil & Isita Talukdar, August 2020
## Inspired by Huang Binbin's approach.
import numpy as np
import numpy.ma as ma
import gdal
import os
import matplotlib.pyplot as plt #visualisation

def make_directory(path):
    if not os.path.exists(path):
        os.makedirs(path)
        
# - - - - - - - - - - -
#     File paths 
# - - - - - - - - - - -

data_dir ="data_UTM/"
input_data_path = os.path.join(data_dir, 'input/')
intermediate_data_path = os.path.join(data_dir, 'intermediate/')
make_directory(intermediate_data_path)
output_data_path = os.path.join(data_dir, 'output/')
make_directory(output_data_path)

dem_filename = 'dem_aoi.tif'
dem_file_path = os.path.join(input_data_path, dem_filename)
    
# Climate Input Data

def temperature_file_path(month_num):
    temperature_filename = 'temperature_degC/' + 'wc2.0_30s_tave_'+ month_num + '.tif'
    return os.path.join(input_data_path, temperature_filename)

def precipitation_file_path(month_num):
    precipitation_filename = 'precipitation_mm/' + 'wc2.0_30s_prec_' + month_num + '.tif'
    return os.path.join(input_data_path, precipitation_filename)

def solar_rad_file_path(month_num):
    solar_rad_filename = 'solar_radiation/'+ 'shortwave_radiation_' + month_num + '.tif'
    return os.path.join(input_data_path, solar_rad_filename)

def rain_days_file_path(month_num):
    prcp_days_filename = 'month_prcp_day/prcp_day_' + month_num + '.tif'
    return os.path.join(input_data_path, prcp_days_filename)

def snow_cover_file_path(month_num):
    snow_cover_filename = 'snow_gm2/' + 'snow_' + month_num + '.tif'
    return os.path.join(input_data_path, snow_cover_filename)

def wind_spd_file_path(month_num):
    wind_speed_filename = 'wind_speed_monthly_clipped/' + 'wind_speed_' + month_num + '.tif'
    return os.path.join(input_data_path, wind_speed_filename)

def weather_factor_out_file_path(month_num):
    wf_file_path = intermediate_data_path + 'WF/'
    make_directory(wf_file_path)
    wf_out_filename = 'wf_' + month_num + '.tif'
    return os.path.join(wf_file_path, wf_out_filename)

# Soil Input Data
sand_filename = 'soil/sand.tif'
sand_file_path = os.path.join(input_data_path, sand_filename)

silt_filename = 'soil/silt.tif'
silt_file_path = os.path.join(input_data_path, silt_filename)

clay_filename = 'soil/clay.tif'
clay_file_path = os.path.join(input_data_path, clay_filename)

org_mat_filename = 'soil/soil_organic_matter_gm2.tif'
org_mat_file_path = os.path.join(input_data_path, org_mat_filename)

scf_filename = 'SCF.tif'
scf_file_path = os.path.join(intermediate_data_path, scf_filename)

ef_filename = 'EF.tif'
ef_file_path = os.path.join(intermediate_data_path, ef_filename)

#Vegetation Factor and K' Input Data
def frac_veg_cov_file_path(month_num):
    frac_veg_cov_filename = 'vegetation_percent_cover/' + 'vegetation_percent_cover_' + month_num + '.tif'
    return os.path.join(input_data_path, frac_veg_cov_filename)

def cog_out_file_path(month_num):
    cog_file_path = intermediate_data_path + "COG/"
    make_directory(cog_file_path)
    cog_out_filename = 'COG_' + month_num + '.tif'
    return os.path.join(cog_file_path, cog_out_filename)

def kk_out_file_path(month_num):
    kk_file_path = intermediate_data_path + "KK/"
    make_directory(kk_file_path)
    kk_out_filename = 'KK_' + month_num + '.tif'
    return os.path.join(kk_file_path, kk_out_filename)

#Actual and Predicted Wind Erosion Output Paths

def sl_actual_out_file_path():
    sl_actual_out_filename = 'SL_actual.tif'
    return os.path.join(output_data_path, sl_actual_out_filename)

def sl_wo_veg_out_file_path():
    sl_wo_veg_out_filename = 'SL_without_veg.tif'
    return os.path.join(output_data_path, sl_wo_veg_out_filename)

def sl_fvc_100_out_file_path():
    sl_fvc_100_out_filename = 'SL_fvc_100.tif'
    return os.path.join(output_data_path, sl_fvc_100_out_filename)

def sand_r_out_file_path():
    sand_re_out_filename = 'sand_re.tif'
    return os.path.join(output_data_path, sand_re_out_filename)

# - - - - - - - - - - -
#      Parameters
# - - - - - - - - - - -
HIST_BINS = 50
Kr_WINDOW_SIZE = 3
en_corr_plot = True
en_hist_plot = False

# days of month 1-12
mondays = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
endday = [30, 58, 89, 119, 150, 180, 211, 242, 272, 303, 333, 364]
#month_id = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'] ## Isita, this variable is not used, is it normal?



# - - - - - - - - - - -
#    Model functions 
# - - - - - - - - - - -

def execute(args):
    ## This is where the steps in "Execute" shoulf eventually go.
    return None

def RasterSave(data, path, ref_file_path):
    """Save Data to Raster File
    Parameters:
        data: array
        path: path where new Raster should be saved
        ref_file_path: file path of reference file used to set dimensions of new raster
        
    Returns:
        --
    """
    ref_raster = gdal.Open(ref_file_path)
    ndv = -9999
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(path, ref_raster.RasterXSize, ref_raster.RasterYSize, 1, gdal.GDT_Float32)
    ds.SetGeoTransform(ref_raster.GetGeoTransform())
    ds.SetProjection(ref_raster.GetProjection())
    data[np.isnan(data)] = ndv
    ds.GetRasterBand(1).WriteArray(data)
    ds.GetRasterBand(1).SetNoDataValue(ndv)
    data[data == ndv] = np.nan
    del ds
    
def read_raster_as_array(file_path):
    """Read Raster Data to an array
    Parameters:
        file_path: path where source Raster is found
        
    Returns:
        GDAL array: raster data as an array, no-data values as nan
    """
    raster = gdal.Open(file_path)
    band = raster.GetRasterBand(1)
    ndv = band.GetNoDataValue()
    array = band.ReadAsArray(0, 0, raster.RasterXSize, raster.RasterYSize)
    array = array.astype('float')
    DataType = band.DataType
    DataType = gdal.GetDataTypeName(DataType)
    # Convert No Data Points to nans
    array[array == ndv] = np.nan
    return array

def calculate_air_pressure(dem_file_path):
    """Calculate air pressure
    Parameters:
        dem_file_path: path to Elevation data, in meters
        
    Returns:
        Air Pressure (P), in kPa
        GDAL array
    """
    dem_a2 = read_raster_as_array(dem_file_path)
    air_pressure = 101.3 * (1 - 0.0255 * dem_a2 / 1000 * (6357 / (6357 + dem_a2 / 1000))) ** 5.256
    return air_pressure

def get_monthly_avg_temp(file_path):
    """Get Average Temperature for a given month from Raster File
    Parameters:
        file_path: path to Raster file with average temperature
    Returns:
        tem_m_a2: Temperature for a given month, in °C
        GDAL array 
    """
    
    # monthly average temperature(℃)
    tem_m = read_raster_as_array(file_path)
    tem_m = (tem_m < -100) * 0.1 + (tem_m >= -100) * tem_m
    return tem_m

def get_monthly_total_precip(file_path):
    """Get Total Precipitation for a given month from Raster File
    Parameters:
        file_path: path to Raster file with total precipitation
    Returns:
        prcp_a: Total Precipitation for a given month, in mm
        GDAL array 
    """
    # monthly total precipitation(mm)
    prcp = read_raster_as_array(file_path)
    prcp = (prcp < 0) * 0 + (prcp >= 0) * prcp
    return prcp

def get_monthly_sol_rad(file_path):
    """Get Total Solar Radiation for a given month from Raster File
    Parameters:
        file_path: path to Raster file with total solar radiation
    Returns:
        SOL_a2: Total Solar Radiation for a given month, in MJ/m2
        GDAL array 
    """

    # monthly total solar radiation (MJ/m2)
    SOL = read_raster_as_array(file_path)
    SOL = (SOL < 0) * 0.1 + (SOL >= 0) * SOL
    return SOL

def get_monthly_num_rain_days(file_path):
    """Get Number of Rain Days for a given month from Raster File
    Parameters:
        file_path: path to Raster file with number of rain days
    Returns:
        prcp_days_a2: Total Number of Rain Days for a given month, in days
        GDAL array 
    """

    # days of rain events in every month
    prcp_days = read_raster_as_array(file_path)
    prcp_days = (prcp_days < 0) * 0 + (prcp_days >= 0) * prcp_days
    return prcp_days

def calculate_evapotranspiration(solar_rad, avg_temp):
    """Calculate the potential evapotranspiration for a given month
    Parameters:
        Temperature for a given month, in °C
        GDAL array
        
        Total Solar Radiation for a given month, in Megajoules/meters squared
        GDAL array
        
    Returns:
        Potential Evapotranspiration (ETp), in millimeters
        GDAL array
    """
    evap_trans = 0.0135 * (solar_rad / 2.54) * (avg_temp + 17.8)
    
    return evap_trans

def calculate_soil_moisture(solar_rad, avg_temp, precip, days_precip):
    """Calculate the Soil Moisture Factor for a given month
    Parameters:
        Potential Evapotranspiration for a given month, in millimeters
        GDAL array
        
        Number of Rain days for a given month
        GDAL array
        
    Returns:
        Soil Moisture Factor (SW)
        GDAL array
    """
    evap_trans = calculate_evapotranspiration(solar_rad, avg_temp)
    soil_moisture = (evap_trans - precip * days_precip) / evap_trans
    soil_moisture = (soil_moisture < 0) * 0 + (soil_moisture >= 0) * soil_moisture
    
    return soil_moisture

def calculate_snow_factor(file_path):
    """Calculate the Snow Factor for a given month
    Parameters:
        file_path: path to Raster file with Snow cover, 
    Returns:
        Snow Factor (SD),
        GDAL array
    """

    SD = read_raster_as_array(file_path)
    SD = (SD < 0) * 0 + (SD >= 0) * SD
    SD = (1 - SD * 0.01)
    return SD

def calculate_air_density(temperature, air_pressure):
    """Calculate air density for a given month.
    Parameters:
        temperature: Temperature for a given month, in °C
        GDAL array
        
        air_pressure: Air Pressure, in kPa
        GDAL array
    Returns:
        Air density (rho), in kg/m^3
        GDAL array
    """
    air_density_rho = 1.293 * (273 / (273 + temperature)) * air_pressure / 101.3
    return air_density_rho

def calculate_wind_factor_monthly(file_path):
    """Calculate the wind factor for a given month
    Parameters:
        file_path: path to Raster file with monthly wind speed
        
    Returns:
        Wind Factor (wf)
        GDAL array
    """
     # daily wind speed (m/s)
    wind_speed = read_raster_as_array(file_path)
    # wind factor
    wf = wind_speed * (wind_speed - 5) ** 2
    
    return wf

def calculate_monthly_weather_factor(wind_file_path, temp, precip, sol_rad, precip_days, snow_factor, pressure):
    """Calculate the Weather Factor for a given month
    Parameters:
        Wind Factor (wf)
        GDAL array
        
        Air density (rho), in kg/m^3
        GDAL array
        
        Soil Moisture Factor (SW)
        GDAL array
        
        Snow Factor (SD)
        GDAL array
        
    Returns:
        Weather Factor for a given month (WF)
        GDAL array
    """
   
    wind_factor =  calculate_wind_factor_monthly(wind_file_path)

    air_density_rho = calculate_air_density(temp, pressure)
    
    soil_moisture = calculate_soil_moisture(sol_rad, temp, precip, precip_days)
    
    weather_factor = wind_factor * air_density_rho / 9.8 * soil_moisture * snow_factor
    
    return weather_factor

def preprocess_soil_nonneg(file_path):
    """Get Soil Data(ratio of: sand, silt, clay, org matter) from Raster File and preprocess for non-negative values
    Parameters:
        file_path: path to Raster file with soil ratio
        
    Returns:
        soil_ratio : soil data preprocessesd
        GDAL array
    """

    #ratio of soil data(%)
    soil_ratio = read_raster_as_array(file_path)
    soil_ratio = (soil_ratio <= 0) * 0.1 + (soil_ratio >0) * soil_ratio
    return soil_ratio



def calculate_soil_crust_factor(clay_ratio, org_mat_ratio):
    """Calculate the Soil Crusting Factor
    Parameters:
        Clay Ratio(%)
        GDAL array
        
        Organic Matter Ratio(%)
        GDAL array
        
    Returns:
        Soil Crusting Factor (SCF)
        GDAL array
    """
    soil_crust_factor = 1 / (1 + 0.0066 * clay_ratio ** 2 + 0.021 * org_mat_ratio ** 2)

    return soil_crust_factor

def calculate_soil_erodibility_factor(sand_ratio, silt_ratio, clay_ratio, org_mat_ratio):
    """Calculate the Soil Erodibility Factor
    Parameters:
        Sand Ratio(%)
        GDAL array
        
        Silt Ratio(%)
        GDAL array
        
        Clay Ratio(%)
        GDAL array
        
        Organic Matter Ratio(%)
        GDAL array
        
    Returns:
        Soil Erodibility Factor (EF)
        GDAL array
    """
    soil_erode_factor = (29.09 + 0.31 * sand_ratio + 0.17 * silt_ratio + 0.33 * sand_ratio / clay_ratio - 2.59 * org_mat_ratio - 0.95 * 0) / 100


    return soil_erode_factor
def preprocess_veg_coverage(file_path):
    """Preprocess the Fraction of Monthly Vegetation Coverage For a Given Month
    Parameters:
        Fraction of Monthly Vegetation Coverage
        file_path: path to Raster file with Fraction of Monthly Vegetation Coverage
        
    Returns:
        fvc: fractional vegetation coverage 
        GDAL array
    """
    fvc = read_raster_as_array(file_path)
    fvc = ((fvc < 0) | (fvc > 100)) * 0 + ((fvc >= 0) & (fvc <= 100)) * fvc
    return fvc

def calculate_vegetation_factor(file_path):
    """Calculate the Fraction of Monthly Vegetation Coverage For a Given Month
    Parameters:
        Fraction of Monthly Vegetation Coverage
        file_path: path to Raster file with Fraction of Monthly Vegetation Coverage
        
    Returns:
        Vegetation Factor(COG)
        GDAL array
    """
    # read fraction of vegetation coverage data (%)
    fvc = preprocess_veg_coverage(file_path)
    vegetation_factor = np.exp(-0.00483*fvc)
    
    return vegetation_factor

def window_delta(data_arr, window_size):
    """Calculate Delta H(max elevation - min elevation)

    Parameters:
        data_arr: Elevation data, in meters
        read from array
        
        window size: number, width of section of pixels
        
    Returns:
        Roughness Length, in cm
        stored in array 
    """
    rows, columns = data_arr.shape
    temp_sum = np.zeros((window_size, window_size))
    deltas = np.zeros((rows, columns))
    
    # create a padded copy window only odd sized
    pad = int(window_size/2)
    matrix = np.pad(data_arr, pad, 'edge')
    window_rows = window_size
    window_cols = window_size
    # Level 2: traversing the window 
    for y in range(rows):
        for x in range(columns):
            # Level 1: handling the matrix 
            # (rows, columns = data.shape !)
            temp_sum = matrix[y : window_rows + y,
                              x : window_cols + x]
            #print(temp_sum.max(), temp_sum.min())
            deltas[y, x] = temp_sum.max() -  temp_sum.min()
              
    return deltas

def calculate_roughness_length(dem_file_path):
    """Calculate Roughness Length 

    Parameters:
        dem: Elevation data, in meters
        file path to elevation data
    Returns:
        Roughness Length, in cm
        stored in array 
    """
    dem_arr = read_raster_as_array(dem_file_path)
    raster = gdal.Open(dem_file_path)
    geotransform = raster.GetGeoTransform()
    
    deltas = window_delta(dem_arr, Kr_WINDOW_SIZE)
    #print(deltas)
    #print_max_min(deltas, 'delats')
    #convert to kilometers
    deltas = deltas/1000
    resolution = (abs(geotransform[1]) + abs(geotransform[5]))/2
    #print('res:',resolution)
    #convert to kilometers
    resolution = resolution/1000
    roughness_length = 0.2* ((deltas **2)/(resolution*(Kr_WINDOW_SIZE-1)))
    #conversion to centimeters
    #roughness_length *= 100
    return roughness_length
    
def calculate_chain_rand_roughness(fvc_path):
    """Calculate the Chain Random Roughness

    Parameters:
        Fraction of Monthly Vegetation Coverage
        fvc_path: name of Raster file with clay ratio
        Raster File: 1 pixel = 5532m
        
        
    Returns:
        Chain Random Roughness(Crr)
        stored in array
    """
    fvc_a = read_raster_as_array(fvc_path)
    fvc_a2 = ((fvc_a < 0) | (fvc_a > 100)) * 0 + ((fvc_a >= 0) & (fvc_a <= 100)) * fvc_a
    fvc_a2=fvc_a2/100
    Crr = 17.46 * (0.025 + 2.464*fvc_a2**3.56)**0.738
    return Crr
    
    
def calculate_surface_terr_rough(dem_file_path, fvc_path):
    """Calculate the Surface Terrain Roughness

    Parameters:
        Fraction of Monthly Vegetation Coverage
        fvc_path: name of Raster file with clay ratio
        Raster File: 1 pixel = 5532m
        
        dem: Elevation data, in meters
        
    Returns:
        Surface Terrain Roughness(K')
        stored in array
    """
    Kr  = calculate_roughness_length(dem_file_path)
    Crr = calculate_chain_rand_roughness(fvc_path)
    exp = 1.86*Kr - 2.41*Kr**0.934 - 0.127*Crr
    KK = np.exp(1.86*Kr - 2.41*(Kr**0.934) - 0.127*Crr)
    KK_a2 = (KK <= 0) * 0.1 + (KK > 0) * KK
    return KK_a2

def read_WF(file_path):
    """Read WF from Raster File to array
    Parameters:
        file_path: path to Raster file with Weather Factor
    Returns:
        weather_factor: Weather Factor data in an array
        GDAL array 
    """
    #reads the output TIF into an array
    weather_factor = read_raster_as_array(file_path)
    weather_factor = (weather_factor <= 0) * 0.1 + (weather_factor > 0) * weather_factor
    return weather_factor

def read_COG(file_path):
    """Read COG from Raster File to array
    Parameters:
        file_path: path to Raster file with Vegetation Factor
    Returns:
        vegetation_factor: Vegetation Factor data in an array
        GDAL array 
    """
    #reads the output TIF into an array
    vegetation_factor = read_raster_as_array(file_path)
    vegetation_factor = (vegetation_factor <= 0) * 0.1 + (vegetation_factor > 0) * vegetation_factor
    return vegetation_factor

def read_Kprime(file_path):
    """Read K' from Raster File to array
    Parameters:
        file_path: path to Raster file with Surface Terrain Roughness Factor
    Returns:
        vegetation_factor: Surface Terrain Roughness Factor data in an array
        GDAL array 
    """
    #reads the output TIF into an array
    KK = read_raster_as_array(file_path)
    return KK


def calculate_monthly_wind_erosion(weather_factor,soil_erode_factor, kprime, soil_crust_factor, cog, fvc_100):
    """Calculate Potential Wind Erosion for a given month
    Parameters:
        Potential RWEQ Maximum Horizontal Flux (Qmax)
        array from function
        
        Potential Critical Field Length (s)
        array from function
        
    Returns:
         Potential Wind Erosion for a given month (SL)
         stored in array
    """
    #Actual Wind Erosion 
    Qmax = 109.8 * (weather_factor * soil_erode_factor * kprime * soil_crust_factor * cog)
    s = 105.71 * (weather_factor * soil_erode_factor * kprime * soil_crust_factor * cog)**-0.3711
    wind_erosion = 100/(s * s +0.01)*Qmax * (np.exp(-(50 / (s + 0.01)) ** 2))

    #Potential Wind Erosion
    Qmaxp = 109.8 * (weather_factor * soil_erode_factor * kprime * soil_crust_factor)
    sp = 105.71 * (weather_factor * soil_erode_factor * kprime * soil_crust_factor)**-0.3711
    wind_erosionp = 100/(sp * sp +0.01)*Qmaxp * (np.exp(-(50 / (sp + 0.01)) ** 2))

    #100 Vegetation Coverage Wind Erosion 
    Qmax = 109.8 * (weather_factor * soil_erode_factor * kprime * soil_crust_factor * fvc_100)
    s = 105.71 * (weather_factor * soil_erode_factor * kprime * soil_crust_factor * fvc_100)**-0.3711
    wind_erosion_cog100 = 100/(s * s +0.01)*Qmax * (np.exp(-(50 / (s + 0.01)) ** 2))
    return wind_erosion, wind_erosionp,wind_erosion_cog100

def compute_corr(arr_x, arr_y):
    a=ma.masked_invalid(arr_x)
    b=ma.masked_invalid(arr_y)
    msk = (~a.mask & ~b.mask)
    return round(ma.corrcoef(a[msk],b[msk])[0,1], 3)

# - - - - - - - - - - -
#        Execute
# - - - - - - - - - - -
# # # # # # # # # # # # # # # # #

# Step 1 : Weather Factor
# # # # # # # # # # # # #

#calculate air pressure
pressure = calculate_air_pressure(dem_file_path)

#Compute Monthly Weather Factor 
for month_id in range(0, 12):
    mont_id_str = str(month_id+1)
    if month_id+1 < 10:
        mont_id_str = '0'+str(month_id+1)
    
    temp_file_path = temperature_file_path(mont_id_str)
    temp = get_monthly_avg_temp(temp_file_path)
    
    prcp_file_path = precipitation_file_path(mont_id_str)  
    precip = get_monthly_total_precip(prcp_file_path)
    
    sol_file_path = solar_rad_file_path(str(month_id + 1))
    sol_rad = get_monthly_sol_rad(sol_file_path)
    
    prcp_days_file_path = rain_days_file_path(str(month_id + 1))
    precip_days = get_monthly_num_rain_days(prcp_days_file_path)
    
    snow_factor_file_path = snow_cover_file_path(str(month_id + 1)) 
    snow_factor = calculate_snow_factor(snow_factor_file_path)
    
    wind_speed_file_path = wind_spd_file_path(mont_id_str) 

    monthly_WF = calculate_monthly_weather_factor(wind_speed_file_path, temp, precip, sol_rad, precip_days, snow_factor, pressure)
    wf_out_file_path =  weather_factor_out_file_path(str(month_id + 1))
    RasterSave(monthly_WF, wf_out_file_path, dem_file_path)
    
# Step 2 : Soil Crusting Factor and Erodibility Factor 
# # # # # # # # # # # # #
sand_ratio = preprocess_soil_nonneg(sand_file_path)
sand_ratio = sand_ratio*100

silt_ratio = preprocess_soil_nonneg(silt_file_path)
silt_ratio = silt_ratio*100

org_mat_ratio = preprocess_soil_nonneg(org_mat_file_path)
org_mat_ratio = org_mat_ratio/1000

clay_ratio = preprocess_soil_nonneg(clay_file_path)
clay_ratio = clay_ratio*100

scf = calculate_soil_crust_factor(clay_ratio,org_mat_ratio)
RasterSave(scf,scf_file_path,dem_file_path)

ef = calculate_soil_erodibility_factor(sand_ratio, silt_ratio, clay_ratio, org_mat_ratio) 
RasterSave(ef,ef_file_path, dem_file_path)


# Step 3 : Vegetation Factor and Step 4 : Surface Terrain Roughness Factor
# # # # # # # # # # # # #
for month_id in range(0, 12):
    fvc_file_path = frac_veg_cov_file_path(str(month_id+1))

    monthly_cog = calculate_vegetation_factor(fvc_file_path)
    COG_out_file_path  = cog_out_file_path(str(month_id+1))
    RasterSave(monthly_cog, COG_out_file_path, fvc_file_path)
    
    monthly_kprime = calculate_surface_terr_rough(dem_file_path, fvc_file_path)
    Kprime_out_file_path = kk_out_file_path(str(month_id+1))
    RasterSave(monthly_kprime, Kprime_out_file_path, fvc_file_path)
    


# Step 5-6 : Actual and Potential Wind Erosion
# # # # # # # # # # # # #
SL_sum = 0.0
SL_wo_veg_sum = 0.0
SL_fvc_100 = 0.0

if en_corr_plot:
    fig_corr_wf, axes_corr_wf = plt.subplots(nrows=4,ncols=3, figsize=(30, 20))
    fig_corr_ef, axes_corr_ef = plt.subplots(nrows=4,ncols=3, figsize=(30, 20))
    fig_corr_k, axes_corr_k = plt.subplots(nrows=4,ncols=3, figsize=(30, 20))
    fig_corr_scf, axes_corr_scf = plt.subplots(nrows=4,ncols=3, figsize=(30, 20))
    fig_corr_cog, axes_corr_cog = plt.subplots(nrows=4,ncols=3, figsize=(30, 20))

if en_hist_plot:
    fig_hist_wf, axes_hist_wf = plt.subplots(nrows=4,ncols=3, figsize=(30, 20))
    fig_hist_ef, axes_hist_ef = plt.subplots(nrows=4,ncols=3, figsize=(30, 20))
    fig_hist_k, axes_hist_k = plt.subplots(nrows=4,ncols=3, figsize=(30, 20))
    fig_hist_scf, axes_hist_scf = plt.subplots(nrows=4,ncols=3, figsize=(30, 20))
    fig_hist_cog, axes_hist_cog = plt.subplots(nrows=4,ncols=3, figsize=(30, 20))

for month_id in range(0, 12):
    wf_file_path = weather_factor_out_file_path(str(month_id + 1))
    wf = read_WF(wf_file_path)
    cog_file_path = cog_out_file_path(str(month_id+1))
    cog = read_COG(cog_file_path)
    fvc_100 = np.zeros((cog.shape[0], cog.shape[1]))
    fvc_100[cog == np.nan] = np.nan
    fvc_100[cog != np.nan] = np.exp(-0.00483*100)
    Kprime_file_path = kk_out_file_path(str(month_id+1))
    kprime = read_Kprime(Kprime_file_path)
    
    wind_erosion_m, wind_erosion_pot_m, wind_erosion_fvc_100 = calculate_monthly_wind_erosion(wf,ef,kprime,scf,cog, fvc_100)
    
    SL_sum +=  wind_erosion_m
    SL_wo_veg_sum += wind_erosion_pot_m
    SL_fvc_100 += wind_erosion_fvc_100
    
    if en_hist_plot:
        row = int(month_id / 3)
        col = month_id % 3
        axes_hist_wf[row, col].hist(wf.ravel(), bins= HIST_BINS)
        axes_hist_ef[row, col].set_title('wf:'+ str(month_id+1), fontdict={'fontsize':'16'});

        axes_hist_ef[row, col].hist(ef.ravel(), bins = HIST_BINS)
        axes_hist_ef[row, col].set_title('ef:'+ str(month_id+1), fontdict={'fontsize':'16'});
        
        axes_hist_k[row, col].hist(kprime.ravel(), bins = HIST_BINS)
        axes_hist_k[row, col].set_title('k:'+ str(month_id+1), fontdict={'fontsize':'16'});
        
        axes_hist_scf[row, col].hist(scf.ravel(), bins = HIST_BINS)
        axes_hist_scf[row, col].set_title('scf:'+ str(month_id+1), fontdict={'fontsize':'16'});
        
        axes_hist_cog[row, col].hist(cog.ravel(), bins = HIST_BINS)
        axes_hist_cog[row, col].set_title('cog:'+ str(month_id+1), fontdict={'fontsize':'16'}); 
        
    if en_corr_plot:
        row = int(month_id / 3)
        col = month_id % 3
        axes_corr_wf[row, col].scatter(wf,wind_erosion_m)
        plot_title = 'wf:'+ str(month_id+1)+ ', corr coeff = '+ str(compute_corr(wf,wind_erosion_m))
        axes_corr_wf[row, col].set_title(plot_title, fontdict={'fontsize':'16'}); 
        
        axes_corr_ef[row, col].scatter(ef,wind_erosion_m)
        plot_title = 'ef:'+ str(month_id+1)+ ', corr coeff = '+ str(compute_corr(ef,wind_erosion_m))
        axes_corr_ef[row, col].set_title(plot_title, fontdict={'fontsize':'16'}); 
        
        axes_corr_k[row, col].scatter(kprime,wind_erosion_m)
        plot_title = 'kprime:'+ str(month_id+1)+ ', corr coeff = '+ str(compute_corr(kprime,wind_erosion_m))
        axes_corr_k[row, col].set_title(plot_title, fontdict={'fontsize':'16'}); 
        
        axes_corr_scf[row, col].scatter(scf,wind_erosion_m)
        plot_title = 'scf:'+ str(month_id+1)+ ', corr coeff = '+ str(compute_corr(scf,wind_erosion_m))
        axes_corr_scf[row, col].set_title(plot_title, fontdict={'fontsize':'16'}); 
        
        axes_corr_cog[row, col].scatter(cog,wind_erosion_m)
        plot_title = 'cog:'+ str(month_id+1)+ ', corr coeff = '+ str(compute_corr(cog,wind_erosion_m))
        axes_corr_cog[row, col].set_title(plot_title, fontdict={'fontsize':'16'}); 

    
if en_hist_plot:
    fig_hist_wf.tight_layout()
    fig_hist_ef.tight_layout()
    fig_hist_k.tight_layout()
    fig_hist_scf.tight_layout()
    fig_hist_cog.tight_layout()
    fig_file_dir = 'data_UTM/output/histograms_plots/'
    make_directory(fig_file_dir)
    fig_file_path = os.path.join(fig_file_dir, 'wf_hist.png')
    fig_hist_wf.savefig(fig_file_path)
    fig_file_path = os.path.join(fig_file_dir, 'ef_hist.png')
    fig_hist_ef.savefig(fig_file_path)
    fig_file_path = os.path.join(fig_file_dir, 'k_hist.png')
    fig_hist_k.savefig(fig_file_path)
    fig_file_path = os.path.join(fig_file_dir, 'scf_hist.png')
    fig_hist_scf.savefig(fig_file_path)
    fig_file_path = os.path.join(fig_file_dir, 'cog_hist.png')
    fig_hist_cog.savefig(fig_file_path)
    
if en_corr_plot:
    fig_corr_wf.tight_layout()
    fig_corr_ef.tight_layout()
    fig_corr_k.tight_layout()
    fig_corr_scf.tight_layout()
    fig_corr_cog.tight_layout()
    fig_file_dir = 'data_UTM/output/correlation_plots/'
    make_directory(fig_file_dir)
    fig_file_path = os.path.join(fig_file_dir, 'wf_vs_sl.png')
    fig_corr_wf.savefig(fig_file_path)
    fig_file_path = os.path.join(fig_file_dir, 'ef_vs_sl.png')
    fig_corr_ef.savefig(fig_file_path)
    fig_file_path = os.path.join(fig_file_dir, 'k_vs_sl.png')
    fig_corr_k.savefig(fig_file_path)
    fig_file_path = os.path.join(fig_file_dir, 'scf_vs_sl.png')
    fig_corr_scf.savefig(fig_file_path)
    fig_file_path = os.path.join(fig_file_dir, 'cog_vs_sl.png')
    fig_corr_cog.savefig(fig_file_path)
plt.show()

SL_out_file_path = sl_actual_out_file_path() 
RasterSave(SL_sum, SL_out_file_path, wf_file_path)

fig_file_dir = 'data_UTM/output/histograms_plots/'
make_directory(fig_file_dir)
fig_file_path = os.path.join(fig_file_dir, 'SL_actual_hist.png')
plt.hist(SL_sum.ravel(), bins=HIST_BINS)
plt.savefig(fig_file_path)
plt.show()

SL_without_veg_out_file_path = sl_wo_veg_out_file_path()
RasterSave(SL_wo_veg_sum, SL_without_veg_out_file_path, wf_file_path)

fig_file_path = os.path.join(fig_file_dir, 'SL_wo_veg_hist.png')
plt.hist(SL_wo_veg_sum.ravel(), bins=HIST_BINS)
plt.savefig(fig_file_path)
plt.show()

SL_fvc_100_file_path = sl_fvc_100_out_file_path()
RasterSave(SL_fvc_100, SL_fvc_100_file_path, wf_file_path)
fig_file_path = os.path.join(fig_file_dir, 'SL_fvc_100_hist.png')
plt.hist(SL_fvc_100.ravel(), bins=HIST_BINS)
plt.savefig(fig_file_path)

sand_re = SL_wo_veg_sum - SL_sum
sand_re = (sand_re < 0) * 0 + (sand_re >= 0) * sand_re

sand_re_out_file_path = sand_r_out_file_path()
RasterSave(sand_re, sand_re_out_file_path, dem_file_path)
