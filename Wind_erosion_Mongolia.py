#!/usr/bin/env python
# coding: utf-8

## Wind Erosion Model
## Written by Charlie Weil & Isita Talukdar, August 2020
## Inspired by Huang Binbin's approach.

import numpy as np
import gdal
import os

# - - - - - - - - - - -
#     File paths 
# - - - - - - - - - - -

data_dir ="data_UTM/"
input_data_path = os.path.join(data_dir, 'input/')

dem_filename = 'dem_clipped.tif'
dem_file_path = os.path.join(input_data_path, dem_filename)

# Climate Input Data

def temperature_file_path(month_num):
    temperature_filename = 'temperature_degC/' + 'wc2.0_30s_tave_'+ month_num + '.tif'
    return os.path.join(input_data_path, temperature_filename)

def precipitation_file_path(month_num):
    precipitation_filename = 'precipitation_mm/' + 'wc2.0_30s_prec_' + month_num + '.tif'
    return os.path.join(input_data_path, precipitation_filename)


## Isita : do the same thing to define here, upfront, prcp_file_path, sol_file_path. prcp_days_file_path, snow_factor_file_path, wind speed, wf_out
## Be very clear about what is "filename" and "filepath". Don't overwrite these variables !!

# Soil Input Data
sand_filename = 'soil/sand.tif'
sand_file_path = os.path.join(input_data_path, sand_filename)

## Isita: same thing with other soil variables, and variables of steps 4,5,6 : TO DO !!
## .....



# - - - - - - - - - - -
#      Parameters
# - - - - - - - - - - -

Kr_WINDOW_SIZE = 3

# days of month 1-12
mondays = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
endday = [30, 58, 89, 119, 150, 180, 211, 242, 272, 303, 333, 364]
month_id = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'] ## Isita, this variable is not used, is it normal?



# - - - - - - - - - - -
#    Model functions 
# - - - - - - - - - - -

def execute(args):
    ## This is where the steps in "Execute" shoulf eventually go.
    return None

# define a class to store raster result
def RasterSave(data, path, d1):
    driver = gdal.GetDriverByName('GTiff')
    water_r_path = path
    ds = driver.Create(water_r_path, d1.RasterXSize, d1.RasterYSize, 1, gdal.GDT_Float32)
    ds.SetGeoTransform(d1.GetGeoTransform())
    ds.SetProjection(d1.GetProjection())
    ds.GetRasterBand(1).WriteArray(data)
    # ds.GetRasterBand(1).SetNoDataValue(0)
    # ds.GetRasterBand(1).DeleteNoDataValue()
    del ds
    
def read_raster_as_array(file_path):
    raster = gdal.Open(file_path)
    array = raster.ReadAsArray(0, 0, raster.RasterXSize, raster.RasterYSize)
    return array

def get_elevation_data(file_path):
    """Get Elevation Data from Raster File
    Parameters:
        file_path: name of Raster file 
        Raster File: 1 pixel = 5532m
        
    Returns:
        dem: Elevation data, in meters
        GDAL Dataset
    """

    dem = gdal.Open(file_path)
    return dem

def read_dem_as_array(dem):
    dem_a = dem.ReadAsArray(0, 0, dem.RasterXSize, dem.RasterYSize)
    #Selct only the first band
    dem_a = dem_a[0,:,:]
    dem_a2 = (dem_a < 0) * 0 + (dem_a >= 0) * dem_a
    return dem_a2

def calculate_air_pressure(dem):
    """Calculate air pressure

    Parameters:
        dem: Elevation data, in meters
        GDAL Dataset
    Returns:
        Air Pressure (P), in kPa
        stored in array 
    """
    dem_a2 = read_dem_as_array(dem)
    air_pressure = 101.3 * (1 - 0.0255 * dem_a2 / 1000 * (6357 / (6357 + dem_a2 / 1000))) ** 5.256
    return air_pressure

def get_monthly_avg_temp(file_path):
    """Get Average Temperature for a given month from Raster File

    Parameters:
        file_path: name of Raster file with average temperature
        Raster File: 1 pixel = 5532m
    Returns:
        tem_m_a2: Temperature for a given month, in °C
        stored in array 
    """
    
    # monthly average temperature(℃)
    tem_m_a = read_raster_as_array(file_path)
    # multiply 0.1 because the unit of source data is 0.1 ℃
    tem_m_a2 = (tem_m_a < -100) * 0.1 + (tem_m_a >= -100) * tem_m_a
    return tem_m_a2

def get_monthly_total_precip(file_path):
    """Get Total Precipitation for a given month from Raster File

    Parameters:
        file_path: name of Raster file with total precipitation
        Raster File: 1 pixel = 5532m

    Returns:
        prcp_a: Total Precipitation for a given month, in mm
        stored in array 
    """
    # monthly total precipitation(mm)
    prcp_a = read_raster_as_array(file_path)
    # multiply 0.1 because the unit of source data is 0.1 mm
    prcp_a2 = (prcp_a < 0) * 0 + (prcp_a >= 0) * prcp_a
    return prcp_a

def get_monthly_sol_rad(file_path):
    """Get Total Solar Radiation for a given month from Raster File

    Parameters:
        file_path: name of Raster file with total solar radiation
        Raster File: 1 pixel = 5532m
    Returns:
        SOL_a2: Total Solar Radiation for a given month, in MJ/m2
        stored in array 
    """

    # monthly total solar radiation (MJ/m2)
    SOL_a = read_raster_as_array(file_path)
    SOL_a2 = (SOL_a < 0) * 0.1 + (SOL_a >= 0) * SOL_a
    return SOL_a2

def get_monthly_num_rain_days(file_path):
    """Get Number of Rain Days for a given month from Raster File

    Parameters:
        file_path: name of Raster file with number of rain days
        Raster File: 1 pixel = 5532m
    Returns:
        prcp_days_a2: Total Number of Rain Days for a given month, in days
        stored in array 
    """

    # days of rain events in every month
    prcp_days_a = read_raster_as_array(file_path)
    prcp_days_a2 = (prcp_days_a < 0) * 0 + (prcp_days_a >= 0) * prcp_days_a
    return prcp_days_a2

def calculate_evapotranspiration(solar_rad, avg_temp):
    """Calculate the potential evapotranspiration for a given month

    Parameters:
        Temperature for a given month, in °C
        array from Raster file
        
        Total Solar Radiation for a given month, in Megajoules/meters squared
        array from Raster file
        
    Returns:
        Potential Evapotranspiration (ETp), in millimeters
        stored in array 
    """
    evap_trans = 0.0135 * (solar_rad / 2.54) * (avg_temp + 17.8)
    
    return evap_trans

def calculate_soil_moisture(solar_rad, avg_temp, precip, days_precip):
    """Calculate the Soil Moisture Factor for a given month

    Parameters:
        Potential Evapotranspiration for a given month, in millimeters
        array generated by function
        
        Number of Rain days for a given month
        array from Raster file
        
    Returns:
        Soil Moisture Factor (SW)
    """
    evap_trans = calculate_evapotranspiration(solar_rad, avg_temp)
    soil_moisture = (evap_trans - precip * days_precip) / evap_trans
    soil_moisture = (soil_moisture < 0) * 0 + (soil_moisture >= 0) * soil_moisture
    
    return soil_moisture

def calculate_snow_factor(file_path):
    """Calculate the Snow Factor for a given month

    Parameters:
        file_path: name of Raster file with Snow cover
        Raster File: 1 pixel = 5532m
    Returns:
        Snow Factor (SD)
        stored in array
    """

    SD_a = read_raster_as_array(file_path)
    SD_a2 = (SD_a < 0) * 0 + (SD_a >= 0) * SD_a
    SD_a2 = (1 - SD_a2 * 0.01)
    return SD_a2

def calculate_air_density(temperature, air_pressure):
    """Calculate air density for a given month.

    Parameters:
        temperature: Temperature for a given month, in °C
        array from Raster file
        
        air_pressure: Air Pressure, in kPa
        array from Raster file

    Returns:
        Air density (rho), in kg/m^3
        stored in array
    """
    air_density_rho = 1.293 * (273 / (273 + temperature)) * air_pressure / 101.3
    return air_density_rho

def calculate_wind_factor_daily(file_path):
    """Calculate the wind factor for a given day

    Parameters:
        file_path: name of Raster file with wind speed
        Raster File: 1 pixel = 5532m
        
    Returns:
        Wind Factor (wf)
        stored in array
    """
     # daily wind speed (m/s)
    wind_a = read_raster_as_array(file_path)
    # multiply 0.1 because the unit of source data is 0.1 m/s
    wind_a2 = (wind_a < 50) * 0 + ((wind_a >= 50) & (wind_a < 2000)) * wind_a * 0.1 + (wind_a >= 2000) * 200

    # wind factor
    wf = wind_a2 * (wind_a2 - 5) ** 2
    
    return wf

def calculate_daily_weather_factor(wind_factor_d, air_density_rho, soil_moisture, snow_factor):
    """Calculate the Weather Factor for a given day 

    Parameters:
        Wind Factor (wf)
        array from function
        
        Air density (rho), in kg/m^3
        array from function
        
        Soil Moisture Factor (SW)
        array from function
        
        Snow Factor (SD)
        array from Raster file
        
    Returns:
        Weather Factor for a given day
        stored in array
    """
    weather_factor_d = wind_factor_d * air_density_rho / 9.8 * soil_moisture * snow_factor
    return weather_factor_d

def calculate_monthly_weather_factor(wind_file_names_prefix, month_index, temp, precip, sol_rad, precip_days, snow_factor, pressure):
    """Calculate the Weather Factor for a given month

    Parameters:
        Wind Factor (wf)
        array from function
        
        Air density (rho), in kg/m^3
        array from function
        
        Soil Moisture Factor (SW)
        array from function
        
        Snow Factor (SD)
        array from Raster file
        
    Returns:
        Weather Factor for a given month (WF)
        stored in array
    """
   
    sw = calculate_soil_moisture(sol_rad, temp, precip, precip_days)

    air_density_rho = calculate_air_density(temp, pressure)
    
   
    weather_factor_m = 0.0
    #start_date = (month_index == 0) * 0 + (month_index > 0) * int(endday[month_index - 1])
    start_date = 0 if (month_index == 0) else int(endday[month_index - 1])
    end_date = int(endday[month_index]) + 1
    #Loop over every day
    for i in range(start_date, end_date):
        #wind_file_path = wind_file_names_prefix + str(i + 1) + '.tif'
        wind_file_path = wind_file_names_prefix
        wind_factor_d =  calculate_wind_factor_daily(wind_file_path)
        weather_factor_d = calculate_daily_weather_factor(wind_factor_d, air_density_rho, sw, snow_factor)
        #print(i, " : ",np.sum(snow_factor))
        weather_factor_m = weather_factor_m + weather_factor_d
        
    return weather_factor_m

def get_sand_ratio(file_path):
    """Get Sand Ratio from Raster File

    Parameters:
        file_path: name of Raster file with sand ratio
        Raster File: 1 pixel = 5532m
        
    Returns:
        sand_a2: Sand Ratio
        stored in array
    """

    #ratio of sand(%)
    sand_a = read_raster_as_array(file_path)
    sand_a2 = (sand_a <= 0) * 0.1 + (sand_a >0) * sand_a
    return sand_a2

def get_silt_ratio(file_path):
    """Get Silt Ratio from Raster File

    Parameters:
        file_path: name of Raster file with silt ratio
        Raster File: 1 pixel = 5532m
    Returns:
        silt_a2: Silt Ratio
        stored in array
    """

    #ratio of silt(%)
    silt_a = read_raster_as_array(file_path)
    silt_a2 = (silt_a <= 0) * 0.1 + (silt_a > 0) * silt_a
    return silt_a2

def get_org_mat_ratio(file_path):
    """Get Organic Matter Ratio from Raster File

    Parameters:
        file_path: name of Raster file with organic matter ratio
        Raster File: 1 pixel = 5532m

    Returns:
        om_a2: Organic Matter Ratio
        stored in array
    """
    #ratio of organic matter(%)
    om_a = read_raster_as_array(file_path)
    om_a2 = (om_a <= 0) * 0.1 + (om_a > 0) * om_a
    return om_a2

def get_clay_ratio(file_path):
    """Get Clay Ratio from Raster File

    Parameters:
        file_path: name of Raster file with clay ratio
        Raster File: 1 pixel = 5532m
        
    Returns:
        clay_a2: Clay Ratio
        stored in array
    """
    #ratio of clay(%)
    clay_a = read_raster_as_array(file_path)
    clay_a2 = (clay_a <= 0) * 0.1 + (clay_a >0) * clay_a
    return clay_a2

def calculate_soil_crust_factor(clay_ratio, org_mat_ratio):
    """Calculate the Soil Crusting Factor

    Parameters:
        Clay Ratio(%)
        array from function
        
        Organic Matter Ratio(%)
        array from function
        
    Returns:
        Soil Crusting Factor (SCF)
        stored in array
    """
    soil_crust_factor = SCF = 1 / (1 + 0.0066 * clay_ratio ** 2 + 0.021 * org_mat_ratio ** 2)

    return soil_crust_factor

def calculate_soil_erodibility_factor(sand_ratio, silt_ratio, clay_ratio, org_mat_ratio):
    """Calculate the Soil Erodibility Factor

    Parameters:
        Sand Ratio(%)
        array from function
        
        Silt Ratio(%)
        array from function
        
        Clay Ratio(%)
        array from function
        
        Organic Matter Ratio(%)
        array from function
        
    Returns:
        Soil Erodibility Factor (EF)
        stored in array
    """
    soil_erode_factor = (29.09 + 0.31 * sand_ratio + 0.17 * silt_ratio + 0.33 * sand_ratio / clay_ratio - 2.59 * org_mat_ratio - 0.95 * 0) / 100


    return soil_erode_factor

def calculate_vegetation_factor(file_path):
    """Calculate the Fraction of Monthly Vegetation Coverage For a Given Month

    Parameters:
        Fraction of Monthly Vegetation Coverage
        file_path: name of Raster file with Fraction of Monthly Vegetation Coverage
        Raster File: 1 pixel = 5532m
        
    Returns:
        Vegetation Factor(COG)
        stored in array
    """
    # read fraction of vegetation coverage data (%)
    fvc_a = read_raster_as_array(file_path)
    fvc_a2 = ((fvc_a < 0) | (fvc_a > 100)) * 0 + ((fvc_a >= 0) & (fvc_a <= 100)) * fvc_a
    vegetation_factor = np.exp(-0.00483*fvc_a2)
    
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
    temp_sum = np.zeros((rows, columns))
    deltas = np.zeros((rows, columns))
    
    # create a padded copy window only odd sized
    pad = int(window_size/2)
    matrix = np.pad(data_arr, pad, 'edge')
    window_rows = window_size
    window_cols = window_size
    # Level 2: traversing the window 
    for y in range(5):
        for x in range(5):

            # Level 1: handling the matrix 
            # (rows, columns = data.shape !)
            temp_sum = matrix[y : window_rows + y,
                              x : window_cols + x]
            deltas[x, y] = temp_sum.max() -  temp_sum.min()
    return deltas

def calculate_roughness_length(dem):
    """Calculate Roughness Length 

    Parameters:
        dem: Elevation data, in meters
        GDAL Dataset
    Returns:
        Roughness Length, in cm
        stored in array 
    """
    dem_arr = read_dem_as_array(dem)
    geotransform = dem.GetGeoTransform()
    
    deltas = window_delta(dem_arr, Kr_WINDOW_SIZE)
    resolution = (abs(geotransform[1]) + abs(geotransform[2]))/2
    roughness_length = 0.2* ((deltas **2)/(resolution*Kr_WINDOW_SIZE))
    #conversion to centimeters
    roughness_length *= 100
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
    Crr = 17.46 * (0.025 + 2.464*fvc_a2**3.56)**0.738
    return Crr
    
    
def calculate_surface_terr_rough(dem, fvc_path):
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
    Kr  = calculate_roughness_length(dem)
    Crr = calculate_chain_rand_roughness(fvc_path)
    
    KK = np.exp(1.86*Kr - 2.41*Kr**0.934 - 0.127*Crr)
    KK_a2 = (KK <= 0) * 0.1 + (KK > 0) * KK
    return KK_a2

def read_WF(file_path):
    #reads the output TIF into an array
    WF_a = read_raster_as_array(file_path)
    weather_factor = (WF_a <= 0) * 0.1 + (WF_a > 0) * WF_a
    return weather_factor

def read_COG(file_path):
    #reads the output TIF into an array
    COG_a = read_raster_as_array(file_path)
    vegetation_factor = (COG_a <= 0) * 0.1 + (COG_a > 0) * COG_a
    return vegetation_factor

def read_Kprime(file_path):
    #reads the output TIF into an array
    KK_a = read_raster_as_array(file_path)
    return KK_a


def calculate_monthly_wind_erosion(weather_factor,soil_erode_factor, kprime, soil_crust_factor, cog):
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

    return wind_erosion, wind_erosionp


# - - - - - - - - - - -
#        Execute
# - - - - - - - - - - -
# # # # # # # # # # # # # # # # #

# Step 1 : Weather Factor
# # # # # # # # # # # # #

#get elevation data
dem = get_elevation_data(dem_file_path)

#calculate air pressure
pressure = calculate_air_pressure(dem)

#Compute Monthly Weather Factor
for k in range(0, 12):
    num = str(k+1)
    if k+1 < 10:
        num = '0'+str(k+1)
    
    temp_file_path = temperature_file_path(num)
    temp = get_monthly_avg_temp(temp_file_path)
    
    prcp_file_path = 'precipitation_mm/' + 'wc2.0_30s_prec_' + num + '.tif' ## Isita, this path is to be defined, with a function, at the begining of the script in the "File paths" section ! 
    prcp_file_path = os.path.join(input_data_path, prcp_file_path)
    precip = get_monthly_total_precip(prcp_file_path)
    
    sol_file_path = 'solar_radiation/'+ 'shortwave_radiation_' + str(k + 1) + '.tif' ## Isita, idem
    sol_file_path = os.path.join(input_data_path, sol_file_path)
    sol_rad = get_monthly_sol_rad(sol_file_path)
    
    prcp_days_file_path = 'month_prcp_day/prcp_day_' + str(k + 1) + '.tif' ## Isita, idem
    prcp_days_file_path = os.path.join(input_data_path, prcp_days_file_path)
    precip_days = get_monthly_num_rain_days(prcp_days_file_path)
    
    snow_factor_file_path = 'snow_gm2/' + 'snow_' + str(k + 1) + '.tif' ## Isita, idem
    snow_factor_file_path = os.path.join(input_data_path, snow_factor_file_path)
    snow_factor = calculate_snow_factor(snow_factor_file_path)
    
    #wind_speed_file_path =  'month_wind_day/wind_'
    wind_speed_file_path = 'wind_speed_clipped.tif' ## Isita, idem
    wind_speed_file_path = os.path.join(input_data_path, wind_speed_file_path)
    monthly_WF = calculate_monthly_weather_factor(wind_speed_file_path,k, temp, precip, sol_rad, precip_days, snow_factor, pressure)
    wf_out_file_path = 'intermediate/WF/wf_' + str(k + 1) + '.tif'
    wf_out_file_path = os.path.join(data_dir, wf_out_file_path)
    RasterSave(monthly_WF, wf_out_file_path, dem)

# Step 2 : Soil Crusting Factor and Erodibility Factor
# # # # # # # # # # # # #
sand_ratio = get_sand_ratio(sand_file_path)

silt_file_path = 'soil/silt.tif'
silt_file_path = os.path.join(input_data_path, silt_file_path) ## Isita: this needs to be up in the "Files name" section, and clean - as I did for sand !
silt_ratio = get_silt_ratio(silt_file_path)

om_file_path = 'soil/soil_organic_matter_gm2.tif' ## Isita: this needs to be up in the "Files name" section, and clean - as I did for sand !
om_file_path = os.path.join(input_data_path, om_file_path)
org_mat_ratio = get_org_mat_ratio(om_file_path)
 
clay_file_path = 'soil/clay.tif' ## Isita: this needs to be up in the "Files name" section, and clean - as I did for sand !
clay_file_path = os.path.join(input_data_path, clay_file_path)
clay_ratio = get_clay_ratio(clay_file_path)

scf = calculate_soil_crust_factor(clay_ratio,org_mat_ratio) ## Isita: this needs to be up in the "Files name" section, and clean - as I did for sand !
scf_out_file_path = 'intermediate/SCF.tif'
scf_out_file_path = os.path.join(data_dir, scf_out_file_path)
RasterSave(scf,scf_out_file_path,dem)

ef = calculate_soil_erodibility_factor(sand_ratio, silt_ratio, clay_ratio, org_mat_ratio) ## Isita: this needs to be up in the "Files name" section, and clean - as I did for sand !
ef_out_file_path = 'intermediate/EF.tif'
ef_out_file_path = os.path.join(data_dir, ef_out_file_path)
RasterSave(ef,ef_out_file_path,dem)


# Step 3 : Vegetation Factor and Step 4 : Surface Terrain Roughness Factor
# # # # # # # # # # # # #
for i in range(0, 12):
    fvc_file_path = 'vegetation_percent_cover/' + 'vegetation_percent_cover_' + str(i+1) + '.tif'
    fvc_file_path = os.path.join(input_data_path, fvc_file_path)
    monthly_cog = calculate_vegetation_factor(fvc_file_path)
    COG_out_file_path  = 'intermediate/COG/COG_' + str(i + 1) + '.tif'
    COG_out_file_path = os.path.join(data_dir, COG_out_file_path)
    fvc = gdal.Open(fvc_file_path)
    RasterSave(monthly_cog, COG_out_file_path, fvc)
    
    monthly_kprime = calculate_surface_terr_rough(dem, fvc_file_path)
    Kprime_out_file_path  = 'intermediate/KK/KK_' + str(i + 1) + '.tif'
    Kprime_out_file_path = os.path.join(data_dir, Kprime_out_file_path)
    RasterSave(monthly_kprime, Kprime_out_file_path, fvc)
    


# Step 5-6 : Actual and Potential Wind Erosion
# # # # # # # # # # # # #
SL_sum = 0.0
SL_p_sum = 0.0


for i in range(0, 12):
    wf_file_path = 'intermediate/WF/WF_' + str(i + 1) + '.tif'
    wf_file_path = os.path.join(data_dir, wf_file_path)
    wf = read_WF(wf_file_path)
    
    cog_file_path ='intermediate/COG/COG_' + str(i + 1) + '.tif'
    cog_file_path = os.path.join(data_dir, cog_file_path)
    cog = read_COG(cog_file_path)
    
    Kprime_file_path  = 'intermediate/KK/KK_' + str(i + 1) + '.tif'
    Kprime_file_path = os.path.join(data_dir, Kprime_file_path)
    kprime = read_Kprime(Kprime_file_path)
    
    wind_erosion_m, wind_erosion_pot_m = calculate_monthly_wind_erosion(wf,ef,kprime,scf,cog)
    SL_sum +=  wind_erosion_m
    SL_p_sum += wind_erosion_pot_m

WF_format = gdal.Open(wf_file_path)
SL_out_file_path ='output/SL_actual.tif'
SL_out_file_path = os.path.join(data_dir, SL_out_file_path)
RasterSave(SL_sum, SL_out_file_path, WF_format)

SL_p_out_file_path ='output/SL_without_veg.tif'
SL_p_out_file_path = os.path.join(data_dir, SL_p_out_file_path)
RasterSave(SL_p_sum, SL_p_out_file_path, WF_format)

sand_re = SL_p_sum - SL_sum
sand_re = (sand_re < 0) * 0 + (sand_re >= 0) * sand_re

sand_re_out_file_path = 'output/sand_re.tif'
sand_re_out_file_path = os.path.join(data_dir, sand_re_out_file_path)
RasterSave(sand_re, sand_re_out_file_path, dem)



# In[ ]:




