#!/usr/bin/env python
# coding: utf-8

# In[33]:


#!/usr/bin/env python
# coding: utf-8

## Wind Erosion Model
## Written by Charlie Weil & Isita Talukdar, August 2020
## Inspired by Huang Binbin's approach.
import numpy as np
import pandas as pd
import numpy.ma as ma
import gdal
import os
import matplotlib.pyplot as plt #visualisation
import warnings
warnings.filterwarnings('ignore')

import sys, getopt

save_output_to_file = False
dopltshow = False
def make_directory(path):
    if not os.path.exists(path):
        os.makedirs(path)
        
# - - - - - - - - - - -
#     File paths 
# - - - - - - - - - - -
data_dir ="Data/"
input_dir = "Inputs/"


sl_output_suffix = 'none'

argv = sys.argv[1:]

input_dir_given = False

dem_dir = ""
dem_file_name = ""
temper_dir = ""
temper_prefix = ""
precip_dir = ""
precip_prefix = ""
sol_dir = ""
sol_prefix = ""
prcp_days_dir = ""
prcp_days_prefix = ""
snow_dir = ""
snow_prefix = ""
wind_speed_dir = ""
wind_speed_prefix = ""
fvc_dir = ""
fvc_prefix = ""
sand_dir = ""
sand_file_name = ""
silt_dir = ""
silt_file_name = ""
clay_dir = ""
clay_file_name = ""
som_dir = ""
som_file_name = ""

try:
    opts, args = getopt.getopt(argv, "h",["data_dir=", "input_data_dir=","dem_dir=", "dem_file_name=",                                           "temper_dir=", "temper_prefix=",                                          "precip_dir=", "precip_prefix=",                                           "sol_dir=", "sol_prefix=",                                         "prcp_days_dir=","prcp_days_prefix=",                                          "snow_dir=", "snow_prefix=",                                          "wind_speed_dir=", "wind_speed_prefix=",                                          "sand_dir=", "sand_file_name=",                                          "silt_dir=","silt_file_name=",                                          "clay_dir=", "clay_file_name=",                                          "som_dir=", "som_file_name=",                                         "fvc_dir=","fvc_prefix=", "sl_output_suffix="                                         ])
except getopt.GetoptError:
    print('mongolia_erosion.py -h --data_dir=<working directory> --dem_file_name=<DEM file name>')
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print('mongolia_erosion.py -h --data_dir=<working directory>')
        sys.exit()
    elif opt in ("--data_dir"):
        data_dir = arg
    elif opt in ("--input_data_dir"):
        input_dir = arg
        input_dir_given = True
    elif opt in ("--dem_dir"):
        dem_dir = arg
    elif opt in ("--dem_file_name"):
        dem_file_name = arg
    elif opt in ("--temper_dir"):
        temper_dir = arg 
    elif opt in ("--temper_prefix"):
        temper_prefix = arg
    elif opt in ("--precip_dir"):
        precip_dir = arg 
    elif opt in ("--precip_prefix"):
        precip_prefix = arg 
    elif opt in ("--sol_dir"):
        sol_dir = arg 
    elif opt in ("--sol_prefix"):
        sol_prefix = arg 
    elif opt in ("--prcp_days_dir"):
        prcp_days_dir = arg 
    elif opt in ("--prcp_days_prefix"):
        prcp_days_prefix = arg 
    elif opt in ("--snow_dir"):
        snow_dir = arg 
    elif opt in ("--snow_prefix"):
        snow_prefix = arg 
    elif opt in ("--wind_speed_dir"):
        wind_speed_dir = arg 
    elif opt in ("--wind_speed_prefix"):
        wind_speed_prefix = arg  
    elif opt in ("--sand_dir"):
        sand_dir = arg 
    elif opt in ("--sand_file_name"):
        sand_file_name = arg 
    elif opt in ("--silt_dir"):
        silt_dir = arg 
    elif opt in ("--silt_file_name"):
        silt_file_name = arg 
    elif opt in ("--clay_dir"):
        clay_dir = arg 
    elif opt in ("--clay_file_name"):
        clay_file_name = arg 
    elif opt in ("--som_dir"):
        som_dir = arg 
    elif opt in ("--som_file_name"):
        som_file_name = arg
    elif opt in ("--fvc_dir"):
        fvc_dir = arg 
    elif opt in ("--fvc_prefix"):
        fvc_prefix = arg
    elif opt in ("--sl_output_suffix"):
        sl_output_suffix = arg
        
        
 
if input_dir_given:
    if dem_dir != "":
        print("DEM directory will be ignored. Input directory has been specified as: ", input_dir) 
    if dem_file_name != "":
        print("DEM filename will be ignored. Input directory has been specified as: ", input_dir) 
    if temper_dir != "":
        print("Temperature directory will be ignored. Input directory has been specified as: ", input_dir) 
    if temper_prefix != "":
        print("Temperature filename will be ignored. Input directory has been specified as: ", input_dir) 
    if precip_dir != "":
        print("Precipitation directory will be ignored. Input directory has been specified as: ", input_dir) 
    if precip_prefix != "":
        print("Precipitation filename will be ignored. Input directory has been specified as: ", input_dir) 
    if sol_dir != "":
        print("Solar Radiation directory will be ignored. Input directory has been specified as: ", input_dir) 
    if sol_prefix != "":
        print("Solar Radiation filename will be ignored. Input directory has been specified as: ", input_dir) 
    if prcp_days_dir != "":
        print("Precipitation Days directory will be ignored. Input directory has been specified as: ", input_dir) 
    if prcp_days_prefix != "":
        print("Precipitation Days filename will be ignored. Input directory has been specified as: ", input_dir) 
    if snow_dir != "":
        print("Snow directory will be ignored. Input directory has been specified as: ", input_dir) 
    if snow_prefix != "":
        print("Snow filename will be ignored. Input directory has been specified as: ", input_dir) 
    if wind_speed_dir != "":
        print("Wind speed directory will be ignored. Input directory has been specified as: ", input_dir) 
    if wind_speed_prefix != "":
        print("Wind speed filename will be ignored. Input directory has been specified as: ", input_dir) 
    if fvc_dir != "":
        print("Fractional Vegetation Coverage directory will be ignored. Input directory has been specified as: ", input_dir) 
    if fvc_prefix != "":
        print("Fractional Vegetation Coverage filename will be ignored. Input directory has been specified as: ", input_dir) 
    if sand_dir != "":
        print("Sand  directory will be ignored. Input directory has been specified as: ", input_dir) 
    if sand_file_name != "":
        print("Sand  filename will be ignored. Input directory has been specified as: ", input_dir) 
    if silt_dir != "":
        print("Silt directory will be ignored. Input directory has been specified as: ", input_dir) 
    if silt_file_name != "":
        print("Silt filename will be ignored. Input directory has been specified as: ", input_dir) 
    if clay_dir != "":
        print("Clay directory will be ignored. Input directory has been specified as: ", input_dir) 
    if clay_file_name != "":
        print("Clay filename will be ignored. Input directory has been specified as: ", input_dir) 
    if som_dir != "":
        print("SOM directory will be ignored. Input directory has been specified as: ", input_dir) 
    if som_file_name != "":
        print("SOM filename will be ignored. Input directory has been specified as: ", input_dir) 

    dem_dir = input_dir
    dem_file_name = "dem_aoi.tif"
    temper_dir = os.path.join(input_dir, "temperature_degC/")
    temper_prefix = "wc2.0_30s_tave_"
    precip_dir = os.path.join(input_dir, "precipitation_mm/")
    precip_prefix = "wc2.0_30s_prec_"
    sol_dir = os.path.join(input_dir, "solar_radiation/")
    sol_prefix = "shortwave_radiation_"
    prcp_days_dir = os.path.join(input_dir, "month_prcp_day/")
    prcp_days_prefix = "prcp_day_"
    snow_dir = os.path.join(input_dir, "snow_gm2/")
    snow_prefix = "snow_"
    wind_speed_dir = os.path.join(input_dir, "wind_speed_daily/")
    wind_speed_prefix = "utm_daily_wind_"
    sand_dir = os.path.join(input_dir, "soil/")
    sand_file_name = "sand.tif"
    silt_dir = os.path.join(input_dir, "soil/")
    silt_file_name = "silt.tif"
    clay_dir = os.path.join(input_dir, "soil/")
    clay_file_name = "clay.tif"
    som_dir = os.path.join(input_dir, "soil/")
    som_file_name = "soil_organic_matter_gm2.tif"
    fvc_dir = os.path.join(input_dir, "vegetation_percent_cover/")
    fvc_prefix = "vegetation_percent_cover_"
else:
    missing = False
    if dem_dir == "":
        print("DEM directory missing, exiting...")
        missing = True
    if dem_file_name == "":
        print("DEM filename missing, exiting...")
        missing = True
    if temper_dir == "":
        print("Temperature directory missing, exiting...")
        missing = True
    if temper_prefix == "":
        print("Temperature filename missing, exiting...")
        missing = True
    if precip_dir == "":
        print("Precipitation directory missing, exiting...")
        missing = True
    if precip_prefix == "":
        print("Precipitation filename missing, exiting...")
        missing = True
    if sol_dir == "":
        print("Solar Radiation directory missing, exiting...")
        missing = True
    if sol_prefix == "":
        print("Solar Radiation filename missing, exiting...")
        missing = True
    if prcp_days_dir == "":
        print("Precipitation Days directory missing, exiting...")
        missing = True
    if prcp_days_prefix == "":
        print("Precipitation Days filename missing, exiting...")
        missing = True
    if snow_dir == "":
        print("Snow directory missing, exiting...")
        missing = True
    if snow_prefix == "":
        print("Snow filename missing, exiting...")
        missing = True
    if wind_speed_dir == "":
        print("Wind speed directory missing, exiting...")
        missing = True
    if wind_speed_prefix == "":
        print("Wind speed filename missing, exiting...")
        missing = True
    if fvc_dir == "":
        print("Fractional Vegetation Coverage directory missing, exiting...")
        missing = True
    if fvc_prefix == "":
        print("Fractional Vegetation Coverage filename missing, exiting...")
        missing = True
    if sand_dir == "":
        print("Sand  directory missing, exiting...")
        missing = True
    if sand_file_name == "":
        print("Sand  filename missing, exiting...")
        missing = True
    if silt_dir == "":
        print("Silt directory missing, exiting...")
        missing = True
    if silt_file_name == "":
        print("Silt filename missing, exiting...")
        missing = True
    if clay_dir == "":
        print("Clay directory missing, exiting...")
        missing = True
    if clay_file_name == "":
        print("Clay filename missing, exiting...")
        missing = True
    if som_dir == "":
        print("SOM directory missing, exiting...")
        missing = True
    if som_file_name == "":
        print("SOM filename missing, exiting...")
        missing = True

    if ( missing ):
        sys.exit(2)

print('########################')
print("\t * "+temper_dir+temper_prefix+ "<n>.tif: monthly average temperature in deg C ")
print("\t * "+precip_dir+precip_prefix+ "<n>.tif: monthly total precipitation in mm ")
print("\t * "+sol_dir+sol_prefix+ "<n>.tif: monthly solar radiation in MJ/m^2 ")
print("\t * "+prcp_days_dir+prcp_days_prefix+ "<n>.tif: monthly number of rain days  ")
print("\t * "+snow_dir+snow_prefix+ "<n>.tif: monthly snow cover factor (probability) ")
print("\t * "+wind_speed_dir+wind_speed_prefix+ "<n>.tif: daily average wind speed in m/s ")
print("\t * "+fvc_dir+fvc_prefix+ "<n>.tif: monthly fractional vegetation coverage in % ")
print("\t * "+sand_dir+sand_file_name+ " : non-temporal sand ratio (filename is a TIF file including extension) ")
print("\t * "+silt_dir+silt_file_name+ " : non-temporal silt ratio (filename is a TIF file including extension) ")
print("\t * "+clay_dir+clay_file_name+ " : non-temporal clay ratio (filename is a TIF file including extension) ")
print("\t * "+som_dir+som_file_name+ " : non-temporal organic matter ratio, data expected to be percentage*100 (filename is a TIF file including extension) ")

print('########################')
print("Creating all outputs in:", data_dir)

intermediate_data_path = os.path.join(data_dir, 'Intermediate/')
make_directory(intermediate_data_path)
output_data_path = os.path.join(data_dir, 'Outputs/')
make_directory(output_data_path)


if(save_output_to_file):
    result_file_name = os.path.join(output_data_path, 'out.txt')
    orig_stdout = sys.stdout
    f = open(result_file_name, 'w')
    sys.stdout = f


#DEM
dem_file_path = os.path.join(dem_dir, dem_file_name)

# Climate Input Data

def temperature_file_path(month_num):
    temperature_filename = temper_prefix + month_num + '.tif'
    return os.path.join(temper_dir, temperature_filename)

def precipitation_file_path(month_num):
    precipitation_filename = precip_prefix + month_num + '.tif'
    return os.path.join(precip_dir, precipitation_filename)

def solar_rad_file_path(month_num):
    solar_rad_filename = sol_prefix + month_num + '.tif'
    return os.path.join(sol_dir, solar_rad_filename)

def rain_days_file_path(month_num):
    prcp_days_filename = prcp_days_prefix + month_num + '.tif'
    return os.path.join(prcp_days_dir, prcp_days_filename)

def snow_cover_file_path(month_num):
    snow_cover_filename = snow_prefix + month_num + '.tif'
    return os.path.join(snow_dir, snow_cover_filename)

def wind_spd_file_path(day_num):
    wind_speed_filename = wind_speed_prefix + day_num + '.tif'
    return os.path.join(wind_speed_dir, wind_speed_filename)

def weather_factor_out_file_path(month_num):
    wf_file_path = intermediate_data_path + 'WF/'
    make_directory(wf_file_path)
    wf_out_filename = 'wf_' + month_num + '.tif'
    return os.path.join(wf_file_path, wf_out_filename)

# Soil Input Data
sand_file_path = os.path.join(sand_dir, sand_file_name)

silt_file_path = os.path.join(silt_dir, silt_file_name)

clay_file_path = os.path.join(clay_dir, clay_file_name)

org_mat_file_path = os.path.join(som_dir, som_file_name)

scf_filename = 'SCF.tif'
scf_file_path = os.path.join(intermediate_data_path, scf_filename)

ef_filename = 'EF.tif'
ef_file_path = os.path.join(intermediate_data_path, ef_filename)

#Vegetation Factor and K' Input Data
def frac_veg_cov_file_path(month_num):
    frac_veg_cov_filename = fvc_prefix + month_num + '.tif'
    return os.path.join(fvc_dir, frac_veg_cov_filename)

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
    sl_actual_out_filename = 'SL_'+ sl_output_suffix + '.tif'
    return os.path.join(output_data_path, sl_actual_out_filename)

def sl_actual_monthly_out_file_path(month_num): 
    sl_actual_out_filename = 'SL_'+sl_output_suffix+'_' + str(month_num + 1)+ '.tif'
    return os.path.join(output_data_path, sl_actual_out_filename)

def sl_wo_veg_out_file_path():
    sl_wo_veg_out_filename = 'SL_without_veg.tif'
    return os.path.join(output_data_path, sl_wo_veg_out_filename)

def sl_percent_diff_out_file_path():
    sl_percent_diff_out_filename = 'SL_percent_diff.tif'
    return os.path.join(output_data_path, sl_percent_diff_out_filename)

def sl_fvc_100_out_file_path():
    sl_fvc_100_out_filename = 'SL_fvc_100.tif'
    return os.path.join(output_data_path, sl_fvc_100_out_filename)

def sl_fvc_50_out_file_path():
    sl_fvc_50_out_filename = 'SL_fvc_50.tif'
    return os.path.join(output_data_path, sl_fvc_50_out_filename)

def sl_fvc_0_out_file_path():
    sl_fvc_0_out_filename = 'SL_fvc_0.tif'
    return os.path.join(output_data_path, sl_fvc_0_out_filename)

def sl_fvc_10p_out_file_path():
    sl_fvc_10p_out_filename = 'SL_fvc_10p.tif'
    return os.path.join(output_data_path, sl_fvc_10p_out_filename)

def sl_fvc_20p_out_file_path():
    sl_fvc_20p_out_filename = 'SL_fvc_20p.tif'
    return os.path.join(output_data_path, sl_fvc_20p_out_filename)

def sand_r_out_file_path():
    sand_re_out_filename = 'sand_re.tif'
    return os.path.join(output_data_path, sand_re_out_filename)

# - - - - - - - - - - -
#      Parameters
# - - - - - - - - - - -
HIST_BINS = 50
Kr_WINDOW_SIZE = 3
en_corr_plot = False
en_hist_plot = False
min_ind = -5.0
max_ind = -10.0
# days of month 1-12
mondays = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
endday = [30, 58, 89, 119, 150, 180, 211, 242, 272, 303, 333, 364]



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
    ndv = -9999.0
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
    if temperature_factor == min_ind:
        tem_m = temperature_min
    elif temperature_factor == max_ind:
        tem_m = temperature_max
    else:
        tem_m = tem_m * temperature_factor
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
    if precip_factor == min_ind:
        prcp = precip_min
    elif precip_factor == max_ind:
        prcp = precip_max
    else:
        prcp = prcp * precip_factor
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
    
    if solar_rad_factor == min_ind:
        SOL = solar_rad_min
    elif solar_rad_factor == max_ind:
        SOL = solar_rad_max
    else:
        SOL = SOL * solar_rad_factor
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
    if prcp_days_factor == min_ind:
        prcp_days = prcp_days_min
    elif prcp_days_factor == max_ind:
        prcp_days = prcp_days_max
    else:
        prcp_days = prcp_days* prcp_days_factor
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
    if snow_cover_factor == min_ind:
        SD = 0
    elif snow_cover_factor == max_ind:
        SD = 0
    else:
        SD = SD * snow_cover_factor
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

def calculate_wind_factor_monthly(month_id):
    """Calculate the wind factor for a given month
    Parameters:
        file_path: path to Raster file with monthly wind speed
        
    Returns:
        Wind Factor (wf)
        GDAL array
    """
    if month_id == 0:
        start_date = 0
        end_date = int(endday[month_id]) + 1
    else:
        start_date = int(endday[month_id - 1]) + 1
        end_date = int(endday[month_id]) + 1
    
    wf_sum = 0.0
    counter = 0
    for i in range(start_date, end_date):
     # daily wind speed (m/s)
        file_path = wind_spd_file_path(str(i+1))
        wind_speed = read_raster_as_array(file_path)
        wind_speed = np.absolute(wind_speed)
        if wind_speed_factor == min_ind:
            wind_speed = wind_speed_min
        elif wind_speed_factor == max_ind:
            wind_speed = wind_speed_max
        else:
            wind_speed = wind_speed * wind_speed_factor 
    
        #Forcing wind speed to zero if it is less than the RWEQ threshold of 5m/s
        wind_speed = (wind_speed < 5) * 0 + (wind_speed >= 5) * wind_speed
    
        # wind factor
        wf = wind_speed * (wind_speed - 5) ** 2
        wf_sum = wf_sum + wf
        counter = counter + 1
    return wf_sum

def calculate_monthly_weather_factor(wind_file_path, temp, precip, sol_rad, precip_days, snow_factor, pressure, month_id):
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
   
    wind_factor =  calculate_wind_factor_monthly(month_id)

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
def preprocess_veg_coverage(file_path, fvc_fact):
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
    fvc1 = fvc
    if fvc_fact == min_ind:
        fvc1[np.invert(np.isnan(fvc))] = fvc_min
    elif fvc_fact == max_ind:
        fvc1[np.invert(np.isnan(fvc))] = fvc_max
    else:
        fvc1 = fvc

    return fvc1

def calculate_vegetation_factor(file_path, fvc_factor):
    """Calculate the Fraction of Monthly Vegetation Coverage For a Given Month
    Parameters:
        Fraction of Monthly Vegetation Coverage
        file_path: path to Raster file with Fraction of Monthly Vegetation Coverage
        
    Returns:
        Vegetation Factor(COG)
        GDAL array
    """
    # read fraction of vegetation coverage data (%)
    fvc = preprocess_veg_coverage(file_path, fvc_factor)
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
    #convert to kilometers
    deltas = deltas/1000
    resolution = (abs(geotransform[1]) + abs(geotransform[5]))/2
    #convert to kilometers
    resolution = resolution/1000
    roughness_length = 0.2* ((deltas **2)/(resolution*(Kr_WINDOW_SIZE-1)))
    return roughness_length
    
def calculate_chain_rand_roughness(fvc_path, fvc_factor):
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
    fvc = ((fvc_a < 0) | (fvc_a > 100)) * 0 + ((fvc_a >= 0) & (fvc_a <= 100)) * fvc_a
    
    if fvc_factor == min_ind:
        fvc = fvc_min
    elif fvc_factor == max_ind:
        fvc = fvc_max
    else:
        fvc = fvc*fvc_factor
    fvc = fvc/100
    Crr = 17.46 * (0.025 + 2.464*fvc**3.56)**0.738
    return Crr
    
    
def calculate_surface_terr_rough(dem_file_path, fvc_path, fvc_factor):
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
    Crr = calculate_chain_rand_roughness(fvc_path, fvc_factor)
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
    #Wind Erosion 
    Qmax = 109.8 * (weather_factor * soil_erode_factor * kprime * soil_crust_factor * cog)
    s = 105.71 * (weather_factor * soil_erode_factor * kprime * soil_crust_factor * cog)**-0.3711
    wind_erosion = 100/(s * s +0.01)*Qmax * (np.exp(-(50 / (s + 0.01)) ** 2))
    

    
    return wind_erosion

def compute_corr(arr_x, arr_y):
    a=ma.masked_invalid(arr_x)
    b=ma.masked_invalid(arr_y)
    msk = (~a.mask & ~b.mask)
    return round(ma.corrcoef(a[msk],b[msk])[0,1], 3)

def corr_max_min(arr_x, arr_y, corr_min, corr_max, old_min_index, old_max_index, month_id):
    corr = compute_corr(arr_x,arr_y)
    if corr > corr_max:
        max = corr
        max_index = month_id
    else: 
        max = corr_max
        max_index = old_max_index
    if corr < corr_min:
        min = corr
        min_index = month_id
    else:
        min = corr_min
        min_index = old_min_index
    return min, max, min_index, max_index,  

# - - - - - - - - - - -
#        Execute
# - - - - - - - - - - -
# # # # # # # # # # # # # # # # #

master_name_list = ['wind_speed', 'temperature', 'precip', 'som', 'fvc', 'solar_rad', 'snow_cover', 'sand', 'silt', 'clay', 'primary']
scale_list = [1.00]
factor_list = [0]
SL_sum_0p5_Exists = False
SL_sum_1p5_Exists = False
doSensitivity = False

for factor_select in factor_list:
    SL_avg_arr = np.zeros(len(scale_list))
    SL_wo_veg_avg_arr = np.zeros(len(scale_list))
    SL_fvc_100_avg_arr = np.zeros(len(scale_list))
    SL_fvc_50_avg_arr = np.zeros(len(scale_list))
    SL_fvc_0_avg_arr = np.zeros(len(scale_list))
    SL_fvc_10p_avg_arr = np.zeros(len(scale_list))
    SL_fvc_20p_avg_arr = np.zeros(len(scale_list))
    fact_ind = 0
    for scale_factor in scale_list:
        wind_speed_factor = 1.0
        temperature_factor = 1.0
        precip_factor = 1.0
        som_factor = 1.0
        fvc_factor = 1.0
        solar_rad_factor = 1.0
        snow_cover_factor = 1.0
        sand_factor = 1.0
        silt_factor = 1.0
        clay_factor = 1.0
        primary_factor = 1.0
        prcp_days_factor = 1.0
        if(factor_select == 0):
            wind_speed_factor = scale_factor 
        if(factor_select == 1):
            temperature_factor = scale_factor 
        if(factor_select == 2):
            precip_factor = scale_factor
        if(factor_select == 3):
            som_factor = scale_factor
        if(factor_select == 4):
            fvc_factor = scale_factor
        if(factor_select == 5):
            solar_rad_factor = scale_factor
        if(factor_select == 6):
            snow_cover_factor = scale_factor
        if(factor_select == 7):
            sand_factor = scale_factor
        if(factor_select == 8):
            silt_factor = scale_factor
        if(factor_select == 9):
            clay_factor = scale_factor
        if(factor_select == 10):
            if scale_factor == min_ind or scale_factor == max_ind:
                primary_factor = 1.0
            else:
                primary_factor = scale_factor   
        if(factor_select == 11):
            prcp_days_factor = scale_factor
            
        # Step 1 : Weather Factor
        # # # # # # # # # # # # #

        #calculate air pressure
        pressure = calculate_air_pressure(dem_file_path)

        #Compute Monthly Weather Factor 
        for month_id in range(0, 12):
            month_id_str = str(month_id+1)

            temp_file_path = temperature_file_path(month_id_str)
            temp = get_monthly_avg_temp(temp_file_path)

            prcp_file_path = precipitation_file_path(month_id_str)  
            precip = get_monthly_total_precip(prcp_file_path)

            sol_file_path = solar_rad_file_path(month_id_str)
            sol_rad = get_monthly_sol_rad(sol_file_path)

            prcp_days_file_path = rain_days_file_path(month_id_str)
            precip_days = get_monthly_num_rain_days(prcp_days_file_path)

            snow_factor_file_path = snow_cover_file_path(month_id_str) 
            snow_factor = calculate_snow_factor(snow_factor_file_path)

            wind_speed_file_path = wind_spd_file_path(month_id_str) 

            monthly_WF = calculate_monthly_weather_factor(wind_speed_file_path, temp, precip, sol_rad, precip_days, snow_factor, pressure, month_id)
            wf_out_file_path =  weather_factor_out_file_path(month_id_str)
            RasterSave(monthly_WF, wf_out_file_path, dem_file_path)

        # Step 2 : Soil Crusting Factor and Erodibility Factor 
        # # # # # # # # # # # # #
        sand_ratio = preprocess_soil_nonneg(sand_file_path)

        if sand_factor == min_ind:
            sand_ratio = sand_min
        elif sand_factor == max_ind:
            sand_ratio = sand_max
        else:
            sand_ratio = sand_ratio*sand_factor
        sand_ratio = sand_ratio*100

        silt_ratio = preprocess_soil_nonneg(silt_file_path)
        if silt_factor == min_ind:
            silt_ratio = silt_min
        elif silt_factor == max_ind:
            silt_ratio = silt_max
        else:
            silt_ratio = silt_ratio*silt_factor
        silt_ratio = silt_ratio*100

        org_mat_ratio = preprocess_soil_nonneg(org_mat_file_path)
        if som_factor == min_ind:
            org_mat_ratio = som_min
        elif som_factor == max_ind:
            org_mat_ratio = som_max
        else:
            org_mat_ratio = org_mat_ratio*som_factor
        org_mat_ratio = org_mat_ratio/1000

        clay_ratio = preprocess_soil_nonneg(clay_file_path)
        if clay_factor == min_ind:
            clay_ratio = clay_min
        elif clay_factor == max_ind:
            clay_ratio = clay_max
        else:
            clay_ratio = clay_ratio*clay_factor 
        clay_ratio = clay_ratio*100

        scf = calculate_soil_crust_factor(clay_ratio,org_mat_ratio)
        RasterSave(scf,scf_file_path,dem_file_path)

        ef = calculate_soil_erodibility_factor(sand_ratio, silt_ratio, clay_ratio, org_mat_ratio) 
        RasterSave(ef,ef_file_path, dem_file_path)


        # Step 3 : Vegetation Factor and Step 4 : Surface Terrain Roughness Factor
        # # # # # # # # # # # # #
        for month_id in range(0, 12):
            month_id_str = str(month_id+1)
            fvc_file_path = frac_veg_cov_file_path(month_id_str)

            monthly_cog = calculate_vegetation_factor(fvc_file_path, 1.0)

            COG_out_file_path  = cog_out_file_path(month_id_str)
            RasterSave(monthly_cog, COG_out_file_path, fvc_file_path)
            monthly_kprime = calculate_surface_terr_rough(dem_file_path, fvc_file_path, 1.0)
            Kprime_out_file_path = kk_out_file_path(month_id_str)
            RasterSave(monthly_kprime, Kprime_out_file_path, fvc_file_path)

        fvc_file_path = frac_veg_cov_file_path(str(100)) 
        kprime_fvc_100 = calculate_surface_terr_rough(dem_file_path, fvc_file_path, 1.0)
        Kprime_out_fvc_100_file_path = kk_out_file_path(str(100))
        RasterSave(kprime_fvc_100, Kprime_out_fvc_100_file_path, fvc_file_path)

        cog_fvc_100 = calculate_vegetation_factor(fvc_file_path, 1.0)
        
        cog_out_fvc_100_file_path = cog_out_file_path(str(100))
        RasterSave(cog_fvc_100, cog_out_fvc_100_file_path, fvc_file_path)

        # Step 5-6 : Actual and Potential Wind Erosion
        # # # # # # # # # # # # #
        SL_sum = 0.0
        SL_wo_veg_sum = 0.0
        SL_fvc_100 = 0.0
        SL_fvc_50 = 0.0
        SL_fvc_0 = 0.0
        SL_fvc_10p = 0.0
        SL_fvc_20p  = 0.0
        


        indices = pd.Index(["fvc", "cog", "wf", "ef", "scf", "kprime", "SL"])
        df_min = pd.DataFrame(data=np.tile(2.0, (len(indices), len(indices))), index=indices, columns=indices)
        df_max = pd.DataFrame(data=np.tile(-2.0, (len(indices),len(indices))), index=indices, columns=indices)
        
        df_min_index = pd.DataFrame(data=np.tile(0, (len(indices), len(indices))), index=indices, columns=indices)
        df_max_index = pd.DataFrame(data=np.tile(0, (len(indices),len(indices))), index=indices, columns=indices)

        
        for month_id in range(0, 12):
            month_id_str = str(month_id+1)
            wf_file_path = weather_factor_out_file_path(month_id_str)
            wf = read_WF(wf_file_path)
            
            fvc_file_path = frac_veg_cov_file_path(month_id_str)
            cog = calculate_vegetation_factor(fvc_file_path, fvc_factor)
            kprime = calculate_surface_terr_rough(dem_file_path, fvc_file_path, fvc_factor)

            fvc_file_path = frac_veg_cov_file_path(month_id_str)
            cog_10p = calculate_vegetation_factor(fvc_file_path, 1.10)
            kprime_10p = calculate_surface_terr_rough(dem_file_path, fvc_file_path, fvc_factor)

            fvc_file_path = frac_veg_cov_file_path(month_id_str)
            cog_20p = calculate_vegetation_factor(fvc_file_path, 1.20)
            kprime_20p = calculate_surface_terr_rough(dem_file_path, fvc_file_path, fvc_factor)

            fvc_file_path = frac_veg_cov_file_path(str(100))
            cog_100 = calculate_vegetation_factor(fvc_file_path, 1.0)
            kprime_100 = calculate_surface_terr_rough(dem_file_path, fvc_file_path, fvc_factor)

            fvc_file_path = frac_veg_cov_file_path(str(50))
            cog_50 = calculate_vegetation_factor(fvc_file_path, 1.0)
            kprime_50 = calculate_surface_terr_rough(dem_file_path, fvc_file_path, fvc_factor)

            #compute actual wind erosion
            wind_erosion_m = calculate_monthly_wind_erosion(primary_factor*wf,ef,kprime,scf,cog)
            #compute potential wind erosion
            cog_1 = cog_50/cog_50
            wind_erosion_pot_m = calculate_monthly_wind_erosion(primary_factor*wf,ef,kprime,scf,cog_1)
            #compute 100% fvc wind erosion
            wind_erosion_fvc_100 = calculate_monthly_wind_erosion(primary_factor*wf,ef,kprime_100,scf,cog_100)
            #compute 50% fvc wind erosion
            wind_erosion_fvc_50 = calculate_monthly_wind_erosion(primary_factor*wf,ef,kprime_50,scf,cog_50)
            #compute +10% fvc wind erosion
            wind_erosion_fvc_10p = calculate_monthly_wind_erosion(primary_factor*wf,ef,kprime_10p,scf,cog_10p)
            #compute +20% fvc wind erosion
            wind_erosion_fvc_20p = calculate_monthly_wind_erosion(primary_factor*wf,ef,kprime_20p,scf,cog_20p)
            
            SL_monthly_out_file_path = sl_actual_monthly_out_file_path(month_id) 
            RasterSave(wind_erosion_m, SL_monthly_out_file_path, wf_file_path)
            
            SL_sum +=  wind_erosion_m
            SL_wo_veg_sum += wind_erosion_pot_m
            SL_fvc_100 += wind_erosion_fvc_100
            SL_fvc_50 += wind_erosion_fvc_50
            SL_fvc_10p += wind_erosion_fvc_10p
            SL_fvc_20p += wind_erosion_fvc_20p
        
            fvc_file_path = frac_veg_cov_file_path(month_id_str)
            fvc_arr = preprocess_veg_coverage(fvc_file_path, 1.0)
       
      
        if scale_factor == 0.5:
            SL_sum_0p5 = SL_sum
            SL_sum_0p5_Exists = True
        if scale_factor == 1.0:
            SL_sum_1p0 = SL_sum
        if scale_factor == 1.5:
            SL_sum_1p5 = SL_sum
            SL_sum_1p5_Exists = True
            
        SL_avg = np.nansum(SL_sum) 
        SL_wo_veg_avg = np.nansum(SL_wo_veg_sum)
        SL_fvc_100_avg = np.nansum(SL_fvc_100)
        SL_fvc_50_avg = np.nansum(SL_fvc_50)
        SL_fvc_10p_avg = np.nansum(SL_fvc_10p)
        SL_fvc_20p_avg = np.nansum(SL_fvc_20p)
     
            
        SL_avg_arr[fact_ind] = SL_avg
        SL_wo_veg_avg_arr[fact_ind] = SL_wo_veg_avg
        SL_fvc_100_avg_arr[fact_ind] = SL_fvc_100_avg
        SL_fvc_50_avg_arr[fact_ind] = SL_fvc_50_avg
        SL_fvc_10p_avg_arr[fact_ind] = SL_fvc_10p_avg
        SL_fvc_20p_avg_arr[fact_ind] = SL_fvc_20p_avg
        
        fact_ind+=1
        SL_out_file_path = sl_actual_out_file_path() 
        RasterSave(SL_sum, SL_out_file_path, wf_file_path)


