#!/usr/bin/env python
# coding: utf-8

# In[1]:


## Pseudo Code for Wind Erosion Model
## Written by Charlie Weil & Isita Talukdar, August 2020
## Inspired by Huang Binbin's approach.

import numpy as np
import gdal
import os

def execute(args):
    ## See typical InVEST Model structure.
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

def get_elevation_data(file_name):
    """Get Elevation Data from Raster File

    Parameters:
        file_name: name of Raster file

    Returns:
        dem: Elevation data, in meters
    """

    dem = gdal.Open(file_name)
    return dem

def calculate_air_pressure(dem):
    """Calculate air pressure

    Parameters:
        dem: Elevation data, in meters

    Returns:
        Air Pressure (P), in UNITS
    """
    dem_a = dem.ReadAsArray(0, 0, dem.RasterXSize, dem.RasterYSize)
    dem_a2 = (dem_a < 0) * 0 + (dem_a >= 0) * dem_a
    air_pressure = 101.3 * (1 - 0.0255 * dem_a2 / 1000 * (6357 / (6357 + dem_a2 / 1000))) ** 5.256
    return air_pressure

def get_monthly_avg_temp(file_name):
    """Get Average Temperature for a given month from Raster File

    Parameters:
        file_name: name of Raster file

    Returns:
        tem_m_a2: Temperature for a given month, in °C
    """
    
    # monthly average temperature(℃)
    tem_m = gdal.Open(file_name)
    tem_m_a = tem_m.ReadAsArray(0, 0, tem_m.RasterXSize, tem_m.RasterYSize)
    # multiply 0.1 because the unit of source data is 0.1 ℃
    tem_m_a2 = (tem_m_a < -100) * 0.1 + (tem_m_a >= -100) * tem_m_a
    return tem_m_a2

def get_monthly_total_precip(file_name):
    """Get Total Precipitation for a given month from Raster File

    Parameters:
        file_name: name of Raster file

    Returns:
        prcp_a: Total Precipitation for a given month, in mm
    """
    # monthly total precipitation(mm)
    prcp = gdal.Open(file_name)
    prcp_a = prcp.ReadAsArray(0, 0, prcp.RasterXSize, prcp.RasterYSize)
    # multiply 0.1 because the unit of source data is 0.1 mm
    prcp_a2 = (prcp_a < 0) * 0 + (prcp_a >= 0) * prcp_a
    return prcp_a

def get_monthly_sol_rad(file_name):
    """Get Total Solar Radiation for a given month from Raster File

    Parameters:
        file_name: name of Raster file

    Returns:
        SOL_a2: Total Solar Radiation for a given month, in MJ/m2
    """

    # monthly total solar radiation (MJ/m2)
    SOL = gdal.Open(file_name)
    SOL_a = SOL.ReadAsArray(0, 0, SOL.RasterXSize, SOL.RasterYSize)
    SOL_a2 = (SOL_a < 0) * 0.1 + (SOL_a >= 0) * SOL_a
    return SOL_a2

def get_monthly_num_rain_days(file_name):
    """Get Number of Rain Days for a given month from Raster File

    Parameters:
        file_name: name of Raster file

    Returns:
        prcp_days_a2: Total Number of Rain Days for a given month, in days
    """

    # days of rain events in every month
    prcp_days = gdal.Open(file_name)
    prcp_days_a = prcp_days.ReadAsArray(0, 0, prcp_days.RasterXSize, prcp_days.RasterYSize)
    prcp_days_a2 = (prcp_days_a < 0) * 0 + (prcp_days_a >= 0) * prcp_days_a
    return prcp_days_a2

def calculate_evapotranspiration(solar_rad, avg_temp):
    """Calculate the potential evapotranspiration for a given month

    Parameters:
        Temperature for a given month, in °C
        Total Solar Radiation for a given month, in Megajoules/meters squared

    Returns:
        Potential Evapotranspiration (ETp), in millimeters
    """
    evap_trans = 0.0135 * (solar_rad / 2.54) * (avg_temp + 17.8)
    
    return evap_trans

def calculate_soil_moisture(solar_rad, avg_temp, precip, days_precip):
    """Calculate the Soil Moisture Factor for a given month

    Parameters:
        Potential Evapotranspiration for a given month, in millimeters
        Number of Rain days for a given month

    Returns:
        Soil Moisture Factor (SW)
    """
    evap_trans = calculate_evapotranspiration(solar_rad, avg_temp)
    soil_moisture = (evap_trans - precip * days_precip) / evap_trans
    soil_moisture = (soil_moisture < 0) * 0 + (soil_moisture >= 0) * soil_moisture
    
    return soil_moisture

def calculate_snow_factor(file_name):
    """Calculate the Snow Factor for a given month

    Parameters:
        Snow Cover, in UNITS

    Returns:
        Snow Factor (SD)
    """
    SD = gdal.Open(file_name)
    SD_a = SD.ReadAsArray(0, 0, SD.RasterXSize, SD.RasterYSize)
    SD_a2 = (SD_a < 0) * 0 + (SD_a >= 0) * SD_a
    SD_a2 = (1 - SD_a2 * 0.01)
    return SD_a2

def calculate_air_density(temperature, air_pressure):
    """Calculate air density for a given month.

    Parameters:
        temperature: Temperature for a given month, in °C
        air_pressure: Air Pressure, in UNITS

    Returns:
        Air density (rho), in kg/m^3
    """
    air_density_rho = 1.293 * (273 / (273 + temperature)) * air_pressure / 101.3

    return air_density_rho

def calculate_wind_factor_daily(file_name):
    """Calculate the wind factor for a given day

    Parameters:
        wind speed for a day, in m/s

    Returns:
        Wind Factor (wf)
    """
     # daily wind speed (m/s)
    wind = gdal.Open(file_name)
    wind_a = wind.ReadAsArray(0, 0, wind.RasterXSize, wind.RasterYSize)

    # multiply 0.1 because the unit of source data is 0.1 m/s
    wind_a2 = (wind_a < 50) * 0 + ((wind_a >= 50) & (wind_a < 2000)) * wind_a * 0.1 + (wind_a >= 2000) * 200

    # wind factor
    wf = wind_a2 * (wind_a2 - 5) ** 2
    
    return wf

def calculate_daily_weather_factor(wind_factor_d, air_density_rho, soil_moisture, snow_factor):
    """Calculate the Weather Factor for a given day 

    Parameters:
        Wind Factor (wf)
        Air density (rho), in kg/m^3
        Soil Moisture Factor (SW)
        Snow Factor (SD)
        
    Returns:
        Weather Factor for a given day
    """
    weather_factor_d = wind_factor_d * air_density_rho / 9.8 * soil_moisture * snow_factor
    
    #print(day_index, " : ",np.sum(snow_factor))
    return weather_factor_d

def calculate_monthly_weather_factor(month_index, temp, precip, sol_rad, precip_days, snow_factor, pressure):
    """Calculate the Weather Factor for a given month

    Parameters:
        Wind Factor (wf)
        Air density (rho), in kg/m^3
        Soil Moisture Factor (SW)
        Snow Factor (SD)
        
    Returns:
        Weather Factor for a given month (WF)
    """
   
    sw = calculate_soil_moisture(sol_rad, temp, precip, precip_days)

    air_density_rho = calculate_air_density(temp, pressure)
    
   
    weather_factor_m = 0.0
    #start_date = (month_index == 0) * 0 + (month_index > 0) * int(endday[month_index - 1])
    start_date = 0 if (month_index == 0) else int(endday[month_index - 1])
    end_date = int(endday[month_index]) + 1
    #print("Start Date: ", start_date)
    #print(" End Date: ", end_date)
    #print('\n')
    #Loop over every day
    for i in range(start_date, end_date):
        wind_file_name = 'month_wind_day/wind_' + str(i + 1) + '.tif'
        wind_factor_d =  calculate_wind_factor_daily(wind_file_name)
        weather_factor_d = calculate_daily_weather_factor(wind_factor_d, air_density_rho, sw, snow_factor)
        weather_factor_m = weather_factor_m + weather_factor_d
        
    return weather_factor_m

def get_sand_ratio(file_name):
    """Get Sand Ratio from Raster File

    Parameters:
        file_name: name of Raster file

    Returns:
        sand_a2: Sand Ratio
    """

    #ratio of sand(%)
    sand = gdal.Open(file_name)
    sand_a = sand.ReadAsArray(0, 0, sand.RasterXSize, sand.RasterYSize)
    sand_a2 = (sand_a <= 0) * 0.1 + (sand_a >0) * sand_a
    return sand_a2

def get_silt_ratio(file_name):
    """Get Silt Ratio from Raster File

    Parameters:
        file_name: name of Raster file

    Returns:
        silt_a2: Silt Ratio
    """

    #ratio of silt(%)
    silt = gdal.Open(file_name)
    silt_a = silt.ReadAsArray(0, 0, silt.RasterXSize, silt.RasterYSize)
    silt_a2 = (silt_a <= 0) * 0.1 + (silt_a > 0) * silt_a
    return silt_a2

def get_org_mat_ratio(file_name):
    """Get Organic Matter Ratio from Raster File

    Parameters:
        file_name: name of Raster file

    Returns:
        om_a2: Organic Matter Ratio
    """
    #ratio of organic matter(%)
    om = gdal.Open(file_name)
    om_a = om.ReadAsArray(0, 0, om.RasterXSize, om.RasterYSize)
    om_a2 = (om_a <= 0) * 0.1 + (om_a > 0) * om_a
    return om_a2

def get_clay_ratio(file_name):
    """Get Clay Ratio from Raster File

    Parameters:
        file_name: name of Raster file

    Returns:
        clay_a2: Clay Ratio
    """
    #ratio of clay(%)
    clay_name = 'soil/clay.tif'
    clay = gdal.Open(file_name)
    clay_a = clay.ReadAsArray(0, 0, clay.RasterXSize, clay.RasterYSize)
    clay_a2 = (clay_a <= 0) * 0.1 + (clay_a >0) * clay_a
    return clay_a2

def calculate_soil_crust_factor(clay_ratio, org_mat_ratio):
    """Calculate the Soil Crusting Factor

    Parameters:
        Clay Ratio(%)
        Organic Matter Ratio(%)

    Returns:
        Soil Crusting Factor (SCF)
    """
    soil_crust_factor = SCF = 1 / (1 + 0.0066 * clay_ratio ** 2 + 0.021 * org_mat_ratio ** 2)

    return soil_crust_factor

def calculate_soil_erodibility_factor(sand_ratio, silt_ratio, clay_ratio, org_mat_ratio):
    """Calculate the Soil Erodibility Factor

    Parameters:
        Sand Ratio(%)
        Silt Ratio(%)
        Clay Ratio(%)
        Organic Matter Ratio(%)

    Returns:
        Soil Erodibility Factor (EF)
    """
    soil_erode_factor = (29.09 + 0.31 * sand_ratio + 0.17 * silt_ratio + 0.33 * sand_ratio / clay_ratio - 2.59 * org_mat_ratio - 0.95 * 0) / 100


    return soil_erode_factor

def calculate_vegetation_factor(file_name):
    """Calculate the Fraction of Monthly Vegetation Coverage For a Given Month

    Parameters:
        Fraction of Monthly Vegetation Coverage

    Returns:
        Vegetation Factor(COG)
    """
    # read fraction of vegetation coverage data (%)
    fvc = gdal.Open(file_name)
    fvc_a = fvc.ReadAsArray(0, 0, fvc.RasterXSize, fvc.RasterYSize)
    fvc_a2 = ((fvc_a < 0) | (fvc_a > 100)) * 0 + ((fvc_a >= 0) & (fvc_a <= 100)) * fvc_a
    vegetation_factor = np.exp(-0.00483*fvc_a2)
    
    return vegetation_factor

def calculate_surface_terr_rough():
    # surface roughness factor it can calculated by the smith-Caarson equation
    KK_path = 'KK/KK_1.tif'
    KK = gdal.Open(KK_path)
    KK_a = KK.ReadAsArray(0, 0, KK.RasterXSize, KK.RasterYSize)
    KK_a2 = (KK_a <= 0) * 0.1 + (KK_a > 0) * KK_a
    return KK_a2

def read_WF(file_name):
    WF = gdal.Open(file_name)
    WF_a = WF.ReadAsArray(0, 0, WF.RasterXSize, WF.RasterYSize)
    weather_factor = (WF_a <= 0) * 0.1 + (WF_a > 0) * WF_a
    return weather_factor

def read_COG(file_name):
    COG = gdal.Open(file_name)
    COG_a = COG.ReadAsArray(0, 0, COG.RasterXSize, COG.RasterYSize)
    vegetation_factor = (COG_a <= 0) * 0.1 + (COG_a > 0) * COG_a
    return vegetation_factor


def calculate_monthly_wind_erosion(weather_factor,soil_erode_factor, kprime, soil_crust_factor, cog):
    """Calculate Potential Wind Erosion for a given month

    Parameters:
        Potential RWEQ Maximum Horizontal Flux (Qmax)
        Potential Critical Field Length (s)
        
    Returns:
         Potential Wind Erosion for a given month (SL)
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


########### MAIN CODE ###########
# # # # # # # # # # # # # # # # #

# days of month 1-12
mondays = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

endday = [30, 58, 89, 119, 150, 180, 211, 242, 272, 303, 333, 364]

m_id = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']

# Step 1 : Weather Factor
# # # # # # # # # # # # #

#get elevation data
dem = get_elevation_data("dem.tif")

#calculate air pressure
pressure = calculate_air_pressure(dem)

#Compute Monthly Weather Factor
for k in range(0, 12):
    temp_file_name = 'month_tem/' + str(k + 1) + '月气温_Pro.tif'
    temp = get_monthly_avg_temp(temp_file_name)
    
    prcp_file_name = 'month_prcp/' + str(k + 1) + '月降水量_Pro.tif'
    precip = get_monthly_total_precip(prcp_file_name)
    
    sol_file_name = 'month_SOL/rad2000' + str(k + 1) + '_grd.tif'
    sol_rad = get_monthly_sol_rad(sol_file_name)
    
    prcp_days_file_name = 'month_prcp_day/prcp_day_' + str(k + 1) + '.tif'
    precip_days = get_monthly_num_rain_days(prcp_days_file_name)
    
    snow_factor_file_name = 'month_SD/SD_' + str(k + 1) + '.tif'
    snow_factor = calculate_snow_factor(snow_factor_file_name)
    
    monthly_WF = calculate_monthly_weather_factor(k, temp, precip, sol_rad, precip_days, snow_factor, pressure)
    wf_out_file_name = 'output/WF/wf_' + str(k + 1) + '.tif'
    RasterSave(monthly_WF, wf_out_file_name, dem)

# Step 2 : Soil Crusting Factor and Erodibility Factor
# # # # # # # # # # # # #
sand_file_name = 'soil/sand.tif'
sand_ratio = get_sand_ratio(sand_file_name)

silt_file_name = 'soil/silt.tif'
silt_ratio = get_silt_ratio(silt_file_name)

om_file_name = 'soil/som.tif'
org_mat_ratio = get_org_mat_ratio(om_file_name)

clay_file_name = 'soil/clay.tif'
clay_ratio = get_clay_ratio(clay_file_name)

scf = calculate_soil_crust_factor(clay_ratio,org_mat_ratio)
scf_out_file_name = 'output/SCF.tif'
RasterSave(scf,scf_out_file_name,dem)

ef = calculate_soil_erodibility_factor(sand_ratio, silt_ratio, clay_ratio, org_mat_ratio)
ef_out_file_name = 'output/EF.tif'
RasterSave(ef,ef_out_file_name,dem)


# Step 3 : Vegetation Factor
# # # # # # # # # # # # #
for i in range(0, 12):
    fvc_file_name = 'FVC2015/A2015' + str(i+1) + '_fc.tif'
    monthly_cog = calculate_vegetation_factor(fvc_file_name)
    COG_out_file_name  = 'output/COG/COG_' + str(i + 1) + '.tif'
    fvc = gdal.Open(fvc_file_name)
    RasterSave(monthly_cog, COG_out_file_name, fvc)
    
# Step 4 : Surface Terrain Roughness Factor
# # # # # # # # # # # # #
kprime = calculate_surface_terr_rough()

# Step 5-6 : Actual and Potential Wind Erosion
# # # # # # # # # # # # #
SL_sum = 0.0
SL_p_sum = 0.0


for i in range(0, 12):
    wf_file_name = 'output/WF/WF_' + str(i + 1) + '.tif'
    wf = read_WF(wf_file_name)
    
    cog_file_name ='output/COG/COG_' + str(i + 1) + '.tif'
    cog = read_COG(cog_file_name)
    
    wind_erosion_m, wind_erosion_pot_m = calculate_monthly_wind_erosion(wf,ef,kprime,scf,cog)
    SL_sum +=  wind_erosion_m
    SL_p_sum += wind_erosion_pot_m

WF_format = gdal.Open(wf_file_name)
SL_out_file_name ='output/SL.tif'
RasterSave(SL_sum, SL_out_file_name, WF_format)

SL_p_out_file_name ='output/SL_p.tif'
RasterSave(SL_p_sum, SL_p_out_file_name, WF_format)

sand_re = SL_p_sum - SL_sum
sand_re = (sand_re < 0) * 0 + (sand_re >= 0) * sand_re

sand_re_out_file_name = 'output/sand_re.tif'
RasterSave(sand_re, sand_re_out_file_name, dem)


