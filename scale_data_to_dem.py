#!/usr/bin/env python
# coding: utf-8


from osgeo import gdal, gdalconst
import numpy as np
import numpy as np
import gdal
import os

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
    
def scaleToDEM(inputfile_path , outputfilepath):
    inputfile = inputfile_path
    inp = gdal.Open(inputfile, gdalconst.GA_ReadOnly)
    inputProj = inp.GetProjection()
    inputTrans = inp.GetGeoTransform()

    referencefile = 'data/dem_clipped.tif'
    reference = gdal.Open(referencefile, gdalconst.GA_ReadOnly)
    referenceProj = reference.GetProjection()
    referenceTrans = reference.GetGeoTransform()
    bandreference = reference.GetRasterBand(1)    
    x = reference.RasterXSize 
    y = reference.RasterYSize


    outputfile = outputfilepath
    driver= gdal.GetDriverByName('GTiff')
    output = driver.Create(outputfile,x,y,1,bandreference.DataType)
    output.SetGeoTransform(referenceTrans)
    output.SetProjection(referenceProj)

    gdal.ReprojectImage(inp,output,inputProj,referenceProj,gdalconst.GRA_Bilinear)
    
    del output
    
#####Resize all files to DEM file resolution and size #####
data_dir  = 'data/'
out_dir = 'data_scaled/'
utm_dir = 'data_UTM/input/'

dem_file_path = 'dem_clipped.tif'
dem_in_file = os.path.join(data_dir, dem_file_path)
dem_out_path = os.path.join(utm_dir, dem_file_path)
gdal.Warp(dem_out_path, dem_in_file, srcSRS='EPSG:4326', dstSRS='EPSG:32648')

##Average Number of Rain Days per Month
rain_days = [0,0,0,0,0,1,3,3,1,0,0,0]

for k in range(0, 12):
    num = str(k+1)
    if k+1 < 10:
        num = '0'+str(k+1)

    temp_file_path = 'temperature_degC/' + 'wc2.0_30s_tave_'+ num + '.tif'
    temp_in_path = os.path.join(data_dir, temp_file_path)
    scaled_temp = os.path.join(out_dir, temp_file_path)
    scaleToDEM(temp_in_path, scaled_temp)
    temp_out_path = os.path.join(utm_dir, temp_file_path)
    gdal.Warp(temp_out_path, scaled_temp, srcSRS='EPSG:4326', dstSRS='EPSG:32648')
    
    prcp_file_path = 'precipitation_mm/' + 'wc2.0_30s_prec_' + num + '.tif'
    prcp_in_path = os.path.join(data_dir, prcp_file_path)
    scaled_prcp = os.path.join(out_dir, prcp_file_path)
    scaleToDEM(prcp_in_path, scaled_prcp)
    prcp_out_path = os.path.join(utm_dir, prcp_file_path)
    gdal.Warp(prcp_out_path, scaled_prcp, srcSRS='EPSG:4326', dstSRS='EPSG:32648')
    
    sol_file_path = 'solar_radiation/'+ 'shortwave_radiation_' + str(k + 1) + '.tif'
    sol_in_path = os.path.join(data_dir, sol_file_path)
    scaled_sol = os.path.join(out_dir, sol_file_path)
    scaleToDEM(sol_in_path, scaled_sol)
    sol_out_path = os.path.join(utm_dir, sol_file_path)
    gdal.Warp(sol_out_path, scaled_sol, srcSRS='EPSG:4326', dstSRS='EPSG:32648')
    
   #prcp_days_file_path = 'month_prcp_day/prcp_day_' + str(k + 1) + '.tif'
   #prcp_days_file_path = os.path.join(data_dir, prcp_days_file_path)
   #precip_days = get_monthly_num_rain_days(prcp_days_file_path)


    ###TEMP PRCP DAYS
    raster = gdal.Open('data_UTM/input/dem_clipped.tif')
    array = raster.ReadAsArray(0, 0, raster.RasterXSize, raster.RasterYSize)
    prcp_temp = (array < 0) * 0 + (array >= 0) * raindays[k]
    prcp_temp = np.tile(raindays[k], (raster.RasterYSize, raster.RasterXSize))
    prcp_days_file_path = 'month_prcp_day/prcp_day_' + str(k + 1) + '.tif'
    prcp_days_file_path = os.path.join(utm_dir, prcp_days_file_path)
    RasterSave(prcp_temp, prcp_days_file_path, raster)
    
    snow_factor_file_path = 'snow_gm2/' + 'snow_' + str(k + 1) + '.tif'
    snow_factor_in_path = os.path.join(data_dir, snow_factor_file_path)
    scaled_snow = os.path.join(out_dir, snow_factor_file_path)
    scaleToDEM(snow_factor_in_path, scaled_snow)
    snow_factor_out_path = os.path.join(utm_dir, snow_factor_file_path)
    gdal.Warp(snow_factor_out_path, scaled_snow, srcSRS='EPSG:4326', dstSRS='EPSG:32648')
    
    fvc_file_path = 'vegetation_percent_cover/' + 'vegetation_percent_cover_' + str(k+1) + '.tif'
    fvc_in_path = os.path.join(data_dir, fvc_file_path)
    scaled_fvc = os.path.join(out_dir, fvc_file_path)
    scaleToDEM(fvc_in_path, scaled_fvc)
    fvc_out_path = os.path.join(utm_dir, fvc_file_path)
    gdal.Warp(fvc_out_path, scaled_fvc, srcSRS='EPSG:4326', dstSRS='EPSG:32648')

    wind_spd_file_path = 'wind_speed_monthly_clipped/' + 'wind_speed_' + num + '.tif'
    wind_spd_in_path = os.path.join(data_dir, wind_spd_file_path)
    scaled_wind_spd = os.path.join(out_dir, wind_spd_file_path)
    scaleToDEM(wind_spd_in_path, scaled_wind_spd)
    wind_spd_out_path = os.path.join(utm_dir, wind_spd_file_path)
    gdal.Warp(wind_spd_out_path, scaled_wind_spd, srcSRS='EPSG:4326', dstSRS='EPSG:32648')

sand_file_path = 'soil/sand.tif'
sand_in_path = os.path.join(data_dir, sand_file_path)
scaled_sand = os.path.join(out_dir, sand_file_path)
scaleToDEM(sand_in_path, scaled_sand)
sand_out_path = os.path.join(utm_dir, sand_file_path)
gdal.Warp(sand_out_path, scaled_sand, srcSRS='EPSG:4326', dstSRS='EPSG:32648')

silt_file_path = 'soil/silt.tif'
silt_in_path = os.path.join(data_dir, silt_file_path)
scaled_silt = os.path.join(out_dir, silt_file_path)
scaleToDEM(silt_in_path, scaled_silt)
silt_out_path = os.path.join(utm_dir, silt_file_path)
gdal.Warp(silt_out_path, scaled_silt, srcSRS='EPSG:4326', dstSRS='EPSG:32648')


om_file_path = 'soil/soil_organic_matter_gm2.tif'
om_in_path = os.path.join(data_dir, om_file_path)
scaled_om = os.path.join(out_dir, om_file_path)
scaleToDEM(silt_in_path, scaled_om)
om_out_path = os.path.join(utm_dir, om_file_path)
gdal.Warp(om_out_path, scaled_om, srcSRS='EPSG:4326', dstSRS='EPSG:32648')


clay_file_path = 'soil/clay.tif'
clay_in_path = os.path.join(data_dir, clay_file_path)
scaled_clay = os.path.join(out_dir, clay_file_path)
scaleToDEM(clay_in_path, scaled_clay)
clay_out_path = os.path.join(utm_dir, clay_file_path)
gdal.Warp(clay_out_path, scaled_clay, srcSRS='EPSG:4326', dstSRS='EPSG:32648')