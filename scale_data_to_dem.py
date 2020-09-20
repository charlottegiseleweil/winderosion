#!/usr/bin/env python
# coding: utf-8

# In[53]:


from osgeo import gdal, gdalconst
import pygeoprocessing as pygp
import numpy as np
import numpy as np
import gdal
import os

def make_directory(path):
    if not os.path.exists(path):
        os.makedirs(path)
        
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

def scaleToRef(inputfile_path , outputfile_path, referencefile_path):
    
    reference = gdal.Open(referencefile_path, gdalconst.GA_ReadOnly)
    referenceTrans = reference.GetGeoTransform()
    ref_pixel_x = referenceTrans[1]
    ref_pixel_y = referenceTrans[5]

    pygp.warp_raster(inputfile_path, [ref_pixel_x,ref_pixel_y],outputfile_path,'bilinear', None, None)  
    
def scaleToRef_and_UTM(inputfile_path , scaledfile_path, utmfile_path, ref_for_scale_file_path, ref_for_utm_file_path):
    
    reference = gdal.Open(ref_for_scale_file_path, gdalconst.GA_ReadOnly)
    referenceTrans = reference.GetGeoTransform()
    ref_pixel_x = referenceTrans[1]
    ref_pixel_y = referenceTrans[5]
    target_bb = pygp.get_raster_info(ref_for_scale_file_path)['bounding_box']
    
    #rescaling to match resoultion of reference
    pygp.warp_raster(inputfile_path, [ref_pixel_x,ref_pixel_y],scaledfile_path,'bilinear', target_bb, None)  

    #reproject to UTM
    gdal.Warp(utmfile_path, scaledfile_path, srcSRS='EPSG:4326', dstSRS='EPSG:32648')

def fvc_ndv(fvc_path, dem_path):
    fvc_arr = read_raster_as_array(fvc_path)
    dem_arr = read_raster_as_array(dem_path)
    fvc_arr[np.isnan(dem_arr)] = -9999                                
    fvc_arr[np.isnan(fvc_arr)] = 27.4                          
    fvc_arr[fvc_arr == -9999] = np.nan
    RasterSave(fvc_arr, fvc_path, dem_path)
    
####################
#### MAIN CODE ####
###################


#####Resize all files to DEM file resolution and size #####
data_dir  = 'data_aoi/'
scaled_dir = 'data_scaled/'
make_directory(scaled_dir)
make_directory(scaled_dir+'temperature_degC/')
make_directory(scaled_dir+'precipitation_mm/')
make_directory(scaled_dir+'solar_radiation/')
#make_directory(scaled_dir+'month_prcp_day/')
make_directory(scaled_dir+'snow_gm2/')
make_directory(scaled_dir+'vegetation_percent_cover/')
make_directory(scaled_dir+'wind_speed_monthly_clipped/')
make_directory(scaled_dir+'soil/')

utm_dir = 'data_UTM/input/'
make_directory(utm_dir)
make_directory(utm_dir+'temperature_degC/')
make_directory(utm_dir+'precipitation_mm/')
make_directory(utm_dir+'solar_radiation/')
make_directory(utm_dir+'month_prcp_day/')
make_directory(utm_dir+'snow_gm2/')
make_directory(utm_dir+'vegetation_percent_cover/')
make_directory(utm_dir+'wind_speed_monthly_clipped/')
make_directory(utm_dir+'soil/')

dem_file_name = 'dem_aoi.tif'
dem_file_path = os.path.join(data_dir, dem_file_name)
dem_utm_path = os.path.join(utm_dir, dem_file_name)
gdal.Warp(dem_utm_path, dem_file_path, srcSRS='EPSG:4326', dstSRS='EPSG:32648')

dem_arr = read_raster_as_array(dem_utm_path)
dem_arr[dem_arr == 0] = np.nan
RasterSave(dem_arr, dem_utm_path, dem_utm_path)

raindays = [0.0,0.0,0.0,0.0,0.0,1.0,3.0,3.0,1.0,0.0,0.0,0.0]


for k in range(0, 12):
    num = str(k+1)
    if k+1 < 10:
        num = '0'+str(k+1)

    temp_file_path = 'temperature_degC/' + 'wc2.0_30s_tave_'+ num + '.tif'
    temp_in_path = os.path.join(data_dir, temp_file_path)
    scaled_temp = os.path.join(scaled_dir, temp_file_path)
    temp_out_path = os.path.join(utm_dir, temp_file_path)
    scaleToRef_and_UTM(temp_in_path, scaled_temp, temp_out_path, dem_file_path, dem_utm_path)

    prcp_file_path = 'precipitation_mm/' + 'wc2.0_30s_prec_' + num + '.tif'
    prcp_in_path = os.path.join(data_dir, prcp_file_path)
    scaled_prcp = os.path.join(scaled_dir, prcp_file_path)
    prcp_out_path = os.path.join(utm_dir, prcp_file_path)
    scaleToRef_and_UTM(prcp_in_path, scaled_prcp, prcp_out_path, dem_file_path, dem_utm_path)
    
    sol_file_path = 'solar_radiation/'+ 'shortwave_radiation_' + str(k + 1) + '.tif'
    sol_in_path = os.path.join(data_dir, sol_file_path)
    scaled_sol = os.path.join(scaled_dir, sol_file_path)
    sol_out_path = os.path.join(utm_dir, sol_file_path)
    scaleToRef_and_UTM(sol_in_path, scaled_sol, sol_out_path,dem_file_path, dem_utm_path)

    raster = gdal.Open(dem_utm_path)
    array = raster.ReadAsArray(0, 0, raster.RasterXSize, raster.RasterYSize)
    prcp_temp = (array < 0) * 0 + (array >= 0) * raindays[k]
    prcp_temp = np.tile(raindays[k], (raster.RasterYSize, raster.RasterXSize))
    prcp_days_file_path = 'month_prcp_day/prcp_day_' + str(k + 1) + '.tif'
    prcp_days_file_path = os.path.join(utm_dir, prcp_days_file_path)
    RasterSave(prcp_temp, prcp_days_file_path, dem_utm_path)
    
    snow_factor_file_path = 'snow_gm2/' + 'snow_' + str(k + 1) + '.tif'
    snow_factor_in_path = os.path.join(data_dir, snow_factor_file_path)
    scaled_snow = os.path.join(scaled_dir, snow_factor_file_path)
    snow_factor_out_path = os.path.join(utm_dir, snow_factor_file_path)
    scaleToRef_and_UTM(snow_factor_in_path, scaled_snow, snow_factor_out_path, dem_file_path, dem_utm_path)

    fvc_file_path = 'vegetation_percent_cover/' + 'vegetation_percent_cover_' + str(k+1) + '.tif'
    fvc_in_path = os.path.join(data_dir, fvc_file_path)
    scaled_fvc = os.path.join(scaled_dir, fvc_file_path)
    fvc_out_path = os.path.join(utm_dir, fvc_file_path)
    scaleToRef_and_UTM(fvc_in_path, scaled_fvc, fvc_out_path, dem_file_path, dem_utm_path)
    fvc_ndv(fvc_out_path, dem_utm_path)
    
    fvc_arr = read_raster_as_array(fvc_out_path)
    
    if(k == 0):
        fvc_sum = fvc_arr
    else:
        fvc_sum = fvc_sum+fvc_arr
    
    wind_spd_file_path = 'wind_speed_monthly_clipped/' + 'wind_speed_' + num + '.tif'
    wind_spd_in_path = os.path.join(data_dir, wind_spd_file_path)
    scaled_wind_spd = os.path.join(scaled_dir, wind_spd_file_path)
    wind_spd_out_path = os.path.join(utm_dir, wind_spd_file_path)
    scaleToRef_and_UTM(wind_spd_in_path, scaled_wind_spd, wind_spd_out_path,dem_file_path, dem_utm_path)

fvc_avg= fvc_sum/12.0
fvc_avg_file_path = 'vegetation_percent_cover/' + 'vegetation_percent_cover_avg.tif'
fvc_avg_out_file_path = os.path.join(utm_dir, fvc_avg_file_path)
RasterSave(fvc_avg, fvc_avg_out_file_path, dem_utm_path)
fvc_ndv(fvc_avg_out_file_path, dem_utm_path)

fvc_avg_arr = read_raster_as_array(fvc_avg_out_file_path)
fvc_100 = fvc_avg_arr
fvc_100[fvc_100 != np.nan] = 100.0
fvc_100_file_path = 'vegetation_percent_cover/' + 'vegetation_percent_cover_100.tif'
fvc_100_out_file_path = os.path.join(utm_dir, fvc_100_file_path)
RasterSave(fvc_100, fvc_100_out_file_path, dem_utm_path)
fvc_ndv(fvc_100_out_file_path, dem_utm_path)

fvc_50 = fvc_avg_arr
fvc_50[fvc_50 != np.nan] = 50.0
fvc_50_file_path = 'vegetation_percent_cover/' + 'vegetation_percent_cover_50.tif'
fvc_50_out_file_path = os.path.join(utm_dir, fvc_50_file_path)
RasterSave(fvc_50, fvc_50_out_file_path, dem_utm_path)
fvc_ndv(fvc_50_out_file_path, dem_utm_path)

sand_file_path = 'soil/sand.tif'
sand_in_path = os.path.join(data_dir, sand_file_path)
scaled_sand = os.path.join(scaled_dir, sand_file_path)
sand_out_path = os.path.join(utm_dir, sand_file_path)
scaleToRef_and_UTM(sand_in_path, scaled_sand, sand_out_path,dem_file_path, dem_utm_path)

silt_file_path = 'soil/silt.tif'
silt_in_path = os.path.join(data_dir, silt_file_path)
scaled_silt = os.path.join(scaled_dir, silt_file_path)
silt_out_path = os.path.join(utm_dir, silt_file_path)
scaleToRef_and_UTM(silt_in_path, scaled_silt, silt_out_path, dem_file_path, dem_utm_path)

om_file_path = 'soil/soil_organic_matter_gm2.tif'
om_in_path = os.path.join(data_dir, om_file_path)
scaled_om = os.path.join(scaled_dir, om_file_path)
om_out_path = os.path.join(utm_dir, om_file_path)
scaleToRef_and_UTM(om_in_path, scaled_om, om_out_path,dem_file_path, dem_utm_path)

clay_file_path = 'soil/clay.tif'
clay_in_path = os.path.join(data_dir, clay_file_path)
scaled_clay = os.path.join(scaled_dir, clay_file_path)
clay_out_path = os.path.join(utm_dir, clay_file_path)
scaleToRef_and_UTM(clay_in_path, scaled_clay, clay_out_path,dem_file_path, dem_utm_path)

