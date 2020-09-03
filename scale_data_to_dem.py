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
    
def scaleToRef(inputfile_path , outputfile_path, referencefile_path):
    
    reference = gdal.Open(referencefile_path, gdalconst.GA_ReadOnly)
    referenceTrans = reference.GetGeoTransform()
    ref_pixel_x = referenceTrans[1]
    ref_pixel_y = referenceTrans[5]

    pygp.warp_raster(inputfile_path, [ref_pixel_x,ref_pixel_y],outputfile_path,'bilinear', None, None)  
    
def scaleToRef_and_UTM(inputfile_path , scaledfile_path, utmfile_path, referencefile_path):
    
    reference = gdal.Open(referencefile_path, gdalconst.GA_ReadOnly)
    referenceTrans = reference.GetGeoTransform()
    ref_pixel_x = referenceTrans[1]
    ref_pixel_y = referenceTrans[5]
    
    #rescaling to match resoultion of reference
    pygp.warp_raster(inputfile_path, [ref_pixel_x,ref_pixel_y],scaledfile_path,'bilinear', None, None)  
    
    #reproject to UTM
    gdal.Warp(utmfile_path, scaledfile_path, srcSRS='EPSG:4326', dstSRS='EPSG:32648')

    #force reprojected raster to same dimensions as reference
    base_raster_path_list = [referencefile_path, utmfile_path]
    target_raster_path_list = [referencefile_path, utmfile_path]
    resample_method_list = ['bilinear', 'bilinear']
    target_pixel_size = [ref_pixel_x, ref_pixel_y]
    bounding_box_mode= 'intersection'
    base_vector_path_list=None 
    raster_align_index = 0
    pygp.align_and_resize_raster_stack(
            base_raster_path_list, target_raster_path_list, resample_method_list,
            target_pixel_size, bounding_box_mode, base_vector_path_list,
            raster_align_index)
    


####################
#### MAIN CODE ####
###################


#####Resize all files to DEM file resolution and size #####
data_dir  = 'data_aoi/'
scaled_dir = 'data_scaled/'
utm_dir = 'data_UTM/input/'

dem_file_name = 'dem_aoi.tif'
dem_file_path = os.path.join(data_dir, dem_file_name)
dem_out_path = os.path.join(utm_dir, dem_file_name)
gdal.Warp(dem_out_path, dem_in_file, srcSRS='EPSG:4326', dstSRS='EPSG:32648')

raindays = [0,0,0,0,0,1,3,3,1,0,0,0]

for k in range(0, 12):
    num = str(k+1)
    if k+1 < 10:
        num = '0'+str(k+1)

    temp_file_path = 'temperature_degC/' + 'wc2.0_30s_tave_'+ num + '.tif'
    temp_in_path = os.path.join(data_dir, temp_file_path)
    scaled_temp = os.path.join(scaled_dir, temp_file_path)
    temp_out_path = os.path.join(utm_dir, temp_file_path)
    scaleToRef_and_UTM(temp_in_path, scaled_temp, temp_out_path, dem_out_path)
    
    prcp_file_path = 'precipitation_mm/' + 'wc2.0_30s_prec_' + num + '.tif'
    prcp_in_path = os.path.join(data_dir, prcp_file_path)
    scaled_prcp = os.path.join(scaled_dir, prcp_file_path)
    prcp_out_path = os.path.join(utm_dir, prcp_file_path)
    scaleToRef_and_UTM(prcp_in_path, scaled_prcp, prcp_out_path, dem_out_path)
    
    sol_file_path = 'solar_radiation/'+ 'shortwave_radiation_' + str(k + 1) + '.tif'
    sol_in_path = os.path.join(data_dir, sol_file_path)
    scaled_sol = os.path.join(scaled_dir, sol_file_path)
    sol_out_path = os.path.join(utm_dir, sol_file_path)
    scaleToRef_and_UTM(sol_in_path, scaled_sol, sol_out_path, dem_out_path)

    raster = gdal.Open(dem_out_path)
    array = raster.ReadAsArray(0, 0, raster.RasterXSize, raster.RasterYSize)
    prcp_temp = (array < 0) * 0 + (array >= 0) * raindays[k]
    prcp_temp = np.tile(raindays[k], (raster.RasterYSize, raster.RasterXSize))
    prcp_days_file_path = 'month_prcp_day/prcp_day_' + str(k + 1) + '.tif'
    prcp_days_file_path = os.path.join(utm_dir, prcp_days_file_path)
    RasterSave(prcp_temp, prcp_days_file_path, raster)
    
    snow_factor_file_path = 'snow_gm2/' + 'snow_' + str(k + 1) + '.tif'
    snow_factor_in_path = os.path.join(data_dir, snow_factor_file_path)
    scaled_snow = os.path.join(scaled_dir, snow_factor_file_path)
    snow_factor_out_path = os.path.join(utm_dir, snow_factor_file_path)
    scaleToRef_and_UTM(snow_factor_in_path, scaled_snow, snow_factor_out_path, dem_out_path)

    fvc_file_path = 'vegetation_percent_cover/' + 'vegetation_percent_cover_' + str(k+1) + '.tif'
    fvc_in_path = os.path.join(data_dir, fvc_file_path)
    scaled_fvc = os.path.join(scaled_dir, fvc_file_path)
    fvc_out_path = os.path.join(utm_dir, fvc_file_path)
    scaleToRef_and_UTM(fvc_in_path, scaled_fvc, fvc_out_path, dem_out_path)

    wind_spd_file_path = 'wind_speed_monthly_clipped/' + 'wind_speed_' + num + '.tif'
    wind_spd_in_path = os.path.join(data_dir, wind_spd_file_path)
    scaled_wind_spd = os.path.join(scaled_dir, wind_spd_file_path)
    wind_spd_out_path = os.path.join(utm_dir, wind_spd_file_path)
    scaleToRef_and_UTM(wind_spd_in_path, scaled_wind_spd, wind_spd_out_path, dem_out_path)

sand_file_path = 'soil/sand.tif'
sand_in_path = os.path.join(data_dir, sand_file_path)
scaled_sand = os.path.join(scaled_dir, sand_file_path)
sand_out_path = os.path.join(utm_dir, sand_file_path)
scaleToRef_and_UTM(sand_in_path, scaled_sand, sand_out_path, dem_out_path)

silt_file_path = 'soil/silt.tif'
silt_in_path = os.path.join(data_dir, silt_file_path)
scaled_silt = os.path.join(scaled_dir, silt_file_path)
silt_out_path = os.path.join(utm_dir, silt_file_path)
scaleToRef_and_UTM(silt_in_path, scaled_silt, silt_out_path, dem_out_path)

om_file_path = 'soil/soil_organic_matter_gm2.tif'
om_in_path = os.path.join(data_dir, om_file_path)
scaled_om = os.path.join(scaled_dir, om_file_path)
om_out_path = os.path.join(utm_dir, om_file_path)
scaleToRef_and_UTM(om_in_path, scaled_om, om_out_path, dem_out_path)

clay_file_path = 'soil/clay.tif'
clay_in_path = os.path.join(data_dir, clay_file_path)
scaled_clay = os.path.join(scaled_dir, clay_file_path)
clay_out_path = os.path.join(utm_dir, clay_file_path)
scaleToRef_and_UTM(clay_in_path, scaled_clay, clay_out_path, dem_out_path)

