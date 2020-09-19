# Soil loss due to Wind Erosion Model

#### Requirements

* Sample data (can be downloaded [here](https://drive.google.com/file/d/1YQNOEnWyTdu1Po_D2klCpXNbAR3hdDgd/view?usp=sharing))
* A python environment with GDAL installed
* Input Data Location
   * Option 1: All files are located in one directory with the fixed directory and filenames as shown in the the usage example option 1 below. This directory should be the --input_data_dir option. All other directory and file name options will be ignored when --input_data_dir is used
   * Option 2: Each input parameter data directory and filename must be specified as shown in usage example option 2 below. All directories and filenames must be specified. For parameters with monthly data, a file_name prefix must be defined. The specific file name will be a concatenation of this prefix and "_<month_index>.tif" where month_index is from 1 to 12
   
#### Usage example: OPTION 1
```python3 ../wind_erosion.py --sl_output_suffix=actual --input_data_dir=../Inputs --data_dir=../test_out1```

  #### Inputs
	 * ../Inputs/temperature_degC/wc2.0_30s_tave_<n>.tif: monthly average temperature in deg C 
	 * ../Inputs/precipitation_mm/wc2.0_30s_prec_<n>.tif: monthly total precipitation in mm 
	 * ../Inputs/solar_radiation/shortwave_radiation_<n>.tif: monthly solar radiation in MJ/m^2 
	 * ../Inputs/month_prcp_day/prcp_day_<n>.tif: monthly number of rain days  
	 * ../Inputs/snow_gm2/snow_<n>.tif: monthly snow cover factor (probability) 
	 * ../Inputs/wind_speed_monthly_clipped/wind_speed_<n>.tif: monthly average wind speed in m/s 
	 * ../Inputs/vegetation_percent_cover/vegetation_percent_cover_<n>.tif: monthly fractional vegetation coverage in % 
	 * ../Inputs/soil/sand.tif : non-temporal sand ratio (filename is a TIF file including extension) 
	 * ../Inputs/soil/silt.tif : non-temporal silt ratio (filename is a TIF file including extension) 
	 * ../Inputs/soil/clay.tif : non-temporal clay ratio (filename is a TIF file including extension) 
	 * ../Inputs/soil/soil_organic_matter_gm2.tif : non-temporal organic matter ratio, data expected to be percentage*100 (filename is a TIF file including extension)
  
  	 where n is is a month index from 1 to 12
  
  #### Outputs 
	  * ../test_out1/Output/SL_actual_<n>.tif: monthly Soil Loss in kg/m^2
	  * ../test_out1/Output/SL_actual.tif: Annual Soil Loss in kg/m^2

  	  where n is is a month index from 1 to 12
  
  #### Intermediate
	  * ../test_out1/Intermediate/model_intermediate/WF/wf_<n>.tif: monthly Weather Factor in kg/m
	  * ../test_out1/Intermediate/model_intermediate/KK/KK_<n>.tif: monthly Surface Terrain Roughness unitless
	  * ../test_out1/Intermediate/model_intermediate/COG/COG_<n>.tif: monthly Vegeation Factor unitless
	  * ../test_out1/Intermediate/model_intermediate/EF.tif: Erodible Fraction unitless
	  * ../test_out1/Intermediate/model_intermediate/SCF.tif: Soil Crusting Factor unitless

	  where n is is a month index from 1 to 12
  
#### Usage example: OPTION 2
 ```python3 ../wind_erosion.py --data_dir=../test_out2 --dem_dir=../Inputs --dem_file_name=dem_aoi.tif --temper_dir=../Inputs/temperature_degC/ --temper_prefix=wc2.0_30s_tave_ --precip_dir=../Inputs/precipitation_mm/ --precip_prefix=wc2.0_30s_prec_ --sol_dir=../Inputs/solar_radiation/ --sol_prefix=shortwave_radiation_ --prcp_days_dir=../Inputs/month_prcp_day/ --prcp_days_prefix=prcp_day_ --snow_dir=../Inputs/snow_gm2/ --snow_prefix=snow_ --wind_speed_dir=../Inputs/wind_speed_monthly_clipped/ --wind_speed_prefix=wind_speed_ --sand_dir=../Inputs/soil/ --sand_file_name=sand.tif --silt_dir=../Inputs/soil/ --silt_file_name=silt.tif --clay_dir=../Inputs/soil/ --clay_file_name=clay.tif --som_dir=../Inputs/soil/ --som_file_name=soil_organic_matter_gm2.tif --fvc_dir=../Inputs/vegetation_percent_cover/ --fvc_prefix=vegetation_percent_cover_```

  #### Inputs
  	 * ../Inputs/temperature_degC/wc2.0_30s_tave_<n>.tif: monthly average temperature in deg C 
	 * ../Inputs/precipitation_mm/wc2.0_30s_prec_<n>.tif: monthly total precipitation in mm 
	 * ../Inputs/solar_radiation/shortwave_radiation_<n>.tif: monthly solar radiation in MJ/m^2 
	 * ../Inputs/month_prcp_day/prcp_day_<n>.tif: monthly number of rain days  
	 * ../Inputs/snow_gm2/snow_<n>.tif: monthly snow cover factor (probability) 
	 * ../Inputs/wind_speed_monthly_clipped/wind_speed_<n>.tif: monthly average wind speed in m/s 
	 * ../Inputs/vegetation_percent_cover/vegetation_percent_cover_<n>.tif: monthly fractional vegetation coverage in % 
	 * ../Inputs/soil/sand.tif : non-temporal sand ratio (filename is a TIF file including extension) 
	 * ../Inputs/soil/silt.tif : non-temporal silt ratio (filename is a TIF file including extension) 
	 * ../Inputs/soil/clay.tif : non-temporal clay ratio (filename is a TIF file including extension) 
	 * ../Inputs/soil/soil_organic_matter_gm2.tif : non-temporal organic matter ratio, data expected to be percentage*100 (filename is a TIF file including extension)
 
	 where n is is a month index from 1 to 12

  #### Outputs 
	  * ../test_out2/Output/SL_actual_<n>.tif: monthly Soil Loss in kg/m^2
	  * ../test_out2/Output/SL_actual.tif: Annual Soil Loss in kg/m^2

          where n is is a month index from 1 to 12
  
  #### Intermediate
	  * ../test_out2/Intermediate/model_intermediate/WF/wf_<n>.tif: monthly Weather Factor in kg/m
	  * ../test_out2/Intermediate/model_intermediate/KK/KK_<n>.tif: monthly Surface Terrain Roughness unitless
	  * ../test_out2/Intermediate/model_intermediate/COG/COG_<n>.tif: monthly Vegeation Factor unitless
	  * ../test_out2/Intermediate/model_intermediate/EF.tif: Erodible Fraction unitless
	  * ../test_out2/Intermediate/model_intermediate/SCF.tif: Soil Crusting Factor unitless

	  where n is is a month index from 1 to 12
  

#### Runtime
Seconds.
  

  
