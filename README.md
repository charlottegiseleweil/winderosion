# Soil loss due to Wind Erosion Model

#### Requirements

* Sample data (can be downloaded [here](https://drive.google.com/file/d/1YQNOEnWyTdu1Po_D2klCpXNbAR3hdDgd/view?usp=sharing))
* A python environment with GDAL installed

#### Usage example
 ```python wind_erosion.py --data_dir=<Generated Result Data Dir> --input_data_dir=<Inputs Data Dir> --sl_output_suffix=<suffix> --dem_file_name=dem_aoi.tif --temper_dir=temperature_degC/ --temper_prefix=wc2.0_30s_tave_ --precip_dir=precipitation_mm/ --precip_prefix=wc2.0_30s_prec_ --sol_dir=solar_radiation/ --sol_prefix=shortwave_radiation_ --prcp_days_dir=month_prcp_day/ --prcp_days_prefix=prcp_day_ --snow_dir=snow_gm2/ --snow_prefix=snow_ --wind_speed_dir=wind_speed_monthly_clipped/ --wind_speed_prefix=wind_speed_ --sand_dir=soil/ --sand_file_name=sand.tif --silt_dir=soil/ --silt_file_name=silt.tif --clay_dir=soil/ --clay_file_name=clay.tif --som_dir=soil/ --som_file_name=soil_organic_matter_gm2.tif --fvc_dir=vegetation_percent_cover/ --fvc_prefix=vegetation_percent_cover_ ```

#### Inputs

  * input_data_dir is the location of the unzipped Sample Data
  * temper_dir/temper_prefix_<n>.tif: monthly average temperature in deg C
  * precip_dir/precip_prefix_<n>.tif: monthly total precipitation in mm
  * sol_dir/sol_prefix_<n>.tif: monthly solar radiation in MJ/m^2
  * prcp_days_dir/prcp_days_prefix_<n>.tif: monthly number of rain days 
  * snow_dir/snow_prefix_<n>.tif: monthly snow cover factor (probability)
  * wind_speed_dir/wind_speed_prefix_<n>.tif: monthly average wind speed in m/s
  * fvc_dir/fvc_prefix_<n>.tif: monthly fractional vegetation coverage in %
  * sand_dir/sand_file_name: non-temporal sand ratio (filename is a TIF file including extension)
  * silt_dir/silt_file_name: non-temporal silt ratio (filename is a TIF file including extension)
  * clay_dir/clay_file_name: non-temporal clay ratio (filename is a TIF file including extension)
  * som_dir/som_file_name: non-temporal organic matter ratio, data expected to be percentage*100 (filename is a TIF file including extension)
 
  where n is is a month index from 1 to 12

#### Outputs 
  * data_dir/Output/SL_(out_suffix}_<n>.tif: monthly Soil Loss in kg/m^2
  * data_dir/Output/SL_{out_suffix}.tif: Annual Soil Loss in kg/m^2
 
 where out_suffix is the argument value for --sl_output_suffix(default is "none")
   
#### Intermediate
  * data_dir/Intermediate/model_intermediate/WF/wf_<n>.tif: monthly Weather Factor in kg/m
  * data_dir/Intermediate/model_intermediate/KK/KK_<n>.tif: monthly Surface Terrain Roughness unitless
  * data_dir/Intermediate/model_intermediate/COG/COG_<n>.tif: monthly Vegeation Factor unitless
  * data_dir/Intermediate/model_intermediate/EF.tif: Erodible Fraction unitless
  * data_dir/Intermediate/model_intermediate/SCF.tif: Soil Crusting Factor unitless

#### Runtime
Seconds.
  

  
