# Soil loss due to Wind Erosion Model

#### Requirements

* Sample data (can be downloaded [here](https://drive.google.com/drive/folders/19JXm307FB47rQ3QBA0cwsKx263z4yW8e?usp=sharing))
* A python environment with GDAL installed

#### Usage example
 ```python wind_erosion.py --data_dir=../Data/sample --dem_file_name=dem_aoi.tif --temper_dir=temperature_degC/ --temper_prefix=wc2.0_30s_tave_ --precip_dir=precipitation_mm/ --precip_prefix=wc2.0_30s_prec_ --sol_dir=solar_radiation/ --sol_prefix=shortwave_radiation_ --prcp_days_dir=month_prcp_day/ --prcp_days_prefix=prcp_day_ --snow_dir=snow_gm2/ --snow_prefix=snow_ --wind_speed_dir=wind_speed_monthly_clipped/ --wind_speed_prefix=wind_speed_ --sand_dir=soil/ --sand_file_name=sand.tif --silt_dir=soil/ --silt_file_name=silt.tif --clay_dir=soil/ --clay_file_name=clay.tif --som_dir=soil/ --som_file_name=soil_organic_matter_gm2.tif --fvc_dir=vegetation_percent_cover/ --fvc_prefix=vegetation_percent_cover_ ```

#### Inputs

  * data_dir is the location of the unzipped Sample Data
  * temper_dir, precip_dir, sol_dir, prcp_days_dir, snow_dir, wind_speed_dir, fvc_dir are the directory names inside the Sample_Data/inputs folder where all the different monthly input data is located
  * <>_prefix are the prefix of each of the filenames, the code assumes files are named: prefix_1.tif(January), prefix_2.tif(February), etc
  * sand_dir, silt_dir, clay_dir, som_dir are the directory names for the non-temporal soil inputs

#### Outputs 

* [ EXPLICIT HERE WHAT ARE all the model outputs (with units!!)
   
#### Runtime
Seconds.
  

  
