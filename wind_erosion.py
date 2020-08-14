# this is the process to calclutate wind storm prevention service
# for more information ,you can refer the acrticle《Spatio-temporal variation of wind erosion in Inner Mongolia of China between 2001 and 2010》

import numpy as np
import gdal
import os

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

# days of month 1-12
mondays = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

endday = [30, 58, 89, 119, 150, 180, 211, 242, 272, 303, 333, 364]

m_id = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']

# elevation data (m)
dem = gdal.Open('dem.tif')
dem_a = dem.ReadAsArray(0, 0, dem.RasterXSize, dem.RasterYSize)
dem_a2 = (dem_a < 0) * 0 + (dem_a >= 0) * dem_a

# air pressure
P = 101.3 * (1 - 0.0255 * dem_a2 / 1000 * (6357 / (6357 + dem_a2 / 1000))) ** 5.256

# factor WF
for k in range(0, 12):
    # monthly average temperature(℃)
    tem_m_path = 'month_tem/' + str(k + 1) + '月气温_Pro.tif'
    tem_m = gdal.Open(tem_m_path)
    tem_m_a = tem_m.ReadAsArray(0, 0, tem_m.RasterXSize, tem_m.RasterYSize)

    # multiply 0.1 because the unit of source data is 0.1 ℃
    tem_m_a2 = (tem_m_a < -100) * 0.1 + (tem_m_a >= -100) * tem_m_a

    # monthly total precipitation(mm)
    prcp_path = 'month_prcp/' + str(k + 1) + '月降水量_Pro.tif'
    prcp = gdal.Open(prcp_path)
    prcp_a = prcp.ReadAsArray(0, 0, prcp.RasterXSize, prcp.RasterYSize)

    # multiply 0.1 because the unit of source data is 0.1 mm
    prcp_a2 = (prcp_a < 0) * 0 + (prcp_a >= 0) * prcp_a

    # monthly total solar radiation (MJ/m2)
    SOL_path = 'month_SOL/rad2000' + str(k + 1) + '_grd.tif'
    SOL = gdal.Open(SOL_path)
    SOL_a = SOL.ReadAsArray(0, 0, SOL.RasterXSize, SOL.RasterYSize)
    SOL_a2 = (SOL_a < 0) * 0.1 + (SOL_a >= 0) * SOL_a

    # days of rain events in every month
    prcp_days_path = 'month_prcp_day/prcp_day_' + str(k + 1) + '.tif'
    prcp_days = gdal.Open(prcp_days_path)
    prcp_days_a = prcp_days.ReadAsArray(0, 0, prcp_days.RasterXSize, prcp_days.RasterYSize)
    prcp_days_a2 = (prcp_days_a < 0) * 0 + (prcp_days_a >= 0) * prcp_days_a

    # potential evatranspiration (mm)
    ETp = 0.0135 * (SOL_a2 / 2.54) * (tem_m_a2 + 17.8)

    # soil moisture factor
    SW = (ETp - prcp_a2*prcp_days_a2) / ETp
    SW = (SW < 0) * 0 + (SW >= 0) * SW

    # snow factor
    SD_path = 'month_SD/SD_' + str(k + 1) + '.tif'
    SD = gdal.Open(SD_path)
    SD_a = SD.ReadAsArray(0, 0, SD.RasterXSize, SD.RasterYSize)
    SD_a2 = (SD_a < 0) * 0 + (SD_a >= 0) * SD_a


    if k == 0:
        end_date = int(endday[k]) + 1
        wf_sum = 0.0
        #print("Start Date: ", start_date)
        #print(" End Date: ", end_date)
        #print('\n')
        for i in range(0, end_date):

            # daily wind speed (m/s)
            wind_path = 'month_wind_day/wind_' + str(i + 1) + '.tif'
            wind = gdal.Open(wind_path)
            wind_a = wind.ReadAsArray(0, 0, wind.RasterXSize, wind.RasterYSize)

            # multiply 0.1 because the unit of source data is 0.1 m/s
            wind_a2 = (wind_a < 50) * 0 + ((wind_a >= 50) & (wind_a < 2000)) * wind_a * 0.1 + (wind_a >= 2000) * 200

            # air density
            air_den = 1.293 * (273 / (273 + tem_m_a2)) * P / 101.3

            # wind factor
            wf = wind_a2 * (wind_a2 - 5) ** 2
            # multiply 0.01 because the unit of source data is %
            SD_a2 = (1 - SD_a2 * 0.01)

            # days climate factor
            wf = wf * air_den / 9.8 * SW * SD_a2
            print(i, " : ",np.sum(SD_a2))
            
            wf_sum = wf_sum + wf
    else:
        start_date = int(endday[k - 1])
        end_date = int(endday[k]) + 1
        #print("Start Date: ", start_date)
        #print(" End Date: ", end_date)
        #print('\n')
        wf_sum = 0.0
        for i in range(start_date, end_date):

            # daily wind speed (m/s)
            wind_path = 'month_wind_day/wind_' + str(i + 1) + '.tif'
            wind = gdal.Open(wind_path)
            wind_a = wind.ReadAsArray(0, 0, wind.RasterXSize, wind.RasterYSize)

            # multiply 0.1 because the unit of source data is 0.1 m/s
            wind_a2 = (wind_a < 50) * 0 + ((wind_a >= 50) & (wind_a < 2000)) * wind_a * 0.1 + (wind_a >= 2000) * 200

            # air density
            air_den = 1.293 * (273 / (273 + tem_m_a2)) * P / 101.3

            # wind factor
            wf = wind_a2 * (wind_a2 - 5) ** 2
            
            # multiply 0.01 because the unit of source data is %
            SD_a2 = (1 - SD_a2 * 0.01)

            # days climate factor
            wf = wf * air_den / 9.8 * SW * SD_a2
            print(i, " : ",np.sum(air_den))
            wf_sum = wf_sum + wf

    wf_out_path = 'output/WF/wf_' + str(k + 1) + '.tif'
    RasterSave(wf_sum, wf_out_path, dem)

#ratio of sand(%)
sand_path = 'soil/sand.tif'
sand = gdal.Open(sand_path)
sand_a = sand.ReadAsArray(0, 0, sand.RasterXSize, sand.RasterYSize)
sand_a2 = (sand_a <= 0) * 0.1 + (sand_a >0) * sand_a


#ratio of silt(%)
silt_path = 'soil/silt.tif'
silt = gdal.Open(silt_path)
silt_a = silt.ReadAsArray(0, 0, silt.RasterXSize, silt.RasterYSize)
silt_a2 = (silt_a <= 0) * 0.1 + (silt_a > 0) * silt_a


#ratio of orgnic matter(%)
om_path = 'soil/som.tif'
om = gdal.Open(om_path)
om_a = om.ReadAsArray(0, 0, om.RasterXSize, om.RasterYSize)
om_a2 = (om_a <= 0) * 0.1 + (om_a > 0) * om_a


#ratio of clay(%)
clay_path = 'soil/clay.tif'
clay = gdal.Open(clay_path)
clay_a = clay.ReadAsArray(0, 0, clay.RasterXSize, clay.RasterYSize)
clay_a2 = (clay_a <= 0) * 0.1 + (clay_a >0) * clay_a

# soil erodibility factor
EF = (29.09 + 0.31 * sand_a2 + 0.17 * silt_a2 + 0.33 * sand_a2 / clay_a2 - 2.59 * om_a2 - 0.95 * 0) / 100
RasterSave(EF,'output/EF.tif',dem)

# soil crust factor
SCF = 1 / (1 + 0.0066 * clay_a2 ** 2 + 0.021 * om_a2 ** 2)
RasterSave(EF,'output/SCF.tif',dem)
#

#
for i in range(0, 12):

    # read fraction of vegetation coverage data (%)
    fvc_path = 'FVC2015/A2015' + str(i+1) + '_fc.tif'
    fvc = gdal.Open(fvc_path)
    fvc_a = fvc.ReadAsArray(0, 0, fvc.RasterXSize, fvc.RasterYSize)
    fvc_a2 = ((fvc_a < 0) | (fvc_a > 100)) * 0 + ((fvc_a >= 0) & (fvc_a <= 100)) * fvc_a


    # vegetation factor
    COG =np.exp(-0.00483*fvc_a2)

    COG_path = 'output/COG/COG_' + str(i + 1) + '.tif'
    RasterSave(COG, COG_path, fvc)

SL_sum = 0.0
SL_p_sum = 0.0

# surface roughness factor it can calculated by the smith-Caarson equation
KK_path = 'KK/KK_1.tif'
KK = gdal.Open(KK_path)
KK_a = KK.ReadAsArray(0, 0, KK.RasterXSize, KK.RasterYSize)
KK_a2 = (KK_a <= 0) * 0.1 + (KK_a > 0) * KK_a

for i in range(0, 12):
    WF_path = 'output/WF/WF_' + str(i + 1) + '.tif'
    WF = gdal.Open(WF_path)
    WF_a = WF.ReadAsArray(0, 0, WF.RasterXSize, WF.RasterYSize)
    WF_a2 = (WF_a <= 0) * 0.1 + (WF_a > 0) * WF_a



    COG_path = 'output/COG/COG_' + str(i + 1) + '.tif'
    COG = gdal.Open(COG_path)
    COG_a = COG.ReadAsArray(0, 0, COG.RasterXSize, COG.RasterYSize)
    COG_a2 = (COG_a <= 0) * 0.1 + (COG_a > 0) * COG_a

    Qmax = 109.8 * (WF_a2 * EF * KK_a2 * SCF * COG_a2)
    S = 105.71 * (WF_a2 * EF * KK_a2 * SCF * COG_a2)**-0.3711
    SL = 100/(S*S+0.01)*Qmax * (np.exp(-(50 / (S + 0.01)) ** 2))

    # actual wind erosion (unit:kg/m2)
    SL=((SL<-100000)|(SL>100000))*0+((SL>=-100000)&(SL<=100000))*SL
    SL_path="output/SL/SL_" + str(i + 1) + '.tif'
    RasterSave(SL,SL_path,WF)
    SL_sum = SL_sum + SL

    Qmax_p = 109.8 * (WF_a2 * EF * KK_a2 * SCF)
    S_p = 105.71 * (WF_a2 * EF * KK_a2 * SCF)**-0.3711

    # potentionl wind erosion (unit:kg/m2)
    SL_p = 100/(S_p*S_p+0.01)*Qmax_p * ( np.exp(-(50 / (S_p+0.01)) ** 2))
    SL_p = ((SL_p < -100000) | (SL_p > 100000)) * 0 + ((SL_p >= -100000)& (SL_p <= 100000)) * SL_p
    SL_p_path = "output/SL_p/SL_p_" + str(i + 1) + '.tif'
    RasterSave(SL_p, SL_p_path, WF)
    SL_p_sum = SL_p_sum + SL_p

SL_path='output/SL.tif'
RasterSave(SL_sum, SL_path, WF)

SL_p_path='output/SL_p.tif'
RasterSave(SL_p_sum, SL_p_path, WF)

sand_re = SL_p_sum - SL_sum
sand_re = (sand_re < 0) * 0 + (sand_re >= 0) * sand_re

out_path = 'output/sand_re.tif'
RasterSave(sand_re, out_path, dem)
