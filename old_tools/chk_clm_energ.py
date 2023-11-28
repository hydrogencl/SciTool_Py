#!/usr/bin/python
# Use python and scipy to check if the CLM spun-up

import gridtrans as gt
from scipy.io import netcdf
import numpy as np

def num_of_time(num_stps, num_dif_h, num_y, num_m, num_d, num_h ,if_ic, if_silence=True):
  # dif_hr: out put interval
  # num_y,m,d,h: start Y,M,D,H for calculation
  if if_ic:
    arr_hourofyear = [0,0,744, 1440, 2184, 2904, 3648, 4368, 5112, 5856, 6576, 7320, 8040, 8784]
  else:
    arr_hourofyear = [0,0,744, 1416, 2160, 2880, 3624, 4344, 5088, 5832, 6552, 7296, 8016, 8760]

  start_hour = arr_hourofyear[num_m] + (num_d-1) * 24 + num_h 
  now_hour = start_hour + num_dif_h * num_stps

  # Give stamp of year:
  if now_hour >= 8760:
    yrs = now_hour/8760
    num_stp_y = num_y + 1
    now_hour_cal2mon = now_hour - yrs *8760
  else:
    yrs = 1 
    num_stp_y = num_y
    now_hour_cal2mon = now_hour
  # Give stamp of mon:
  for m in range(1,13):
    mhup = now_hour_cal2mon - arr_hourofyear[m]
    mhlo = arr_hourofyear[m+1] - now_hour_cal2mon
    #print mhup,mhlo
    if mhup >= 0 and mhlo > 0:
      num_stp_m = m 
    elif mhlo == 0:
      num_stp_m = 1 
  # Give stamp of day (hod = hour of day):

  now_hour_hod = now_hour_cal2mon - arr_hourofyear[num_stp_m]
  #print now_hour_hod ,now_hour_cal2mon ,arr_hourofyear[num_stp_m],num_stp_m
  num_stp_d = 1
  for d in range(31):
    mdup= now_hour_hod - d*24  
    mdlo= (d+1)*24 - now_hour_hod 
    #print mdup,mdlo
    if mdup >=0 and mdlo > 0:
      num_stp_d = d + 1

  # Give stamp of hour
  num_stp_h_tmp = now_hour - (num_stp_y - num_y) * 8760  - arr_hourofyear[num_stp_m] - (num_stp_d-1)*24
  num_stp_h = num_stp_h_tmp * 3600
  if if_silence==False:
    print now_hour_hod ,now_hour,( num_stp_y - num_y ) * 8760 ,arr_hourofyear[num_stp_m-1],(num_stp_d-1)*24, num_stp_m, num_stp_d, num_stp_h

  return num_stp_y, num_stp_m, num_stp_d, num_stp_h

def read_nfiles2avg(name_folder, name_filetitle, num_filesteps, arr_nzyx, arr_not, arr_read2dvar, if_mask ):
  # arr_not is the array of input for num_of_time
  # example for clm name: clmoas_jm_yr02_sp.clm2.h0.2012-06-21-00000.nc 
  # clmoas_jm_yr02_sp.clm2.h0.2011-07-26-00000.nc
  #total_arr = np.zeros((len(arr_read2dvar),num_filesteps,ny,nx))
  # old = record every var, steps, ny and nx 
  #tmp_avg_arr = np.zeros((len(arr_read2dvar),num_filesteps,arr_nzyx[1],arr_nzyx[2]))
  #avg_arr = np.zeros((len(arr_read2dvar),num_filesteps,arr_nzyx[1],arr_nzyx[2]))


  tmp_avg_arr = np.zeros((len(arr_read2dvar),arr_nzyx[1],arr_nzyx[2]))
  avg_arr = np.zeros((len(arr_read2dvar),arr_nzyx[1],arr_nzyx[2]))

  for n in range(num_filesteps):
    date_file = num_of_time(n,arr_not[0],arr_not[1],arr_not[2],arr_not[3],arr_not[4],arr_not[5])
    name_file = "{0:s}/{1:s}.clm2.h0.{2:04d}-{3:02d}-{4:02d}-{5:05d}.nc".format(name_folder, name_filetitle, date_file[0], date_file[1], date_file[2], date_file[3] )
    ncdfile = netcdf.netcdf_file(name_file)

    if n == 0 :
      if if_mask == True:
        arr_mask = ncdfile.variables['landmask']
    for var in range(len(arr_read2dvar)):
      arr_tmp_var = ncdfile.variables[arr_read2dvar[var]]
      #print("open_step:{0:02d}".format(n))
      for j in range(int(arr_nzyx[1])):
        for i in range(int(arr_nzyx[2])):
          if if_mask == True:
            if arr_mask[j][i] == 1 :
              tmp_avg_arr[var][j][i] = tmp_avg_arr[var][j][i] + abs(arr_tmp_var[0][j][i])            
            else:
              tmp_avg_arr[var][j][i] = tmp_avg_arr[var][j][i] + 0 
          else:  
            tmp_avg_arr[var][j][i] = tmp_avg_arr[var][j][i] + abs(arr_tmp_var[0][j][i])  
  for var in range(len(arr_read2dvar)):         
    for j in range(arr_nzyx[1]):
      for i in range(arr_nzyx[2]):
        avg_arr[var][j][i] = tmp_avg_arr[var][j][i] / float(num_filesteps)
      
  return avg_arr, arr_mask

def total_avg_mask(arr_var_yr0, arr_var_yr1, arr_mask, ny, nx):
  arr_diff = np.zeros((ny,nx))
  num_sumall = 0
  num_count_grids = 0
  if arr_mask == "":
    num_count_grids = num_count_grids + 1
    num_sumall = num_sumall + arr_var_yr1[j][i] - arr_var_yr0[j][i] 
    arr_diff[j][i] = arr_var_yr1[j][i] - arr_var_yr0[j][i]
        
    print "There is no mask"
  else:
    print "Masked average"
    for j in range(ny):
      for i in range(nx): 
        if arr_mask[j][i] == 1:
          num_count_grids = num_count_grids + 1
          num_sumall = num_sumall + abs(arr_var_yr1[j][i] - arr_var_yr0[j][i]) 
          arr_diff[j][i] = arr_var_yr1[j][i] - arr_var_yr0[j][i]
    num_value_avg = num_sumall / float(num_count_grids)
  return num_value_avg, arr_diff, num_sumall, num_count_grids 
