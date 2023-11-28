#!/usr/bin/python
# A tool kit to do many things:
# 1. READING CLM into 1d
# 2. Resorting: array into different minutes or hours or daily
#               fix time frames. 
# 3. comparing: comparing between observation and simulation data, and then ready for calculation of correlation or regression.
# 4. 
# Z4: means using Z4's 2011-07 to 2012-07. That means there is 2 July and hence has length of 13 months.


import gridtrans as gt
import math,re
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.io import netcdf 
from scipy import stats
import numpy as np


# READING FILES AND ARRAY

def read_ncfile(str_path,str_file):
  nc_file = netcdf.netcdf_file("{0:s}/{1:s}".format (str_path, str_file) )
  return nc_file

def read_var_1d(str_varname, time, nc_file, num_y, num_x, if_multilayers ):
  var = nc_file.variables[str_varname]
  if if_multilayers:
    level = len(var[0])
    arr_new = np.zeros((level,time))
    for t in range(time):
      for l in range(level):
        arr_new[l][t] = var[t][l][num_y][num_x]
  else: 
    arr_new = np.zeros((time))
    for t in range(time):
      arr_new[t] = var[t][num_y][num_x]
  return arr_new

def chk_timestamp_tr32z4(in_string,ini_yr,if_inc_ic,gmt):
  ## Turning time into 10 digit: YYYYMMDDhhmm ##
  ## and providing timeslot in var_arr
  ## !! Customing for TR32 monitoring time format !! ##
  tmp_date = in_string

  tmp_date_arr=re.split('[\s-:T\/]',tmp_date)
  year   = int(tmp_date_arr[0])
  mon    = int(tmp_date_arr[1])
  day    = int(tmp_date_arr[2])
  timehh = int(tmp_date_arr[3])
  timemm = int(tmp_date_arr[4])
  timestamp = year*1E8 + mon * 1E6 + day * 1E4 + timehh * 1E2 + timemm
  timeslot  = (year-ini_yr) * 8760*6 + if_leap_arr(year,if_inc_ic)[mon]*6 + (day-1) * 24*6+ timehh*6 + int(timemm / 10) - gmt *6
  return timestamp,timeslot,year,mon,day

def time_slot_tr32z4(arr_time_read, arr_time_start):
  acc_month_h = [0, 0, 744, 1440, 2160, 2880, 3624, 4344, 5088, 5832, 6552, 7296, 8016, 8760]
  pos_slot =  \
    ( arr_time_read[0] - arr_time_start[0] ) * 8760 *2 + \
    ( acc_month_h [arr_time_read[1]] - acc_month_h [arr_time_start[1]]) * 2 + \
    ( arr_time_read[2] - arr_time_start[2] ) * 24 * 2 + \
    ( arr_time_read[3] - arr_time_start[3] ) * 2  + \
    ( arr_time_read[4] - arr_time_start[4] ) / 30
  
  return pos_slot

def read_csv_tr32_30min(str_fn , digi_start_yr , num_readyrs , num_gmt_shift , if_inc_intercalary):
  day_mon = [31,28,31,30,31,30,31,31,30,31,30,31]
  readf = open(str_fn, 'r')
  num_total_stp = num_readyrs * 8760 * 2
  #print num_total_stp
  arr_U  = np.zeros ( num_total_stp )
  arr_H  = np.zeros ( num_total_stp )
  arr_LE = np.zeros ( num_total_stp )
  arr_readall = readf.readlines()
  # year,DOY,Hr,Mn,DOY1,Hr1,Mn1,UStar_qc,dUStar_qc,UStar_flag,H_qc,dH_qc,H_flag,LvE_qc,dLvE_qc,LvE_flag,FCO2_qc,dFCO2,FCO2_flag,FmolCO2_qc,dFmolCO2,Dir,Mean
  # Read Numbers
  for j in range(1,len(arr_readall)-1 ):
    tmp_line = re.split(',', arr_readall[j].strip())
    year = int(   tmp_line[0]  )
    DOY  = int(   tmp_line[1]  )
    HH   = int(   tmp_line[2]  )
    MM   = int(   tmp_line[3]  )
    U    = float( tmp_line[7]  )
    H    = float( tmp_line[10] )
    LE   = float( tmp_line[13] )
    if if_inc_intercalary == True:
      num_timestep = ( year - digi_start_yr )*8760*2 + (DOY -1) * 24 * 2 + HH * 2  + int(float( MM) / 30.0) + num_gmt_shift * 2
    else:
      if DOY <= 58:
        num_timestep = ( year - digi_start_yr )*8760*2 + (DOY -1) * 24 * 2 + HH * 2  + int(float( MM) / 30.0) + num_gmt_shift * 2 

      elif DOY > 59:
        num_timestep = ( year - digi_start_yr )*8760*2 + (DOY -2) * 24 * 2 + HH * 2  + int(float( MM) / 30.0) + num_gmt_shift * 2
      else:
        pass
      # debug mode
      #print year,DOY,HH,MM
      #print j,num_timestep

      if num_timestep < 0:
        pass
      elif num_timestep - num_total_stp >= 0:
        pass
      else:
        #print num_timestep
        arr_U  [ num_timestep ] = U
        arr_H  [ num_timestep ] = H
        arr_LE [ num_timestep ] = LE
      #print HH,MM,H,LE
  return arr_U, arr_H, arr_LE

def read_val_z4(val_fn , arr_start , arr_read):
  day_mon = [31,28,31,30,31,30,31,31,30,31,30,31]
  
  file_01=netcdf.netcdf_file(val_fn)
  time_ntime = file_01.dimensions['ntime']
  arr_val_sh = file_01.variables['SH']
  arr_val_lh = file_01.variables['LH']
  arr_out = np.zeros((2,8760*2))
  pos_start = time_slot_tr32z4( arr_read, arr_start )
  print pos_start
  print arr_val_lh[time_ntime-1]
  for i in range(8760*2):
    print pos_start+i
    arr_out[0][i] = arr_val_lh[ pos_start + i ] 
    arr_out[1][i] = arr_val_sh[ pos_start + i ] 

  return arr_out

def read_sim_soil1loc(sim_fn,loc):
  # TEST!
  #wrt=open("test_lh_sim.log","w")
  #wrt.write("fcev fgev fctr lh         ")
  file_01=netcdf.netcdf_file(sim_fn)
  ttl_time    = file_01.dimensions['time']
  ttl_loc_num = file_01.dimensions['lat']  
  #print ttl_time,ttl_loc_num
  tsoi    = file_01.variables['TSOI']
  #tg      = file_01.variables['TG']
  h2osoi  = file_01.variables['H2OSOI']
  #print(len(fcev))
 
  arr_tsoi   = np.zeros((10 , ttl_time))
  #arr_tg     = np.zeros((ttl_time))
  arr_h2osoi = np.zeros((10 , ttl_time))
  i = loc
  for j in range(ttl_time):
    #arr_tg[j] = tg[j][i]
    for l in range(10):
      arr_tsoi[l][j] = tsoi[j][l][i]
      arr_h2osoi[l][j] = h2osoi[j][l][i]
  return arr_tsoi,arr_h2osoi

def read_sim_soil(sim_fn):
  # TEST!
  #wrt=open("test_lh_sim.log","w")
  #wrt.write("fcev fgev fctr lh         ")
  file_01=netcdf.netcdf_file(sim_fn)
  ttl_time    = file_01.dimensions['time']
  ttl_loc_num = file_01.dimensions['lat']  
  #print ttl_time,ttl_loc_num
  tsoi    = file_01.variables['TSOI']
  #tg      = file_01.variables['TG']
  h2osoi  = file_01.variables['H2OSOI']
  #print(len(fcev))
 
  arr_tsoi   = np.zeros((10, ttl_loc_num, ttl_time))
  #arr_tg     = np.zeros((ttl_time))
  arr_h2osoi = np.zeros((10, ttl_loc_num, ttl_time))

  for j in range(ttl_time):
    for i in range(ttl_loc_num):
    #arr_tg[j] = tg[j][i]
      for l in range(10):
        arr_tsoi[l][i][j] = tsoi[j][l][i]
        arr_h2osoi[l][i][j] = h2osoi[j][l][i]
  return arr_tsoi,arr_h2osoi

def extract_obs_30min_wu(arr_in,num_start_hour,num_hour_steps):
  new_arr = np.zeros(( num_hour_steps*2 ))
  for i in range(num_hour_steps*2):
    new_arr[i] = arr_in[ num_start_hour * 2 + i ]
  return new_arr

def extract_obs_30min(arr_in,num_start_hour,num_hour_steps):
  new_arr = np.zeros(( num_hour_steps*2 ))
  for i in range(num_hour_steps*2):
    new_arr[i] = arr_in[i*3 + num_start_hour*6]
  return new_arr

def extract_strict_obs_30min_wu(arr_all_in,num_spc_value,num_start_hour,num_hour_steps, num_empty):
  # Strictly extract values, excludes when forcing is -999.999
  # ex: [1,5,6,8,7,4]
  new_arr = np.zeros(( num_hour_steps*2 ))
  for i in range(num_hour_steps*2):
    # check if not valid
    chk = 1
    if arr_all_in[num_spc_value][num_start_hour*2 + i] == -9999:
      new_arr[i] = -999.999
    elif arr_all_in[num_spc_value][num_start_hour*2 + i] >= 2000:
      new_arr[i] = -999.999
    else:
      new_arr[i] = arr_all_in[num_spc_value][num_start_hour*2 + i]
  return new_arr

def extract_strict_obs_30min(arr_all_in,num_spc_value,num_start_hour,num_hour_steps, arr_chk_plc, num_empty):
  # Strictly extract values, excludes the ones with -999.999
  # Strict array means if one of the array is error or empty, then don't use it
  # ex: [1,5,6,8,7,4]
  new_arr = np.zeros(( num_hour_steps*2 ))
  num_plc = len(arr_chk_plc)
  for i in range(num_hour_steps*2):
    # check if not valid
    chk = 0
    for plc in range(num_plc):
      if arr_all_in[arr_chk_plc[plc]][i*3 + num_start_hour] == num_empty :
        chk = chk + 1
    if chk == 0:
      new_arr[i] = arr_all_in[num_spc_value][i*3 + num_start_hour*6]
    else:
      new_arr[i] = -999.999
  return new_arr

# RESORTING ARRAYS 

def shrarr(arr_old, num_start, num_end):
  num_total = num_end - num_start + 1
  arr_new = np.zeros((num_total))
  for i in range(num_total):
    arr_new[i] = arr_old[num_start + i]
  return arr_new

def resort_arr(arr_in, date_start, if_inc):
  # To resort array into Z4 definition
  arr_new = [[],[],[],[],[],[],[],[],[],[],[],[],[]]
  for m in range(13):
    for h in range(24):
      arr_new[m].append([])
 
  #print len(arr_new), len(arr_new[0]),len(arr_new[0][0])
  tmp_stt = re.findall(".",str(date_start))
  arr_num_stt = []
  for i in tmp_stt:
    arr_num_stt.append(int(i))

  arr_date_start  = [ arr_num_stt[0]*1000 + arr_num_stt[1]*100 + arr_num_stt[2]*10 + arr_num_stt[3],arr_num_stt[4]*10 + arr_num_stt[5],arr_num_stt[6] * 10 + arr_num_stt[7]    ]

  if if_inc:
    arr_dinmon  = [31,29,31,30,31,30,31,31,30,31,30,31]
    arr_daccmon = [0,31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366] 
    arr_haccmon = [0, 744, 1440, 2184, 2904, 3648, 4368, 5112, 5856, 6576, 7320, 8040, 8784]
  else:
    arr_dinmon  = [31,28,31,30,31,30,31,31,30,31,30,31]
    arr_daccmon = [0,31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]
    arr_haccmon = [0, 744, 1416, 2160, 2880, 3624, 4344, 5088, 5832, 6552, 7296, 8016, 8760]
  arr_here_dinmon = np.zeros((13)) 
  arr_here_daccmon = np.zeros((14)) 
  for m in range(13):
    if m == 0:
      arr_here_dinmon [ m ] = arr_dinmon[ arr_date_start[1] -1 ] - arr_date_start[2] +1
    elif m + arr_date_start[1]-1 < 12:
      arr_here_dinmon [ m ] = arr_dinmon [ m + arr_date_start[1]-1  ] 
    elif m == 12:
      arr_here_dinmon [ m ] = arr_date_start[2] -1
    else:
      arr_here_dinmon [ m ] = arr_dinmon [ m + arr_date_start[1]-1-12  ]  
  for m in range(13):
    arr_here_daccmon[m+1] = arr_here_daccmon[m] + arr_here_dinmon [ m ]

#  print arr_here_daccmon
#  print arr_here_dinmon
  for m in range(13):
    for d in range(int(arr_here_dinmon[m])):
      for h in range(24):
#      print ( int (arr_here_daccmon[m]) + d )*24 + h
        var_in = arr_in[ ( int (arr_here_daccmon[m]) + d )*24 + h  ]
        if var_in != -999.999:
          arr_new[m][h].append ( var_in ) 

 
  return arr_new

def resort_arr_z4(arr_in, date_start, if_inc):
  arr_new = [[],[],[],[]]
  for c in range(4):
    for m in range(13):
      arr_new[c].append([])
      for h in range(24):
        arr_new[c][m].append([])
  #print len(arr_new), len(arr_new[0]),len(arr_new[0][0])
  tmp_stt = re.findall(".",str(date_start))
  arr_num_stt = []
  for i in tmp_stt:
    arr_num_stt.append(int(i))

  arr_date_start  = [ arr_num_stt[0]*1000 + arr_num_stt[1]*100 + arr_num_stt[2]*10 + arr_num_stt[3],arr_num_stt[4]*10 + arr_num_stt[5],arr_num_stt[6] * 10 + arr_num_stt[7]    ]

  if if_inc:
    arr_dinmon  = [31,29,31,30,31,30,31,31,30,31,30,31]
    arr_daccmon = [0,31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366] 
    arr_haccmon = [0, 744, 1440, 2184, 2904, 3648, 4368, 5112, 5856, 6576, 7320, 8040, 8784]
  else:
    arr_dinmon  = [31,28,31,30,31,30,31,31,30,31,30,31]
    arr_daccmon = [0,31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]
    arr_haccmon = [0, 744, 1416, 2160, 2880, 3624, 4344, 5088, 5832, 6552, 7296, 8016, 8760]
  arr_here_dinmon = np.zeros((13)) 
  arr_here_daccmon = np.zeros((14)) 
  for m in range(13):
    if m == 0:
      arr_here_dinmon [ m ] = arr_dinmon[ arr_date_start[1] -1 ] - arr_date_start[2] +1
    elif m + arr_date_start[1]-1 < 12:
      arr_here_dinmon [ m ] = arr_dinmon [ m + arr_date_start[1]-1  ] 
    elif m == 12:
      arr_here_dinmon [ m ] = arr_date_start[2] -1
    else:
      arr_here_dinmon [ m ] = arr_dinmon [ m + arr_date_start[1]-1-12  ]  
  for m in range(13):
    arr_here_daccmon[m+1] = arr_here_daccmon[m] + arr_here_dinmon [ m ]

  print arr_here_daccmon
#  print arr_here_dinmon
  for c in range(4):
    for m in range(13):
      for d in range(int(arr_here_dinmon[m])):
        for h in range(24):
#        print ( int (arr_here_daccmon[m]) + d )*24 + h
          arr_new[c][m][h].append ( arr_in[c][ ( int (arr_here_daccmon[m]) + d )*24 + h  ] )
  return arr_new

def arr_shrink2timeframe(arr_in, time_start, time_end, if_multilayer=False):
  time_ttl = time_end - time_start + 1
  if if_multilayer == True:
    len_v = len(arr_in)
    arr_out = [[0 for t in range(time_ttl)] for v in range(len_v)]
    for v in range(len_v):
      for t in range(time_ttl):
        arr_out[v][t] = arr_in[v][time_start+t]
  else:
    arr_out = [0 for t in range(time_ttl)]
    for t in range(time_ttl):
      arr_out[t] = arr_in[time_start+t]
  return arr_out 

def chk_num_null(arr_in, num_out=float("NaN"), num_null=-999.999):
  arr_out = [num_out for t in range(len(arr_in))]
  for t in range(len(arr_in)):
    if arr_in[t] != num_null:
      arr_out[t]=arr_in[t]
  return arr_out

def arr_ec_hrly_csv(arr_in_csv,arr_var_in):
  # only read in R_sw, LH, SH, G
  leng_ttl = len(arr_in_csv[0])
  leng_var = len(arr_var_in)
  print "length of hr",leng_ttl
  arr_out1 = [[0 for h in range(leng_ttl)] for var in range(leng_var)]
  arr_out2 = [[0 for h in range(8760)] for var in range(leng_var)]
  for v in range(leng_var):
    num_v = arr_var_in[v]
    for t in range(leng_ttl):
      arr_out1[v][t]=arr_in_csv[num_v][t]
  for v in range(leng_var):
    arr_out2[v] = arr10m1hr(arr_out1[v])
  return arr_out1, arr_out2

# CALCULATING ARRAY

def avg_arr(arr_in):
  arr_out = np.zeros((13,24))
  for m in range(13):
    for h in range(24):
      arr_out[m][h]= sum(arr_in[m][h]) / float( len(arr_in[m][h]))
  return arr_out

def avg_arr_z4(arr_in):
  arr_out = np.zeros((4,13,24))
  for c in range(4):
    for m in range(13):
      for h in range(24):
        arr_out[c][m][h]= sum(arr_in[c][m][h]) / float( len(arr_in[c][m][h]))
  return arr_out

def avg30min1hr(arr_old):
  length = len(arr_old)/2
  arr_new = np.zeros((length))
  
  for i in range(length):
    if arr_old[i*2] + arr_old[i*2+1] < - 100:
      arr_new[i] = 0.0
    elif arr_old[i*2]  < -100 : 
      arr_new[i] = arr_old[i*2 +1] 
  
    elif arr_old[i*2+1] < -100 :
      arr_new[i] = arr_old[i*2] 
    else:
      arr_new[i] = ( arr_old[i*2] + arr_old[i*2 +1] ) / 2.0
  return arr_new  

def arr10m1hr(arr_in,if_avg=True, num_null=-999.999):
  length_in  = len(arr_in)
  length_out = length_in / 6
  
  arr_out=[0 for t in range(length_out)]
  for t in range(length_out):
    tmp1 = arr_in[ t*6   ] 
    tmp2 = arr_in[ t*6+3 ]
    weight1 = 1
    weight2 = 1
    if tmp1 == num_null : weight1= 0
    if tmp2 == num_null : weight2= 0
    if (weight1 + weight2) == 0:
      tmp3 = 0
    else:
      tmp3 = (tmp1 * weight1 + tmp2 * weight2 ) / (weight1 + weight2)
    if if_avg == True: tmp3=tmp3/2.0
    if tmp3 == 0:
      arr_out[t] = num_null
    else:
      arr_out[t] = tmp3
  return arr_out 

def hourly2daily(arr_in, mode_trans="avg", if_strik=False, num_null=-999.999):

  length_out = len(arr_in) / 24
  arr_out    = [0 for d in range(length_out)]
  for t in range(length_out):
    num_tmp = 0
    num_wei = 0
    for h in range(24):
      tmpin = arr_in[t*24 + h]
      if tmpin != num_null:
        num_tmp += tmpin
        num_wei += 1.0
    if if_strik==False:
      if mode_trans=="avg":
        if num_wei==0.0:
          arr_out[t]=num_null
        else:
          arr_out[t]=num_tmp/num_wei
      elif mode_trans=="cum":
        arr_out[t]=num_tmp
    else:
      if num_wei <= 24:
        arr_out[t] = num_null
      else:
        if mode_trans=="avg":
          arr_out[t]=num_tmp/num_wei
        elif mode_trans=="cum":
          arr_out[t]=num_tmp
  return arr_out

def fix_arr_sp_WC2SAT(arr_in,num_level,num_loc,num_SWC ):
  if num_level > 1:
    num_arr_len = len(arr_in[0][0])
    arr_out = np.zeros((num_level,num_loc,num_arr_len))
    for l in range(num_level):
      for j in range(num_loc):
        for i in range(num_arr_len):
          arr_out[l][j][i] = arr_in[l][j][i] / float(num_SWC) * 100
  else:
    num_arr_len = len(arr_in)
    arr_out = np.zeros((num_loc,num_arr_len))
    for i in range(num_arr_len):
      for j in range(num_loc):
        arr_out[j][i] = arr_in[j][i] / float(num_SWC) * 100
  return arr_out

def fix_obs_arr_c2k (arr_in , num_target_col):
  for i in range(len(arr_in[num_target_col])):
    val_in = arr_in[num_target_col][i]
    if val_in == -999.999:
      arr_in[num_target_col][i] = arr_in[num_target_col][i]
      #print "find 999 empty"
    elif val_in == 0.0:
      arr_in[num_target_col][i] = arr_in[num_target_col][i]
      #print "find 000 empty at:",arr_in[0][i]," origin: ",val_in
    elif val_in > 100.0:
      arr_in[num_target_col][i] = arr_in[num_target_col][i]
      #print "find over boil temperature:",arr_in[0][i]," origin: ",val_in 
    else: 
      arr_in[num_target_col][i] = arr_in[num_target_col][i] + 273.15
  return arr_in

def fix_arr_sp_c2k(arr_in,num_level,num_loc ):
  if num_level > 0:
    #print "level:{0:d} ,loc:{1:d} ,time_len:{2:d}".format(len(arr_in),len(arr_in[0]),len(arr_in[0][0]))
    num_arr_len = len(arr_in[0][0])
    arr_out = np.zeros((num_level,num_loc,num_arr_len))
    for l in range(num_level):
      for j in range(num_loc):
        for i in range(num_arr_len):
          arr_out[l][j][i] = arr_in[l][j][i] + 273.15
  else:
    num_arr_len = len(arr_in)
    arr_out = np.zeros((num_loc,num_arr_len))
    for j in range(num_loc):
      for i in range(num_arr_len):
        arr_out[j][i] = arr_in[j][i] + 273.15
  return arr_out

def fix_arr_linear_interpola(arr_in1, arr_in2, a , b , c , num_loc):
  # ( Ar - Br ) / ( Ar - Cr ) = ( a - b ) / ( a - c)
  # Br = Ar - ( a-b)/ (a-c) * (Ar - Cr )

  num_arr_len = len(arr_in1) 
  num_arr_len_chk = len (arr_in2)
  if num_arr_len != num_arr_len_chk:
    print("Error! two arrays do not match their length")
  else:
    if num_loc >1:
      num_arr_len = len(arr_in1[0]) 
      arr_out = np.zeros((num_loc, num_arr_len))
      for j in range(num_loc):
        for i in range(num_arr_len):
          arr_out[j][i] = arr_in1[j][i] + (arr_in2[j][i] - arr_in1[j][i]) * (float((b-a)) / float((c-a)))
          #print arr_out[j][i]
    else:
      arr_out = np.zeros((num_arr_len))
      for i in range(num_arr_len):
        arr_out[i] = arr_in1[i] - (arr_in1[i] - arr_in2[i]) * ((a-b)/(a-c))
  return arr_out

# COMPARING and Calculation

def compare_obs_sim_mon(arr_in_obs, arr_in_sim, year=2012, num_hour_intvl=1, if_ic=False):
# Monthly array
  #num_start_plc = (year - 2011) * 8760*2
  #num_end_stp   = 8760*2
  arr_sim = 0
  arr_obs =  0
  arr_out = [[[],[]]for m in range(12)]
  arr_valid = [0 for m in range(12)]
  if if_ic == True:
    arr_step_sta = [0,744, 1440, 2184, 2904, 3648, 4368, 5112, 5856, 6576, 7320, 8040, 8784]
    arr_step_sto = [744, 1440, 2184, 2904, 3648, 4368, 5112, 5856, 6576, 7320, 8040, 8784]
  else:
    arr_step_sta  = [0,744, 1416, 2160, 2880, 3624, 4344, 5088, 5832, 6552, 7296, 8016, 8760]
    arr_step_sto  = [744, 1416, 2160, 2880, 3624, 4344, 5088, 5832, 6552, 7296, 8016, 8760]

  for m in range(12):
    for i in range(arr_step_sta[m]* num_hour_intvl,arr_step_sto[m]* num_hour_intvl):
      if arr_in_obs[i] == -999.999:
        pass
      elif arr_in_obs[i] == -9999:
        pass
      elif arr_in_obs[i] == 0.00:
        pass
      else:
        arr_out[m][0].append(arr_in_obs[i])
        arr_out[m][1].append(arr_in_sim[i])
        arr_valid[m] = arr_valid[m] + 1
    if len(arr_out[m])==0:
      arr_out[m][0].append(-999.999)
      arr_out[m][1].append(-999.999)
  #format of arr_out: [month][obs/sim][valid_value]
  return arr_out,arr_valid

def data_screen(arr_in,arr_bound_range,num_null):
  num_size = len(arr_in)
  arr_out = np.zeros(num_size)
  for i in range(num_size):
    val_in = arr_in[i]
    if val_in < bound_range[0] :
      pass
    elif val_in > bound_range[1] :
      pass
    elif val_in == num_null :
      pass
    else:
      arr_out[i] = val_in

  return arr_out

def compare_obs_sim_mon_z4(arr_in_sim, arr_in_obs, if_z4=False):
# Monthly array
  if_ic = False
  arr_sim = 0
  arr_obs =  0
  arr_out = [[[],[]],[[],[]],[[],[]],[[],[]],[[],[]],[[],[]],[[],[]],[[],[]],[[],[]],[[],[]],[[],[]],[[],[]]]
  arr_valid = np.zeros((12))
  if if_ic :
    if if_z4:
      arr_step_sta = [  0,  11,  42,  72, 103, 133, 164, 195, 223, 254, 284, 315, 345, 365]
      arr_step_sto = [ 11,  42,  72, 103, 133, 164, 195, 223, 254, 284, 315, 345, 365, 0]
    else:
      arr_step_sta = [0,744, 1440, 2184, 2904, 3648, 4368, 5112, 5856, 6576, 7320, 8040, 8784]
      arr_step_sto = [744, 1440, 2184, 2904, 3648, 4368, 5112, 5856, 6576, 7320, 8040, 8784]

  else:
    if if_z4:
      arr_step_sta = [  0,  11,  42,  72, 103, 133, 164, 195, 223, 254, 284, 315, 345, 365]
      arr_step_sto = [ 11,  42,  72, 103, 133, 164, 195, 223, 254, 284, 315, 345, 365, 0]
    else:
      arr_step_sta  = [0,744, 1416, 2160, 2880, 3624, 4344, 5088, 5832, 6552, 7296, 8016, 8760]
      arr_step_sto  = [744, 1416, 2160, 2880, 3624, 4344, 5088, 5832, 6552, 7296, 8016, 8760]

  for m in range(12):
    for i in range(arr_step_sta[m],arr_step_sto[m]):
      if arr_in_obs[i] == -999.999:
        pass
      elif arr_in_obs[i] == -9999:
        pass
      elif arr_in_obs[i] == 0.00:
        pass
      elif arr_in_obs[i] > 1000:
        print("pass_log: over  1000: {0:d} {1:5.8f} {2:5.8f}".format(i,arr_in_obs[i],arr_in_sim[i]))
        pass
      elif arr_in_obs[i] < -500:
        pass
      else:
        arr_out[m][0].append(arr_in_sim[i])
        arr_out[m][1].append(arr_in_obs[i])
        arr_valid[m] = arr_valid[m] + 1
        #print "compare_log:",i
    if len(arr_out[m])==0:
      arr_out[m][0].append(-999.999)
      arr_out[m][1].append(-999.999)
  #format of arr_out: [month][obs/sim][valid_value]
  return arr_out,arr_valid

def compare_obs_sim_yr (arr_in_sim, arr_in_obs ):
  arr_sim = 0
  arr_obs =  0
  arr_out = [[],[]]
  for i in range(8760):
    if arr_in_obs[i] == -999.999:
      pass
    elif arr_in_obs[i] == -9999:
      pass
    elif arr_in_obs[i] == 0.00:
      pass
    elif arr_in_obs[i] > 1000:
      print("pass_log: over  1000: {0:d} {1:5.8f} {2:5.8f}".format(i,arr_in_obs[i],arr_in_sim[i]))

    elif arr_in_obs[i] < -500:
      print("pass_log: under -500: {0:d} {1:5.8f} {2:5.8f}".format(i,arr_in_obs[i],arr_in_sim[i]))

    else:
      arr_out[0].append(arr_in_sim[i])
      arr_out[1].append(arr_in_obs[i])
  #    if arr_in_obs[i] > 700:
  #      print("compare_log:{0:d} {1:5.8f} {2:5.8f}".format(i,arr_in_obs[i],arr_in_sim[i]))
  #    elif arr_in_sim[i] > 700:
  #      print("compare_log:{0:d} {1:5.8f} {2:5.8f}".format(i,arr_in_obs[i],arr_in_sim[i]))
  #print "end of scheme"
  return arr_out

def cal_balance_in_3var(arr_RSW, arr_LH, arr_SH, arr_G, arr_RLW=[], if_silence=True, num_null=-999.999):
  # check arr_length:
  length_RSW = len( arr_RSW )
  length_RLW = len( arr_RLW )
  length_LH  = len( arr_LH  )
  length_SH  = len( arr_SH  )
  length_G   = len( arr_G   )
  
  chk1 = (length_RSW - length_LH)  
  chk2 = (length_SH -length_G) 

  if chk1 == 0 and chk2 == 0:
    if if_silence: print("Check length of Array: OK")
    # Start working on it
    arr_out1= [0 for t in range(length_RSW) ]
    arr_out2= [0 for t in range(length_RSW) ]
    arr_out3= [0 for t in range(length_RSW) ]
    for t in range(length_RSW):
      if arr_RSW[t] == num_null:
        if if_silence==False: print("arr_RSW is null")
      elif arr_LH[t] == num_null:
        if if_silence==False: print("arr_LH is null")
      elif arr_SH[t] == num_null:
        if if_silence==False: print("arr_SH is null")
      elif arr_G[t] == num_null:
        if if_silence==False: print("arr_G is null")
      else:
        if length_RLW == 0:
          arr_out1[t] = arr_RSW[t] - arr_LH[t] - arr_SH[t] - arr_G[t]
          if arr_RSW[t] == 0:
            arr_out2[t] = 0
          else:
            arr_out2[t] = arr_out1[t] / (arr_RSW[t])
          if t == 0:
            arr_out3[t] = arr_RSW[t] - arr_LH[t] - arr_SH[t] - arr_G[t]
          else:
            arr_out3[t] = arr_out3[t-1] + arr_RSW[t] - arr_LH[t] - arr_SH[t] - arr_G[t] 
          R_total = arr_RSW[t] 
        else:
          arr_out1[t] = arr_RSW[t] + arr_RLW[t] - arr_LH[t] - arr_SH[t] - arr_G[t]
          arr_out2[t] = arr_out1[t] / (arr_RSW[t] + arr_RLW[t])
          if t == 0:
            arr_out3[t] = arr_RSW[t] + arr_RLW[t] - arr_LH[t] - arr_SH[t] - arr_G[t]
          else:
            arr_out3[t] = arr_out3[t-1] + arr_RSW[t] + arr_RLW[t] - arr_LH[t] - arr_SH[t] - arr_G[t] 
          R_total = arr_RSW[t] + arr_RLW[t]
        if arr_out2[t] > 1.0 : 
          print t, arr_out1[t], R_total
          arr_out2[t] = 1.0
        if arr_out2[t] < -1.0:
          print t, arr_out1[t], R_total
          arr_out2[t] = -1.0

    return arr_out1,arr_out2,arr_out3

  else:
    if if_silence: 
      print("Check length of Array: Not ok")
      print("length of R, LH, SH, G:{0:d},{1:d},{2:d},{3:d}".format(length_RSW, length_LH, length_SH, length_G))

def statistic_2arr(arr_in1, arr_in2):
  #locals()[arr_1_in]
  # Find: Mean, median, and Standard Deviation
  arr_1_mean = np.mean( arr_in1 )
  arr_2_mean = np.mean( arr_in2 )
  arr_1_med  = np.median( arr_in1 ) 
  arr_2_med  = np.median( arr_in2 )
  arr_1_std  = np.std( arr_in1 )
  arr_2_std  = np.std( arr_in2 )

  # Find Correlation Coefficient
  arr_1_corrcoef = np.corrcoef( arr_in1, arr_in2 )
  arr_2_corrcoef = np.corrcoef( arr_in2, arr_in1 )
  if len(arr_1_corrcoef) == 0:
    return [arr_1_mean, arr_1_med, arr_1_std, 0],[arr_2_mean, arr_2_med, arr_2_std, 0]
  else:
    return [[arr_1_mean, arr_1_med, arr_1_std, arr_1_corrcoef[0][1]],[arr_2_mean, arr_2_med, arr_2_std, arr_2_corrcoef[0][1]]]
  

def regression_obs_sim(arr_1_in,arr_2_in):
  arr_line = [ [0,0] , [0,0] ]
  arr_line=np.zeros((2,2))
  slope, intercept, r_value, p_value, std_err = stats.linregress(arr_1_in,arr_2_in)
  arr_line[1][0]=min(arr_1_in) * slope + intercept
  arr_line[1][1]=max(arr_1_in) * slope + intercept
  arr_line[0][0]=min(arr_1_in)
  arr_line[0][1]=max(arr_1_in)
  rsqr = r_value ** 2
  arr_eqn = str("y = {0:2.4f} x + {1:2.4f}".format(slope,intercept))
  return arr_line, arr_eqn, rsqr, [slope, intercept, r_value, p_value, std_err]

def cal_modelperform (arr_obs , arr_sim , num_empty=-999.999):
  # Based on Vazquez et al. 2002 (Hydrol. Process.)
  num_arr = len(arr_obs)
  num_n_total = num_arr
  num_sum = 0
  num_obs_sum = 0
  if num_n_total == 0:
    EF = 0.0
    CD = 0.0
    RRMSE = 0.0
  else:
    for n in range( num_arr ):
      if arr_obs[n] == num_empty:
        num_n_total = num_n_total - 1
      else:
        num_sum = num_sum + ( arr_sim[n] - arr_obs[n] )  ** 2 
        num_obs_sum = num_obs_sum +  arr_obs[n]
    RRMSE = ( num_sum / num_n_total ) ** 0.5 *  ( num_n_total / num_obs_sum )
    obs_avg = num_obs_sum / num_n_total
    num_n_total = num_arr
    oo_sum = 0
    po_sum = 0
    for nn in range( num_arr ):
      if arr_obs[nn] == num_empty:
        num_n_total = num_n_total - 1
      else:
        oo_sum = oo_sum + ( arr_obs[nn] - obs_avg )  ** 2
        po_sum = po_sum + ( arr_sim[nn] - arr_obs[nn] ) ** 2
    EF = ( oo_sum - po_sum ) / oo_sum
    CD = oo_sum / po_sum
  return RRMSE,EF,CD, num_arr

def cal_modelperform_mon(str_ef,str_type,year,name_site,name_short,result_ef_1,result_ef_1_all,result_ef_v1,num_empty,if_print_mod):

  arr_result_yearly = np.zeros((3))
  arr_result_monthly = np.zeros((3,12))

  #if if_print_mod:
    #print("-------------------------------------")
    #print("{0:d} || {1:s} ||| {2:s} |||".format(year,name_site,str_ef))
    #print("yearly result ({0:s}):".format(str_type))

  arr_result_yearly  = cal_modelperform (result_ef_1_all[0],result_ef_1_all[1],num_empty)
  if if_print_mod:
    #print(" RRMSE      ;    EF       ;    CD        ")
    print(" {0:3.3f}   ;    {1:3.3f} ;    {2:3.3f}  ".format(arr_result_yearly[0],arr_result_yearly[1],arr_result_yearly[2]))
    #print("Monthly result:")
  for m in range(12):
    if result_ef_v1[m] == 0:
      arr_result_monthly[0][m],arr_result_monthly[1][m],arr_result_monthly[2][m] = [float("NaN"),float("NaN"),float("NaN")]
      #arr_result_monthly[0][m],arr_result_monthly[1][m],arr_result_monthly[2][m] = [0,0,0]
    else:
      arr_result_monthly[0][m],arr_result_monthly[1][m],arr_result_monthly[2][m]   = cal_modelperform(result_ef_1[m][0],result_ef_1[m][1] , num_empty)
    #if if_print_mod:
      #print("MON: {0:d} , RRMSE: {1:3.3f} , EF: {2:3.3f}, CD: {3:3.3f} ".format(m+1,arr_result_monthly[0][m],arr_result_monthly[1][m],arr_result_monthly[2][m]))

  return arr_result_yearly, arr_result_monthly

def print_regres_mon(str_ef,str_type,year,name_site,name_short,result_ef_1,result_ef_1_all,result_ef_v1):
  #print("{0:d} || {1:s} ||| {2:s} |||".format(year,name_site,str_ef))
  print("yearly result:")

  arr_line01,arr_eqn01,r  = regression_obs_sim(result_ef_1_all[0],result_ef_1_all[1])
  print("EQN for {1:s}: {0:s}".format(arr_eqn01,str_type))
  #print("-------------------------------------")
  #print("Monthly result:")
  for m in range(12):
    if result_ef_v1[m] == 0:
      arr_line001 = [0,0]
      arr_eqn001  = " NONE "
    else:
      arr_line001,arr_eqn001,r  = regression_obs_sim(result_ef_1[m][0],result_ef_1[m][1])

def print_regres(str_ef,year,name_site,name_short,result_ef_1,result_ef_2,result_ef_3,result_ef_4,result_ef_1_all,result_ef_2_all,result_ef_3_all,result_ef_4_all,result_ef_v1,result_ef_v2,result_ef_v3,result_ef_v4):
  #print("{0:d} || {1:s} ||| {2:s} |||".format(year,name_site,str_ef))
  #print("yearly result:")

  arr_line01,arr_eqn01  = regression_obs_sim(result_ef_1_all[0],result_ef_1_all[1])
  arr_line02,arr_eqn02  = regression_obs_sim(result_ef_2_all[0],result_ef_2_all[1])
  arr_line03,arr_eqn03  = regression_obs_sim(result_ef_3_all[0],result_ef_3_all[1])
  arr_line04,arr_eqn04  = regression_obs_sim(result_ef_4_all[0],result_ef_4_all[1])
  print("EQN for BB-MO: {0:s}".format(arr_eqn01))
  print("EQN for JV-MO: {0:s}".format(arr_eqn02))
  print("EQN for BB-TG: {0:s}".format(arr_eqn03))
  print("EQN for JV-TG: {0:s}".format(arr_eqn04))
  #print("-------------------------------------")
  #print("Monthly result:")
  for m in range(12):
    if result_ef_v1[m] == 0:
      arr_line001 = [0,0]
      arr_eqn001  = " NONE "
    else:
      arr_line001,arr_eqn001  = regression_obs_sim(result_ef_1[m][0],result_ef_1[m][1])
    if result_ef_v2[m] == 0:
      arr_line002 = [0,0]
      arr_eqn002  = " NONE "
    else:
      arr_line002,arr_eqn002  = regression_obs_sim(result_ef_2[m][0],result_ef_2[m][1])
    if result_ef_v3[m] == 0:
      arr_line003 = [0,0]
      arr_eqn003  = " NONE "
    else:
      arr_line003,arr_eqn003  = regression_obs_sim(result_ef_3[m][0],result_ef_3[m][1])
    if result_ef_v4[m] == 0:
      arr_line004 = [0,0]
      arr_eqn004  = " NONE "
    else:
      arr_line004,arr_eqn004  = regression_obs_sim(result_ef_4[m][0],result_ef_4[m][1])

# PLOTING

def subplot_ohne_obs(fig, arr_in,title,unit_x,unit_y,arr_pos):
  ax01=fig.add_subplot(arr_pos[0],arr_pos[1],arr_pos[2])
  #if arr_pos[2] == 2 :
  #  ax01.set_title("".foramt(title,))
  #else: arr_pos[2] == 2 :
  ax01.set_title(arr_mon[arr_pos[2]-1])
#ax01.set_xlabel("Time (Hours)")
  ax01.plot(arr_in[0][arr_pos[2]-1] ,'g-')
  ax01.plot(arr_in[1][arr_pos[2]-1] ,'b-')
  ax01.plot(arr_in[2][arr_pos[2]-1] ,'r-')
  ax01.plot(arr_in[3][arr_pos[2]-1] ,'k-')
  if arr_pos[2] ==  4: ax01.set_ylabel(unit_y)
  if arr_pos[2] == 11: ax01.set_xlabel(unit_x)
  
  if arr_pos[2] ==  6: ax01.legend(('BB-Collatz','BB-Leuning','Jarvis-Stewart' ),bbox_to_anchor=(1.1, 0.85),prop={'size':10})

def subplot_mit_obs(fig, arr_in, arr_obs, title, unit_x, unit_y, arr_pos):
  ax01=fig.add_subplot(arr_pos[0],arr_pos[1],arr_pos[2])
  if arr_pos[2] == 2 :ax01.set_title(title)
#ax01.set_xlabel("Time (Hours)")
  ax01.plot(arr_in[0][arr_pos[2]-1] ,'g-')
  ax01.plot(arr_in[1][arr_pos[2]-1] ,'b-')
  ax01.plot(arr_in[2][arr_pos[2]-1] ,'r-')
  ax01.plot(arr_in[3][arr_pos[2]-1] ,'k-')
  ax01.plot(arr_obs[arr_pos[2]-1] ,'yo')
  if arr_pos[2] ==  4: ax01.set_ylabel(unit_y)
  if arr_pos[2] == 11: ax01.set_xlabel(unit_x)
  if arr_pos[2] ==  6: ax01.legend(('BB-Collatz','BB-Leuning','Jarvis-Stewart','Obs.' ),bbox_to_anchor=(1.1, 0.85),prop={'size':10})

def subplot_mon_regre(fig, arr_in, arr_in2,arr_line,arr_eqn,reg_r2,suptitle,title,unit_x,unit_y,arr_pos):
  #arr_in[..][0] = sim, [..][1] = obs
  #arr_mon = ["JUL.","AUG.","SEP.","OCT.","NOV.","DEC.","JAN.","FEB.","MAR.","APR.","MAY.","JUN."]
  arr_mon = ["JAN.","FEB.","MAR.","APR.","MAY.","JUN.","JUL.","AUG.","SEP.","OCT.","NOV.","DEC."]
  arr_color = ['#00FFFF','#0489B1','#0B2161','#2EFE2E','#3ADF00','#4B8A08','#886A08','#DF7401','#FF0000','#F4FA58','#F7D358','#FE9A2E']

  ax01=fig.add_subplot(arr_pos[0],arr_pos[1],arr_pos[2])
  ax01.set_title(title)
  ax01_max = np.amax(arr_in2)
  ax01_min = np.amin(arr_in2)
  #ax01_max = 900
  #ax01_min = -50
  for m in range(12):
  #  print len(arr_in[m][1]),len(arr_in[m][0])
    ax01.plot(arr_in[m][1],arr_in[m][0], 'o',color=arr_color[m],ms=5.0,markeredgewidth=0.1)

  ax01.plot(arr_line[0],arr_line[1],'r-')
  ax01.plot([ax01_min,ax01_max],[ax01_min,ax01_max],'m--')
  ax01.text(0.03,0.03,"Yearly: {0:s}, {2:s}={1:4.2f}".format(arr_eqn,reg_r2,r'$R^2$'),horizontalalignment='left',verticalalignment='center', transform = ax01.transAxes, fontsize=11)

  if arr_pos[2] == 1: ax01.set_ylabel(unit_y)
  if arr_pos[2] == 4: ax01.set_ylabel(unit_y)
  if arr_pos[2] == 4: ax01.set_xlabel(unit_x)
  if arr_pos[2] == 5: ax01.set_xlabel(unit_x)
  #if arr_pos[2] == 4: ax01.set_xlabel(unit_x)
  if arr_pos[2] == 6: ax01.set_xlabel(unit_x)
  
  if arr_pos[2] == 3: ax01.legend(arr_mon,bbox_to_anchor=(1.3, 1.00),prop={'size':10})

def plot_line_energy(fig, arr_pos, arr_fluxes, arr_fluxname, str_case, str_ylabel="(W/m**2)", str_xlabel="hours (hr)", arr_lnfmt=["r-","b-","g-","y-"]):  
  ax = fig.add_subplot (arr_pos[0],arr_pos[1], arr_pos[2])
  
  for n in range(len(arr_fluxes)):
    ax.plot(chk_num_null(arr_fluxes[n]),arr_lnfmt[n])
  ax.set_xlim((0,len(arr_fluxes[0])))
  ax.set_title(str_case)
  ax.set_ylabel(str_ylabel)
  ax.set_xlabel(str_xlabel)
  ax.legend(arr_fluxname,bbox_to_anchor=(1.0, 1.00),prop={'size':10} ) 

def plot_line_energy_twx(fig, arr_pos, arr_fluxes, arr_fluxes2, arr_fluxname, arr_fluxname2, str_case, str_ylabel="(W/m**2)", str_ylabel2="%", str_xlabel="hours (hr)", arr_lnfmt=["r-","b-","g-","y-"], arr_lnfmt2=["r-","b-","g-","y-"]):  
  ax   = fig.add_subplot (arr_pos[0],arr_pos[1], arr_pos[2])
  axy2 = ax.twinx()
  for n in range(len(arr_fluxes)):
    print n, len(arr_fluxes)
    ax.plot(chk_num_null(arr_fluxes[n]),arr_lnfmt[n])
  for nn in range(len(arr_fluxes2)):
    print nn, len(arr_fluxes2)
    axy2.plot(chk_num_null(arr_fluxes2[n]),arr_lnfmt2[n])
  ax.set_title(str_case)
  ax.set_ylabel(str_ylabel)
  ax.set_xlabel(str_xlabel)
  ax.legend(arr_fluxname,bbox_to_anchor=(1.0, 1.00),prop={'size':10} ) 




