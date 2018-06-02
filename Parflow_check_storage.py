import numpy
import gridtrans
import math

# 
# Simply check the storage difference for spin-up
# Compare by the begining and the ending
# Auther: YS
# Must lib: gridtrans, numpy


#fld      = "/data/CRUR_out/pfclm_bt/"
#job_name = "crur_bt_yr01_sp"

def read_pf_para(file_specific_storage, file_porosity, file_slope_x, file_slope_y, str_path=".", arr_thickness=[]):
  specific_store = gridtrans.readpfb("{0:s}/{1:s}".format(str_path, file_specific_storage))
  Porosity_Field = gridtrans.readpfb("{0:s}/{1:s}".format(str_path, file_porosity))
  
  slope_x   = gridtrans.readpfb("{0:s}/{1:s}".format(str_path, file_slopex))[0]
  slope_y   = gridtrans.readpfb("{0:s}/{1:s}".format(str_path, file_slopey))[0]

  Domain_Size = [specific_store[1],specific_store[2],specific_store[3]]
  Domain_d    = [specific_store[4],specific_store[5],specific_store[6]]
  if len(arr_thickness) == 0:
    # Default thickness
    Thickness   = [1.35, 1.35, 1.35, 1.35, 1.35,\
                   1.35, 1.35, 1.35, 1.35, 1.35,\
                   1.35, 1.35, 1.35, 1.35, 1.35,\
                   1.35, 1.35, 1.35, 1.35, 1.35,\
                   1.00, 0.70, 0.50, 0.30, 0.20,\
                   0.13, 0.07, 0.05, 0.03, 0.02]
  return specific_store,Porosity_Field,Domain_Size,Domain_d,Thickness 

####### Storage #######
def storage_cal(Press,Satur,Porosity_Field,Domain_d,Ss_Field,Thickness):

  total_store = 0
  t_surf_store =0
  t_compr_store=0
  t_poro_store =0
  print Domain_d,Domain_Size
  
  for k in range(Domain_Size[2]):
    for j in range(Domain_Size[1]):
      for i in range(Domain_Size[0]):
        #print k,j,i
        #print Satur[k][j][i],
        #print Porosity_Field[k][j][i]
        #print Thickness[k]
        poro_store     = Satur[k][j][i] * Porosity_Field[k][j][i] *  Domain_d[0] * Domain_d[1] * Thickness[k]
        compr_store    = Press[k][j][i] * Ss_Field[k][j][i] * Satur[k][j][i] * Domain_d[0] *Domain_d[1]  * Thickness[k]
        t_poro_store = t_poro_store  + poro_store
        t_compr_store= t_compr_store + compr_store

        if k == Domain_Size[2] - 1: 
          surface_press = Press[Domain_Size[2] - 1][j][i]
          if surface_press > 0.0: 
            surf_store     = surface_press * Domain_d[0] * Domain_d[1]
            t_surf_store = t_surf_store + surf_store
  total_store = t_poro_store + t_surf_store + t_compr_store
  print "********** START OF INDIVIDUAL **********"
  print "Surface storage :",t_surf_store
  print "Pore storage :",t_poro_store
  print "Compressed storage :",t_compr_store
  print "Total Subsurface store :",t_compr_store + t_poro_store
  print "Total storage :",total_store
  print "********** END OF INDIVIDUAL **********"
  return total_store,t_poro_store,t_surf_store,t_compr_store

####### RUN-OFF  #######

def RUN_OFF(press,slope,Manning,Domain_d):
  runoff =  1.0 / Manning * ( abs(slope)   ** 0.5) * (press**(5.0/3)) * Domain_d[0]
  #print runoff
  return runoff

#print RUN_OFF(10,0.05,5.52E-6,1.0)

def RUN_OFF_cal(press,Domain_Size,x_slope,y_slope,time):
  k=Domain_Size[2]-1
  R = 0
  for j in range(Domain_Size[1]):
    for i in range(Domain_Size[0]):
        thi = press[Domain_Size[2] - 1][j][i] 
        if i == 0 and thi > 0:
          if x_slope[0][j][i] >0:
            RUNOFF_X=RUN_OFF(thi,x_slope[0][j][i],Manning, Domain_d) * time
            R=R+RUNOFF_X 
        if i == (Domain_Size[0] -1) and thi > 0:
          if x_slope[0][j][i] <0:
            RUNOFF_X=RUN_OFF(thi,x_slope[0][j][i],Manning, Domain_d) * time
            R=R+RUNOFF_X 
        if j == 0 and thi > 0:
          if y_slope[0][j][i] >0:
            RUNOFF_Y=RUN_OFF(thi,y_slope[0][j][i],Manning, Domain_d) * time
            R=R+ RUNOFF_Y
        if j == (Domain_Size[1] -1) and thi > 0:
          if y_slope[0][j][i] <0:
            RUNOFF_Y=RUN_OFF(thi,y_slope[0][j][i],Manning, Domain_d) * time
            R=R+ RUNOFF_Y
        
  print "total runoff:",R
  return R

# Calculate change of storage

def chk_delta_store(file_end_press,file_end_satur, file_start_press, file_start_satur, specific_store, Porosity_Field, Domain_d, Thickness, str_path="."):
  start_value_p = gridtrans.readpfb("{0:s}/{1:s}".format(str_path,file_start_press))
  start_value_s = gridtrans.readpfb("{0:s}/{1:s}".format(str_path,file_start_satur))
  stop_value_p  = gridtrans.readpfb("{0:s}/{1:s}".format(str_path,file_stop_press))
  stop_value_s  = gridtrans.readpfb("{0:s}/{1:s}".format(str_path,file_stop_satur))

  start_total_store,start_t_poro_store,start_t_surf_store,start_t_compr_store = storage_cal(start_value_p[0],start_value_s[0],Porosity_Field[0],Domain_d,specific_store[0],Thickness)
  stop_total_store, stop_t_poro_store, stop_t_surf_store, stop_t_compr_store  = storage_cal(stop_value_p[0],stop_value_s[0],Porosity_Field[0],Domain_d,specific_store[0],Thickness)
  delta_total =  stop_total_store - start_total_store
  relative_diff = delta_total / ((start_total_store+stop_total_store) / 0.5)
  return delta_total,relative_diff

def cal_in_fld(fld,job_name):
  print job_name 
  
  stop_value_yr01_p,stop_value_yr01_s,start_value_yr01_p,start_value_yr01_s  = read_pf_for_sp(fld,job_name)
  result01 = chk_delta_store(stop_value_yr01_p,stop_value_yr01_s,start_value_yr01_p,start_value_yr01_s,specific_store,Porosity_Field,Domain_d,Thickness)
  
  print "********** START OF CALCULATION **********"
  print "DELTA:",result01[0]
  print "Relative DELTA:",result01[1]
  print "********** END OF CALCULATION **********"
