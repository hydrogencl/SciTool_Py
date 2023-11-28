#!/usr/bin/python
# Check the spin-up of CLM, using 8760 hr history fileG.

import numpy as np
from scipy.io import netcdf
import math , re

# Read File
def read_vars(str_path, str_filename1, arr_2dvar): 

  ncdf_file01 = "{0:s}/{1:s}".format(str_path,str_filename1)
  file01 = netcdf.netcdf_file( ncdf_file01 )

  arr_in = file01.variables["LONGXY"]
#  arr_in = np.zeros(( len(arr_2dvar),file01.dimensions[v1.dimensions[0]] , file01.dimensions[v1.dimensions[1]]))
  
  arr_out = np.zeros(( len(arr_2dvar),file01.dimensions[arr_in.dimensions[0]] , file01.dimensions[arr_in.dimensions[1]]))
#  print file01.dimensions[arr_in.dimensions[0]],file01.dimensions[arr_in.dimensions[1]],file01.dimensions[arr_in.dimensions[2]]
  num_lu = 0
  for k in range( len(arr_2dvar) ):
    arr_in  = file01.variables[ arr_2dvar[k] ]
    for j in range( file01.dimensions[arr_in.dimensions[0]] ):
      for i in range( file01.dimensions[arr_in.dimensions[1]]):
        arr_out[k][j][i] = arr_in[j][i] 
  return arr_out

def find_latlon2gp(arr_maps, arr_coor, nx, ny):
  # ["LONGXY","LATIXY","LONE","LONW","LATN","LATS"]
  tg_ix = arr_coor[0]
  tg_iy = arr_coor[1]

# Find LATIXY 
  tmp_iy = -1
  tmp_ix = -1
  for j in range(1,ny):
    if tg_iy -arr_maps[1][ j-1][0] >=0  and tg_iy - arr_maps[1][ j][0] < 0: tmp_iy = j
  if tmp_iy == -1: print("ERROR! IY not found")
    
# Find LONGXY based on LATIXY
  for i in range(1,nx):
    if tg_ix -arr_maps[0][tmp_iy][ i-1] >=0  and tg_ix - arr_maps[0][tmp_iy][ i] < 0: tmp_ix = i
  if tmp_ix == -1: print("ERROR! IX not found")

# Check edge

  tmp_edg_n = arr_maps[4][tmp_iy][tmp_ix]
  tmp_edg_s = arr_maps[5][tmp_iy][tmp_ix]
  tmp_edg_w = arr_maps[3][tmp_iy][tmp_ix]
  tmp_edg_e = arr_maps[2][tmp_iy][tmp_ix]
 
  # Check edge again....
  #check NS
  if tmp_edg_n - tg_iy > 0 and tmp_edg_s - tg_iy > 0:
    tmp_iy -= 1
    print("Modified due to edge check (too S)") 
  if tmp_edg_n - tg_iy < 0 and tmp_edg_s - tg_iy < 0:
    tmp_iy += 1
    print("Modified due to edge check (too N)") 
  #check NS
  if tmp_edg_w - tg_ix > 0 and tmp_edg_e - tg_ix > 0:
    tmp_ix -= 1
    print("Modified due to edge check (too E)") 
  if tmp_edg_w - tg_ix < 0 and tmp_edg_e - tg_ix < 0:
    tmp_ix += 1
    print("Modified due to edge check (too W)") 
  
  tmp_edg_n = arr_maps[4][tmp_iy][tmp_ix]
  tmp_edg_s = arr_maps[5][tmp_iy][tmp_ix]
  tmp_edg_w = arr_maps[3][tmp_iy][tmp_ix]
  tmp_edg_e = arr_maps[2][tmp_iy][tmp_ix]

# Draw Edge and center points

  print("""    
   W:{7:2.6f}: 
   | 
   | --------------------------|-- N:{4:2.6f}: 
   | center:({0:2.6f},{1:2.6f})    |
   | grid  :({2:4d},{3:4d})        |
   | --------------------------|-- S:{5:2.6f}: 
                               |
                               E:{6:2.6f}: 
  """.format(tg_ix, tg_iy, tmp_ix, tmp_iy, tmp_edg_n, tmp_edg_s, tmp_edg_e, tmp_edg_w ))


# Example #

gridinfo = read_vars("/data/sva_inputs/cropNRW/","surfdata_lai11_0160x0090.nc",["LONGXY","LATIXY","LONE","LONW","LATN","LATS"])
arr_coor=[\
[6.4473198,50.8658521],\
[6.2969924,50.9297879],\
[6.3041256,50.6219142],\
[6.3396,50.50305] ]
arr_plz = ["SE","ME","RO","WU"]
for n in range(4):
  print arr_plz[n]
  find_latlon2gp(gridinfo, arr_coor[n], 90, 160)




