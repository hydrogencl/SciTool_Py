import re,os
from  WRF_TOOLS import NamelistCreater as WNC
from  WRF_TOOLS import Executors as RUN_EXECS
from  WRF_TOOLS import Tools 

import LightweightEnsembleFramework as LEF

from  subprocess import run
import argparse

###############
# Parsers
###############

parser = argparse.ArgumentParser(description='Reading input date')
parser.add_argument("-wg","--wps-geogrid", action=argparse.BooleanOptionalAction, 
                    dest='if_exec_wps_geogrid', default=False, help='if the script execute the geogrid.exe')
parser.add_argument("-wu","--wps-ungrib",  action=argparse.BooleanOptionalAction, 
                    dest='if_exec_wps_ungrib' , default=False, help='if the script execute the ungrib.exe ')
parser.add_argument("-wm","--wps-metgrid", action=argparse.BooleanOptionalAction, 
                    dest='if_exec_wps_metgrid', default=False, help='if the script execute the metgrid.exe')

parser.add_argument("-lg","--link_grib",   action=argparse.BooleanOptionalAction, 
                    dest='if_link_grib'       , default=False, help='if link the gribecute the ungrib.exe ')
parser.add_argument("-id","--input-data-set", type=str, 
                    dest='str_input_dataset'    , default='ECMWF', help='The choosing dataset, default: ECMWF ')

parser.add_argument("-r","--wrf-real", action=argparse.BooleanOptionalAction, 
                    dest='if_exec_wrf_real', default=False, help='if the script execute the real.exe')
parser.add_argument("-w","--wrf-wrf",  action=argparse.BooleanOptionalAction, 
                    dest='if_exec_wrf_wrf' , default=False, help='if the script execute the wrf.exe ')
parser.add_argument("-nlw","--namelist-wrf",  action=argparse.BooleanOptionalAction, 
                    dest='if_namelist_wrf' , default=False, help='if the script only create the namelist.input')
parser.add_argument("-ew","--ensemble-wrf",  action=argparse.BooleanOptionalAction, 
                    dest='if_ensemble_wrf' , default=False, help='if the script execute the wrf.exe ')
args = parser.parse_args()

###############
#   Namelist
###############

WNC = WNC()
WNC.debug = True

arrStartDateIn = [2022,  3,  3, 0, 0, 0]
arrEndDateIn   = [2022,  3,  3,12, 0, 0]

WNC.read_user_specific(run_time      = [0, 0, 0, 12, 0, 0],
                       start_time    = [2022,  3,  3, 0, 0, 0], 
                       end_time      = [2022,  3,  3,12, 0, 0],
                       max_dom       = 2, 
                       e_we          = [221, 341],
                       e_sn          = [221, 341],
                       dx            = 15000,
                       dy            = 15000,
                       grid_id       = [1,2],
                       parent_id     = [1,1],
                       i_parent_start = [1,94],
                       j_parent_start =  [1,54],
                       parent_grid_ratio = [1,5])

arrStartDate = [ Tools.make_wrf_datetime(arrStartDateIn), Tools.make_wrf_datetime(arrStartDateIn) ]
arrEndDate   = [ Tools.make_wrf_datetime(arrEndDateIn), Tools.make_wrf_datetime(arrEndDateIn) ]

# Parameters:

DIR_GEODATA="/p/project/cjiek80/ESIAS-MET/DATA/GEODATA"
DIR_GEFS   ="/gpfs/arch/jiek80/jiek8002/GEFS"
DIR_ECMWF  ="/p/fastdata/slmet/slmet111/met_data/ecmwf/era5/grib"

# Folders:

DIR_WPS    = "/p/scratch/cjiek80/jiek8010/WRFV4/WPS"
DIR_WRF    = "/p/scratch/cjiek80/jiek8010/WRFV4/WRF"
DIR_RUN    = "{0:s}/run".format(DIR_WRF)
# Running the geogrid
# Namelist of share

WNC.DIC_share_common_para["start_date"]["VALUE"] = arrStartDate
WNC.DIC_share_common_para["end_date"]["VALUE"]   = arrEndDate
WNC.DIC_share_common_para["wrf_core"]["VALUE"]   = "ARW"
WNC.DIC_share_common_para["max_dom" ]["VALUE"]   = 2

# Namelist of geogrid

WNC.DIC_geogrid_common_para["geog_data_res"]["VALUE"] = \
        [ 'modis_30s+2m', 'modis_30s', 'modis_30s' ]
WNC.DIC_geogrid_common_para["dx"]["VALUE"]                =  15000.0
WNC.DIC_geogrid_common_para["dx"]["VALUE"]                =  15000.0

WNC.DIC_geogrid_common_para["parent_id"]["VALUE"]         = [   1,   1 ]
WNC.DIC_geogrid_common_para["parent_grid_ratio"]["VALUE"] = [   1,   3 ]

WNC.DIC_geogrid_common_para["i_parent_start"]["VALUE"]    = [   1,  94 ]
WNC.DIC_geogrid_common_para["j_parent_start"]["VALUE"]    = [   1,  54 ]
WNC.DIC_geogrid_common_para["e_we"]["VALUE"]              = [ 221, 341 ] 
WNC.DIC_geogrid_common_para["e_sn"]["VALUE"]              = [ 221, 341 ]

WNC.DIC_geogrid_common_para["ref_lat"]["VALUE"]           = 54.0 
WNC.DIC_geogrid_common_para["ref_lon"]["VALUE"]           =  8.5
WNC.DIC_geogrid_common_para["truelat1"]["VALUE"]          = 30.0
WNC.DIC_geogrid_common_para["truelat2"]["VALUE"]          = 60.0
WNC.DIC_geogrid_common_para["stand_lon"]["VALUE"]         =  8.5

WNC.DIC_geogrid_common_para["geog_data_path"]["VALUE"]    = DIR_GEODATA  

##########################
#     WRF namelist
##########################

WNC.DIC_domains_common_para["num_metgrid_levels"]["VALUE"]    = 138  # ECMWF
WNC.DIC_domains_common_para["interp_type"]["VALUE"]           = 1  # ECMWF
WNC.DIC_domains_common_para["sfcp_to_sfcp"]["VALUE"]          = False
WNC.DIC_time_control_common_para["override_restart_timers"]["VALUE"] = True

WNC.DIC_time_control_common_para["history_outname"]["STR_FMT"] = "\'{0:s}_d<domain>{2:s}{1:s}\',"

WNC.DIC_domains_common_para["starting_time_step"]["VALUE"]     = [ 24, 12  ]
WNC.DIC_domains_common_para["max_time_step"]["VALUE"]          = [ 144, 24 ]
WNC.DIC_domains_common_para["time_step"]["VALUE"]              = [  60, 12] 
WNC.DIC_domains_common_para["time_step"]["ARR_TYPE"]           = "N"
WNC.DIC_domains_common_para["use_adaptive_time_step"]["VALUE"] = False

WNC.DIC_time_control_common_para["interval_seconds"]["VALUE"] = 10800

WNC.DIC_physics_common_para["mp_physics"]["VALUE"]         = 7
WNC.DIC_physics_common_para["cu_physics"]["VALUE"]         = [ 5 , 0, 0 ]
WNC.DIC_physics_common_para["bl_pbl_physics"]["VALUE"]     = 6 
WNC.DIC_physics_common_para["ra_sw_physics"]["VALUE"]      = 5 
WNC.DIC_physics_common_para["ra_lw_physics"]["VALUE"]      = 5 
WNC.DIC_physics_common_para["sf_sfclay_physics"]["VALUE"]  = 1 
WNC.DIC_physics_common_para["sf_surface_physics"]["VALUE"] = 2 


##########################
# Preparing the exec
##########################

RUN_EXE = RUN_EXECS()
RUN_EXE.strFolder_WPS = DIR_WPS
RUN_EXE.strFolder_WRF = DIR_WRF

# Execute geogrid.exe 
if args.if_exec_wps_geogrid:
    WNC.create_wps_namelist()
    RUN_EXE.run_geogrid()

# Execute the link of grib input
if args.if_link_grib:
    arrECMWF_Files = \
    [ "/p/fastdata/slmet/slmet111/met_data/ecmwf/era5/grib/2022/03/2022030300_ml.grb",
      "/p/fastdata/slmet/slmet111/met_data/ecmwf/era5/grib/2022/03/2022030303_ml.grb",     
      "/p/fastdata/slmet/slmet111/met_data/ecmwf/era5/grib/2022/03/2022030306_ml.grb",
      "/p/fastdata/slmet/slmet111/met_data/ecmwf/era5/grib/2022/03/2022030309_ml.grb",
      "/p/fastdata/slmet/slmet111/met_data/ecmwf/era5/grib/2022/03/2022030312_ml.grb"]
    RUN_EXE.link_ungrib(arrFiles=arrECMWF_Files)

#Execute ungrib.exe for ecmwf-era5
if args.if_exec_wps_ungrib and args.str_input_dataset == 'ECMWF':
    arrECMWF_Files = \
    [ "/p/fastdata/slmet/slmet111/met_data/ecmwf/era5/grib/2022/03/2022030300_ml.grb",
      "/p/fastdata/slmet/slmet111/met_data/ecmwf/era5/grib/2022/03/2022030303_ml.grb",     
      "/p/fastdata/slmet/slmet111/met_data/ecmwf/era5/grib/2022/03/2022030306_ml.grb",
      "/p/fastdata/slmet/slmet111/met_data/ecmwf/era5/grib/2022/03/2022030309_ml.grb",
      "/p/fastdata/slmet/slmet111/met_data/ecmwf/era5/grib/2022/03/2022030312_ml.grb"]
    RUN_EXE.link_ungrib(arrFiles=arrECMWF_Files)
    RUN_EXE.strTargetVtableName = 'Vtable.ECMWF.ml.grib2' 
    WNC.DIC_ungrib_common_para["prefix"]["VALUE"] = 'ECMWF_ML'
    WNC.DIC_metgrid_common_para["fg_name"]["VALUE"] = ["ECMWF_ML','ECMWF_SFC','PRES"]
    WNC.create_wps_namelist()
    RUN_EXE.run_ungrib()

    arrECMWF_Files = \
    [ "/p/fastdata/slmet/slmet111/met_data/ecmwf/era5/grib/2022/03/2022030300_sf.grb",
      "/p/fastdata/slmet/slmet111/met_data/ecmwf/era5/grib/2022/03/2022030303_sf.grb",     
      "/p/fastdata/slmet/slmet111/met_data/ecmwf/era5/grib/2022/03/2022030306_sf.grb",
      "/p/fastdata/slmet/slmet111/met_data/ecmwf/era5/grib/2022/03/2022030309_sf.grb",
      "/p/fastdata/slmet/slmet111/met_data/ecmwf/era5/grib/2022/03/2022030312_sf.grb"]
    RUN_EXE.link_ungrib(arrFiles=arrECMWF_Files)
    RUN_EXE.strTargetVtableName = 'Vtable.ERA-interim.ml' 
    WNC.DIC_ungrib_common_para["prefix"]["VALUE"] = 'ECMWF_SFC'
    WNC.create_wps_namelist()
    RUN_EXE.run_ungrib()

    WNC.create_wps_namelist()
    WNC.DIC_ungrib_common_para["prefix"]["VALUE"] = 'ECMWF_SFC'
    RUN_EXE.run_calc_ecmwf_p()

# Running the metgrid

if args.if_exec_wps_metgrid and args.str_input_dataset == 'ECMWF':
    WNC.DIC_metgrid_common_para["fg_name"]["VALUE"] = ["ECMWF_ML','ECMWF_SFC','PRES"]
    WNC.create_wps_namelist()
    RUN_EXE.run_metgrid()

if args.if_exec_wrf_real:
    WNC.create_wrf_namelist(STR_DIR = DIR_RUN)
    RUN_EXE.run_real()

if args.if_exec_wrf_wrf: 
    WNC.create_wrf_namelist(STR_DIR = DIR_RUN)
    RUN_EXE.run_wrf() 

if args.if_namelist_wrf: 
    WNC.create_wrf_namelist(STR_DIR = DIR_RUN)

if args.if_ensemble_wrf:
    SC = LEF.SlurmController(
            strProject      = "Test",
            strPartition    = "devel",
            numCoresPerNode = 48,
            strJobname      = "TestFEL",
            strRootdir      = "/p/scratch/cjiek80/jiek8010/WRFV4/WRF",
            ifServerLog     = True
            )
    
    SC.InitEnsemble(members=4, UsingNodes=1)
    print(SC.corespermember)
    print(SC.arr_hostpermember) 
    SC.CreateMembers()

    # Remove the 
    SC.FileControl("./wrfout_d01","remove")
    SC.FileControl("./wrfout_d02","remove")

    # run all the executor with the member
    SC.RunMembers("wrf.exe")

    # To keep this Server alive 
    SC.CheckMembersWRF()


