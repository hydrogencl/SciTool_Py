import math,re,sys,os,time
import random as RD
import time
import struct

try:
    import netCDF4 as NC
except:
    print("You no install netCDF4 for python")
    print("So I do not import netCDF4")
try:
    import numpy as NP
except:
    print("You no install numpy")
    print("Do not import numpy")

class GRIDINFORMATER:
    """
    This object is the information of the input gridcells/array/map. 
    Using
    .add_an_element          to add an element/gridcell
    .add_an_geo_element      to add an element/gridcell
    .create_resample_lat_lon to create a new map of lat and lon for resampling
    .create_resample_map     to create resample map as ARR_RESAMPLE_MAP
    .create_reference_map    to create ARR_REFERENCE_MAP to resample target map.
    .export_reference_map    to export ARR_REFERENCE_MAP into netCDF4 format
    
    """
    STR_VALUE_INIT = "None"
    NUM_VALUE_INIT = -9999.9
    NUM_NULL       = float("NaN")

    ARR_RESAMPLE_X_LIM = []
    ARR_RESAMPLE_Y_LIM = []

    # FROM WRF: module_cam_shr_const_mod.f90
    NUM_CONST_EARTH_R = 6.37122E6 
    NUM_CONST_PI      = 3.14159265358979323846
 
    def __init__(self, name="GRID", ARR_LAT=[], ARR_LON=[], NUM_NT=1, DIMENSIONS=2 ):
        self.STR_NAME       = name
        self.NUM_DIMENSIONS = DIMENSIONS
        self.NUM_LAST_INDEX = -1
        self.ARR_GRID       = []
        self.NUM_NT         = NUM_NT       
        
        self.ARR_LAT  = ARR_LAT
        self.ARR_LON  = ARR_LON
        self.ARR_RESAMPLE_MAP_PARA = { "EDGE": {"N" :-999, "S":-999, "E":-999, "W":-999 } }

        if len(ARR_LAT) != 0 and len(ARR_LON) != 0:
            NUM_ARR_NY_T1 = len(ARR_LAT)
            NUM_ARR_NY_T2 = len(ARR_LON)
            Y_T2 = len(ARR_LON)
            NUM_ARR_NX_T1 = len(ARR_LAT[0])
            NUM_ARR_NX_T2 = len(ARR_LON[0])
            self.NUM_NX = NUM_ARR_NX_T1
            self.NUM_NY = NUM_ARR_NY_T1
            if NUM_ARR_NY_T1 - NUM_ARR_NY_T2 + NUM_ARR_NX_T1 - NUM_ARR_NX_T2 != 0:
                print("The gridcell of LAT is {0:d}&{1:d}, and LON is {2:d}&{3:d} are not match"\
                      .format(NUM_ARR_NY_T1,NUM_ARR_NY_T2,NUM_ARR_NX_T1,NUM_ARR_NX_T2))

    def index_map(self, ARR_IN=[], NUM_IN_NX=0, NUM_IN_NY=0):
        if len(ARR_IN) == 0:
            self.INDEX_MAP = [[ self.NUM_NULL for i in range(self.NUM_NX)] for j in range(self.NUM_NY)]
            NUM_ALL_INDEX = len(self.ARR_GRID)
            for n in range(NUM_ALL_INDEX):
                self.INDEX_MAP[self.ARR_GRID[n]["INDEX_J"]][self.ARR_GRID[n]["INDEX_I"]] =\
                self.ARR_GRID[n]["INDEX"]
        else:
            MAP_INDEX = [[ self.NUM_NULL for i in range(NUM_IN_NX)] for j in range(NUM_IN_NY)]
            NUM_ALL_INDEX = len(ARR_IN)
            for n in range(NUM_ALL_INDEX):
                MAP_INDEX[ARR_IN[n]["INDEX_J"]][ARR_IN[n]["INDEX_I"]] = ARR_IN[n]["INDEX"]
            return MAP_INDEX
                                                     
    def add_an_element(self, ARR_GRID, NUM_INDEX=0, STR_VALUE=STR_VALUE_INIT, NUM_VALUE=NUM_VALUE_INIT ):
        """ Adding an element to an empty array """
        OBJ_ELEMENT = {"INDEX" : NUM_INDEX, \
                       STR_VALUE : NUM_VALUE}
        ARR_GRID.append(OBJ_ELEMENT)
        
    def add_an_geo_element(self, ARR_GRID, NUM_INDEX=-999, NUM_J=0, NUM_I=0, \
                           NUM_NX = 0, NUM_NY = 0, NUM_NT=0, \
                           ARR_VALUE_STR=[], ARR_VALUE_NUM=[] ):
        """ Adding an geological element to an empty array 
            The information for lat and lon of center, edge, and vertex will
            be stored for further used. 
        """
        NUM_NVAR      = len(ARR_VALUE_STR)
        if NUM_NX == 0 or NUM_NY == 0:
            NUM_NX = self.NUM_NX 
            NUM_NY = self.NUM_NY
        if NUM_NT == 0:
            NUM_NT = self.NUM_NT
        NUM_CENTER_LON = self.ARR_LON[NUM_J][NUM_I] 
        NUM_CENTER_LAT = self.ARR_LAT[NUM_J][NUM_I]

        if NUM_I == 0:
            NUM_WE_LON =      ( self.ARR_LON[NUM_J][NUM_I] - self.ARR_LON[NUM_J][NUM_I + 1] ) * 0.5
            NUM_EW_LON = -1 * ( self.ARR_LON[NUM_J][NUM_I] - self.ARR_LON[NUM_J][NUM_I + 1] ) * 0.5
        elif NUM_I == NUM_NX - 1:
            NUM_WE_LON = -1 * ( self.ARR_LON[NUM_J][NUM_I] - self.ARR_LON[NUM_J][NUM_I - 1] ) * 0.5
            NUM_EW_LON =      ( self.ARR_LON[NUM_J][NUM_I] - self.ARR_LON[NUM_J][NUM_I - 1] ) * 0.5
        else:
            NUM_WE_LON = ( self.ARR_LON[NUM_J][NUM_I] - self.ARR_LON[NUM_J][NUM_I + 1] ) * 0.5
            NUM_EW_LON = ( self.ARR_LON[NUM_J][NUM_I] - self.ARR_LON[NUM_J][NUM_I - 1] ) * 0.5
        if NUM_J == 0:
            NUM_SN_LAT = -1 * ( self.ARR_LAT[NUM_J][NUM_I] - self.ARR_LAT[NUM_J + 1][NUM_I ] ) * 0.5
            NUM_NS_LAT =      ( self.ARR_LAT[NUM_J][NUM_I] - self.ARR_LAT[NUM_J + 1][NUM_I ] ) * 0.5
        elif NUM_J == NUM_NY - 1:
            NUM_SN_LAT =      ( self.ARR_LAT[NUM_J][NUM_I] - self.ARR_LAT[NUM_J - 1][NUM_I ] ) * 0.5 
            NUM_NS_LAT = -1 * ( self.ARR_LAT[NUM_J][NUM_I] - self.ARR_LAT[NUM_J - 1][NUM_I ] ) * 0.5
        else:       
            NUM_SN_LAT = ( self.ARR_LAT[NUM_J][NUM_I] - self.ARR_LAT[NUM_J - 1][NUM_I ] ) * 0.5
            NUM_NS_LAT = ( self.ARR_LAT[NUM_J][NUM_I] - self.ARR_LAT[NUM_J + 1][NUM_I ] ) * 0.5
        ARR_NE = [ NUM_CENTER_LON + NUM_EW_LON , NUM_CENTER_LAT + NUM_NS_LAT ]
        ARR_NW = [ NUM_CENTER_LON + NUM_WE_LON , NUM_CENTER_LAT + NUM_NS_LAT ] 
        ARR_SE = [ NUM_CENTER_LON + NUM_EW_LON , NUM_CENTER_LAT + NUM_SN_LAT ] 
        ARR_SW = [ NUM_CENTER_LON + NUM_WE_LON , NUM_CENTER_LAT + NUM_SN_LAT ] 
        if NUM_INDEX == -999: 
            NUM_INDEX = self.NUM_LAST_INDEX +1
            self.NUM_LAST_INDEX += 1
        OBJ_ELEMENT = {"INDEX"      : NUM_INDEX,\
                       "INDEX_I"    : NUM_I,\
                       "INDEX_J"    : NUM_J,\
                       "CENTER" : {"LAT" : NUM_CENTER_LAT, "LON" : NUM_CENTER_LON},\
                       "VERTEX" : {"NE": ARR_NE, "SE": ARR_SE, "SW": ARR_SW, "NW": ARR_NW},\
                       "EDGE"   : {"N": NUM_CENTER_LAT + NUM_NS_LAT,"S": NUM_CENTER_LAT + NUM_SN_LAT,\
                                   "E": NUM_CENTER_LON + NUM_EW_LON,"W": NUM_CENTER_LON + NUM_WE_LON}}
        if len(ARR_VALUE_STR) > 0:
            for I, VAR in enumerate(ARR_VALUE_STR):
                OBJ_ELEMENT[VAR] = [{ "VALUE" : 0.0} for t in range(NUM_NT)  ]
                if len(ARR_VALUE_NUM) == NUM_NVAR: 
                    for T in range(NUM_NT): 
                        OBJ_ELEMENT[VAR][T]["VALUE"] = ARR_VALUE_NUM[I][T]
        ARR_GRID.append(OBJ_ELEMENT)
        
    def add_an_geo_variable(self, ARR_GRID, NUM_INDEX=-999, NUM_J=0, NUM_I=0, NUM_NT=0,\
                           STR_VALUE=STR_VALUE_INIT, NUM_VALUE=NUM_VALUE_INIT ):
        if NUM_INDEX == -999:
            NUM_INDEX = self.INDEX_MAP[NUM_J][NUM_I]
        if NUM_NT == 0:
            NUM_NT = self.NUM_NT
        ARR_GRID[NUM_INDEX][STR_VALUE]  = {{"VALUE": NUM_VALUE } for t in range(NUM_NT)}
        
    def create_resample_lat_lon(self, ARR_RANGE_LAT=[0,0],NUM_EDGE_LAT=0,\
                                      ARR_RANGE_LON=[0,0],NUM_EDGE_LON=0 ):
        self.NUM_GRIDS_LON = round((ARR_RANGE_LON[1] - ARR_RANGE_LON[0])/NUM_EDGE_LON)
        self.NUM_GRIDS_LAT = round((ARR_RANGE_LAT[1] - ARR_RANGE_LAT[0])/NUM_EDGE_LAT)
        self.ARR_LAT = [[ 0 for i in range(self.NUM_GRIDS_LON)] for j in range(self.NUM_GRIDS_LAT) ]
        self.ARR_LON = [[ 0 for i in range(self.NUM_GRIDS_LON)] for j in range(self.NUM_GRIDS_LAT) ]
        for j in range(self.NUM_GRIDS_LAT):
            for i in range(self.NUM_GRIDS_LON):
                NUM_LAT = ARR_RANGE_LAT[0] + NUM_EDGE_LAT * j
                NUM_LON = ARR_RANGE_LON[0] + NUM_EDGE_LON * i
                self.ARR_LON[j][i] = ARR_RANGE_LON[0] + NUM_EDGE_LON * i
                self.ARR_LAT[j][i] = ARR_RANGE_LAT[0] + NUM_EDGE_LAT * j
                
    def create_reference_map(self, MAP_TARGET, MAP_RESAMPLE, STR_TYPE="FIX", NUM_SHIFT=0.001, IF_PB=False):
        """Must input with OBJ_REFERENCE
           WARNING: The edge of gridcells may not be included due to the unfinished algorithm

        """
        self.ARR_REFERENCE_MAP = []
        if STR_TYPE=="GRIDBYGEO":
            NUM_OBJ_G_LEN = len(MAP_TARGET)
            for OBJ_G in MAP_TARGET:
                NUM_G_COOR        = [OBJ_G["CENTER"]["LAT"], OBJ_G["CENTER"]["LON"]]
                for OBJ_R in MAP_RESAMPLE:
                    NUM_CHK_IN_EW = (OBJ_R["EDGE"]["E"] - OBJ_G["CENTER"]["LON"]) *\
                                    (OBJ_R["EDGE"]["W"] - OBJ_G["CENTER"]["LON"])
                    NUM_CHK_IN_SN = (OBJ_R["EDGE"]["N"] - OBJ_G["CENTER"]["LAT"]) *\
                                    (OBJ_R["EDGE"]["S"] - OBJ_G["CENTER"]["LAT"])
                    if NUM_CHK_IN_EW == 0: NUM_CHK_IN_EW = (OBJ_R["EDGE"]["E"] + NUM_SHIFT - OBJ_G["CENTER"]["LON"]) *\
                                                           (OBJ_R["EDGE"]["W"] + NUM_SHIFT - OBJ_G["CENTER"]["LON"])
                    if NUM_CHK_IN_SN == 0: NUM_CHK_IN_SN = (OBJ_R["EDGE"]["E"] + NUM_SHIFT - OBJ_G["CENTER"]["LON"]) *\
                                                           (OBJ_R["EDGE"]["W"] + NUM_SHIFT - OBJ_G["CENTER"]["LON"])
                    if NUM_CHK_IN_EW < 0 and NUM_CHK_IN_SN < 0:
                        OBJ_ELEMENT = {"INDEX"       : OBJ_G["INDEX"],\
                                       "CENTER"      : OBJ_G["CENTER"],\
                                       "INDEX_REF"   : OBJ_R["INDEX"],\
                                       "INDEX_REF_I" : OBJ_R["INDEX_I"],\
                                       "INDEX_REF_J" : OBJ_R["INDEX_J"],\
                                       "CENTER_REF"  : OBJ_R["CENTER"],\
                                             }          
                        self.ARR_REFERENCE_MAP.append(OBJ_ELEMENT)
                        break
                if IF_PB: TOOLS.progress_bar(TOOLS.cal_loop_progress([OBJ_G["INDEX"]], [NUM_OBJ_G_LEN]), STR_DES="CREATING REFERENCE MAP")
        elif STR_TYPE=="FIX":
            NUM_OBJ_G_LEN = len(MAP_TARGET)
            for OBJ_G in MAP_TARGET:
                NUM_G_COOR    = [OBJ_G["CENTER"]["LAT"], OBJ_G["CENTER"]["LON"]]
                if self.ARR_RESAMPLE_MAP_PARA["EDGE"]["W"] == -999 or self.ARR_RESAMPLE_MAP_PARA["EDGE"]["E"] == -999:
                    NUM_CHK_EW_IN = -1
                else:
                    NUM_CHK_EW_IN = (NUM_G_COOR[1] - self.ARR_RESAMPLE_MAP_PARA["EDGE"]["W"] ) * ( NUM_G_COOR[1] - self.ARR_RESAMPLE_MAP_PARA["EDGE"]["E"] )

                if self.ARR_RESAMPLE_MAP_PARA["EDGE"]["N"] == -999 or self.ARR_RESAMPLE_MAP_PARA["EDGE"]["S"] == -999:
                    NUM_CHK_SN_IN = -1
                else:
                    NUM_CHK_SN_IN = (NUM_G_COOR[0] - self.ARR_RESAMPLE_MAP_PARA["EDGE"]["S"] ) * ( NUM_G_COOR[0] - self.ARR_RESAMPLE_MAP_PARA["EDGE"]["N"] )
                if NUM_CHK_EW_IN < 0 and NUM_CHK_SN_IN < 0:            
                    for OBJ_R in MAP_RESAMPLE:
                        NUM_CHK_IN_EW = (OBJ_R["EDGE"]["E"] - OBJ_G["CENTER"]["LON"]) *\
                                        (OBJ_R["EDGE"]["W"] - OBJ_G["CENTER"]["LON"])
                        NUM_CHK_IN_SN = (OBJ_R["EDGE"]["N"] - OBJ_G["CENTER"]["LAT"]) *\
                                        (OBJ_R["EDGE"]["S"] - OBJ_G["CENTER"]["LAT"])
                        if NUM_CHK_IN_EW == 0: NUM_CHK_IN_EW = (OBJ_R["EDGE"]["E"] + NUM_SHIFT - OBJ_G["CENTER"]["LON"]) *\
                                                             (OBJ_R["EDGE"]["W"] + NUM_SHIFT - OBJ_G["CENTER"]["LON"])
                        if NUM_CHK_IN_SN == 0: NUM_CHK_IN_SN = (OBJ_R["EDGE"]["E"] + NUM_SHIFT - OBJ_G["CENTER"]["LON"]) *\
                                                             (OBJ_R["EDGE"]["W"] + NUM_SHIFT - OBJ_G["CENTER"]["LON"])
                        if NUM_CHK_IN_EW < 0 and NUM_CHK_IN_SN < 0:
                            OBJ_ELEMENT = {"INDEX"       : OBJ_G["INDEX"],\
                                           "INDEX_I"     : OBJ_G["INDEX_I"],\
                                           "INDEX_J"     : OBJ_G["INDEX_J"],\
                                           "CENTER"      : OBJ_G["CENTER"],\
                                           "INDEX_REF"   : OBJ_R["INDEX"],\
                                           "INDEX_REF_I" : OBJ_R["INDEX_I"],\
                                           "INDEX_REF_J" : OBJ_R["INDEX_J"],\
                                           "CENTER_REF"  : OBJ_R["CENTER"],\
                                                 }          
                            self.ARR_REFERENCE_MAP.append(OBJ_ELEMENT)
                            break
                if IF_PB: TOOLS.progress_bar(TOOLS.cal_loop_progress([OBJ_G["INDEX"]], [NUM_OBJ_G_LEN]), STR_DES="CREATING REFERENCE MAP")

    def export_grid_map(self, ARR_GRID_IN, STR_DIR, STR_FILENAME, ARR_VAR_STR=[],\
            ARR_VAR_ITEM=["MEAN", "MEDIAN", "MIN", "MAX", "P95", "P75", "P25", "P05"],\
            NUM_NX=0, NUM_NY=0, NUM_NT=0, STR_TYPE="netCDF4", IF_PB=False ):
        TIME_NOW = time.gmtime()
        STR_DATE_NOW = "{0:04d}-{1:02d}-{2:02d}".format(TIME_NOW.tm_year, TIME_NOW.tm_mon, TIME_NOW.tm_mday)
        STR_TIME_NOW = "{0:04d}:{1:02d}:{2:02d}".format(TIME_NOW.tm_hour, TIME_NOW.tm_min, TIME_NOW.tm_sec)
    
        if NUM_NX==0: NUM_NX = self.NUM_NX
        if NUM_NY==0: NUM_NY = self.NUM_NY
        if NUM_NT==0: NUM_NT = self.NUM_NT
    
        if STR_TYPE == "netCDF4":
            NCDF4_DATA = NC.Dataset("{0:s}/{1:s}".format(STR_DIR, STR_FILENAME), 'w', format="NETCDF4")
            # CREATE ATTRIBUTEs:
            NCDF4_DATA.description = \
            "The grid information in netCDF4"
            NCDF4_DATA.history = "Create on {0:s} at {1:s}".format(STR_DATE_NOW, STR_TIME_NOW)
    
            # CREATE DIMENSIONs:
            NCDF4_DATA.createDimension("Y"     , NUM_NY )
            NCDF4_DATA.createDimension("X"     , NUM_NX )
            NCDF4_DATA.createDimension("Time"  , NUM_NT )
            NCDF4_DATA.createDimension("Values", None   )
    
            # CREATE BASIC VARIABLES:
            NCDF4_DATA.createVariable("INDEX",          "i4", ("Y", "X"))
            NCDF4_DATA.createVariable("INDEX_J",        "i4", ("Y", "X"))
            NCDF4_DATA.createVariable("INDEX_I",        "i4", ("Y", "X"))
            NCDF4_DATA.createVariable("CENTER_LON",     "f8", ("Y", "X"))
            NCDF4_DATA.createVariable("CENTER_LAT",     "f8", ("Y", "X"))
           # CREATE GROUP for Variables: 
            for VAR in ARR_VAR_STR:
                NCDF4_DATA.createGroup(VAR)
                for ITEM in ARR_VAR_ITEM:
                    if ITEM == "VALUE" :
                        NCDF4_DATA.groups[VAR].createVariable(ITEM, "f8", ("Time", "Y", "X", "Values"))
                    else:
                        NCDF4_DATA.groups[VAR].createVariable(ITEM, "f8", ("Time", "Y", "X"))
            # WRITE IN VARIABLE
            
            for V in ["INDEX", "INDEX_J", "INDEX_I"]:
                map_in = self.convert_grid2map(ARR_GRID_IN, V, NX=NUM_NX, NY=NUM_NY, NC_TYPE="INT")
                for n in range(len(map_in)):
                    NCDF4_DATA.variables[V][n] = map_in[n]
            for V1 in ["CENTER"]:
                for V2 in ["LON", "LAT"]:
                    map_in = self.convert_grid2map(ARR_GRID_IN, V1, V2, NX=NUM_NX, NY=NUM_NY, NC_TYPE="FLOAT")
                    for n in range(len(map_in)):
                        NCDF4_DATA.variables["{0:s}_{1:s}".format(V1, V2)][n] = map_in[n]
            
            for V1 in ARR_VAR_STR:
                for V2 in ARR_VAR_ITEM:
                    map_in = self.convert_grid2map(ARR_GRID_IN, V1, V2, NX=NUM_NX, NY=NUM_NY, NT=NUM_NT)
                    for n in range(len(map_in)):
                        NCDF4_DATA.groups[V1].variables[V2][n] = map_in[n]
        NCDF4_DATA.close()

    def export_grid(ARR_GRID_IN, STR_DIR, STR_FILENAME, ARR_VAR_STR=[], STR_GRID_NAME="",\
                    STR_TYPE="netCDF4", IF_PB=False ):
        TIME_NOW = time.gmtime()
        STR_DATE_NOW = "{0:04d}-{1:02d}-{2:02d}".format(TIME_NOW.tm_year, TIME_NOW.tm_mon, TIME_NOW.tm_mday)
        STR_TIME_NOW = "{0:04d}:{1:02d}:{2:02d}".format(TIME_NOW.tm_hour, TIME_NOW.tm_min, TIME_NOW.tm_sec)
        if STR_GRID_NAME == "":
            STR_GRID_NAME = [ k for k,v in locals().items() if v == ARR_GRID_IN][0]
        if STR_TYPE == "netCDF4":
            NCDF4_DATA = Dataset("{0:s}/{1:s}".format(STR_DIR, STR_FILENAME), 'a', format="NETCDF4")
            # CREATE ATTRIBUTEs:
            NCDF4_DATA.description = \
            "The grid information in netCDF4, the grid information will be based on the input grids. "
            NCDF4_DATA.history = "Create on {0:s} at {1:s}".format(STR_DATE_NOW, STR_TIME_NOW)
    
            # CREATE THE DIMENSIONS:
    
            #numGrids = len(ARR_GRID_IN)
            try:
                NCDF4_DATA.createDimension("index", None    )
            except:
                print("Dimension existed")
            # CREATE GROUP for Variables: 
            try:
                group = NCDF4_DATA.createGroup(STR_GRID_NAME)
            except:
                print("Group existed, not creating the group: {0:s}".format(STR_GRID_NAME))
                group = NCDF4_DATA.group[STR_GRID_NAME]
    
            # CREATE BASIC VARIABLES:
            arrVarOut = []
            arrVars   = [ k for k in ARR_GRID_IN[0].keys() ]
            for key in arrVars:
                try:
                    arrSubVar = ARR_GRID_IN[0][key].keys()
                    for subkey in arrSubVar:
                        arrVarOut.append("{0:s}_sub_{1:s}".format(key,subkey))
                except:
                    arrVarOut.append(key)
    
            try:
                for key in arrVarOut:
                    group.createVariable( key, "f8", ("index"))
            except:
                print("Variables existed")
    
            # WRITE IN VARIABLE
            NUM_LEN_GRID = len(ARR_GRID_IN)
            for IND, OBJ in enumerate(ARR_GRID_IN):
                for key in arrVarOut:
                    if len(re.findall("_sub_", key)) == 0:
                        v_in  = group.variables[key]
                        v_in[IND] = OBJ[key]
                    else:
                        v_in  = group.variables[key]
                        key1, key2 = re.split("_sub_", key)
                        v_in[IND] = OBJ[key1][key2]
                        #print(key1, key2, IND, OBJ[key1][key2] )
    
                if IF_PB: TOOLS.progress_bar(IND/(NUM_LEN_GRID-1), STR_DES="WRITING PROGRESS")
            NCDF4_DATA.close()

    def import_grid(STR_DIR, STR_FILENAME, STR_GRID_NAME="",\
                STR_TYPE="netCDF4", IF_PB=False ):
        NC_IN = Dataset("{0:s}/{1:s}".format(STR_DIR, STR_FILENAME), 'r', format="NETCDF4")
        group = NC_IN.groups[STR_GRID_NAME]
        ARR_GRIDS = []
        numGrid_size = 0
        for ind in range(NC_IN.dimensions["index"].size):
            chk   = group.variables["INDEX"][ind].mask
            value = int(group.variables["INDEX"][ind])
            if chk:
                break
            else:
                numGrid_size+=1
                Obj_In = {"INDEX" : value   }
    
                for Key in group.variables:
                    v_in = group.variables[Key]
                    if len(re.split("_sub_",Key)) == 2:
                        key1, key2 = re.split("_sub_",Key)
                        Obj_In[key1] = {key2 : float(v_in[value])}                            
                    elif len(re.split("_sub_",Key)) == 0:
                        Obj_In[Key] = float(v_in[value])
                
                ARR_GRIDS.append(Obj_In)
            TOOLS.progress_bar(ind/NC_IN.dimensions["index"].size)
        return ARR_GRIDS


    def export_reference_map(self, STR_DIR, STR_FILENAME, STR_TYPE="netCDF4", IF_PB=False, IF_PARALLEL=False ):
        TIME_NOW = time.gmtime()
        self.STR_DATE_NOW = "{0:04d}-{1:02d}-{2:02d}".format(TIME_NOW.tm_year, TIME_NOW.tm_mon, TIME_NOW.tm_mday) 
        self.STR_TIME_NOW = "{0:02d}:{1:02d}:{2:02d}".format(TIME_NOW.tm_hour, TIME_NOW.tm_min, TIME_NOW.tm_sec)
        STR_INPUT_FILENAME = "{0:s}/{1:s}".format(STR_DIR, STR_FILENAME)
        if STR_TYPE == "netCDF4":
            IF_FILECHK = os.path.exists(STR_INPUT_FILENAME) 
            if IF_FILECHK:
                NCDF4_DATA = NC.Dataset(STR_INPUT_FILENAME, 'a', format="NETCDF4", parallel=IF_PARALLEL)
                INDEX          = NCDF4_DATA.variables["INDEX"          ]
                INDEX_J        = NCDF4_DATA.variables["INDEX_J"        ]
                INDEX_I        = NCDF4_DATA.variables["INDEX_I"        ]
                CENTER_LON     = NCDF4_DATA.variables["CENTER_LON"     ]
                CENTER_LAT     = NCDF4_DATA.variables["CENTER_LAT"     ]
                INDEX_REF      = NCDF4_DATA.variables["INDEX_REF"      ]
                INDEX_REF_J    = NCDF4_DATA.variables["INDEX_REF_J"    ]
                INDEX_REF_I    = NCDF4_DATA.variables["INDEX_REF_I"    ]
                CENTER_REF_LON = NCDF4_DATA.variables["CENTER_REF_LON" ]
                CENTER_REF_LAT = NCDF4_DATA.variables["CENTER_REF_LAT" ]
            else: 
                NCDF4_DATA = NC.Dataset(STR_INPUT_FILENAME, 'w', format="NETCDF4", parallel=IF_PARALLEL)
                # CREATE ATTRIBUTEs:
                NCDF4_DATA.description = \
                "The netCDF4 version of reference map which contains grid information for resampling"
                NCDF4_DATA.history = "Create on {0:s} at {1:s}".format(self.STR_DATE_NOW, self.STR_TIME_NOW)
                # CREATE DIMENSIONs:
                NCDF4_DATA.createDimension("Y",self.NUM_NY)
                NCDF4_DATA.createDimension("X",self.NUM_NX)
                # CREATE_VARIABLES:
                INDEX          = NCDF4_DATA.createVariable("INDEX",          "i4", ("Y", "X"))
                INDEX_J        = NCDF4_DATA.createVariable("INDEX_J",        "i4", ("Y", "X"))
                INDEX_I        = NCDF4_DATA.createVariable("INDEX_I",        "i4", ("Y", "X"))
                CENTER_LON     = NCDF4_DATA.createVariable("CENTER_LON",     "f8", ("Y", "X"))
                CENTER_LAT     = NCDF4_DATA.createVariable("CENTER_LAT",     "f8", ("Y", "X"))
                INDEX_REF      = NCDF4_DATA.createVariable("INDEX_REF",      "i4", ("Y", "X"))
                INDEX_REF_J    = NCDF4_DATA.createVariable("INDEX_REF_J",    "i4", ("Y", "X"))
                INDEX_REF_I    = NCDF4_DATA.createVariable("INDEX_REF_I",    "i4", ("Y", "X"))
                CENTER_REF_LON = NCDF4_DATA.createVariable("CENTER_REF_LON", "f8", ("Y", "X"))
                CENTER_REF_LAT = NCDF4_DATA.createVariable("CENTER_REF_LAT", "f8", ("Y", "X"))
            NUM_TOTAL_OBJ = len(self.ARR_REFERENCE_MAP)
            NUM_MAX_I     = self.NUM_NX
            for OBJ in self.ARR_REFERENCE_MAP:
                j = OBJ["INDEX_J"]
                i = OBJ["INDEX_I"]
                INDEX[j,i]            = OBJ["INDEX"]
                INDEX_J[j,i]          = OBJ["INDEX_J"]
                INDEX_I[j,i]          = OBJ["INDEX_I"]
                INDEX_REF[j,i]        = OBJ["INDEX_REF"]
                INDEX_REF_J[j,i]      = OBJ["INDEX_REF_J"]
                INDEX_REF_I[j,i]      = OBJ["INDEX_REF_I"]
                CENTER_LON [j,i]      = OBJ["CENTER"]["LON"]
                CENTER_LAT [j,i]      = OBJ["CENTER"]["LAT"]
                CENTER_REF_LON [j,i]  = OBJ["CENTER_REF"]["LON"]
                CENTER_REF_LAT [j,i]  = OBJ["CENTER_REF"]["LAT"]
                if IF_PB: TOOLS.progress_bar((i+j*NUM_MAX_I)/float(NUM_TOTAL_OBJ), STR_DES="Exporting")
            NCDF4_DATA.close()
        
    def import_reference_map(self, STR_DIR, STR_FILENAME, ARR_X_RANGE=[], ARR_Y_RANGE=[], STR_TYPE="netCDF4", IF_PB=False):
        self.ARR_REFERENCE_MAP = []
        self.NUM_MAX_INDEX_RS = 0
        self.NUM_MIN_INDEX_RS = 999
        if len(ARR_X_RANGE) != 0:
            self.I_MIN = ARR_X_RANGE[0]
            self.I_MAX = ARR_X_RANGE[1]
        else:
            self.I_MIN = 0 
            self.I_MAX = self.REFERENCE_MAP_NX

        if len(ARR_Y_RANGE) != 0:
            self.J_MIN = ARR_Y_RANGE[0]
            self.J_MAX = ARR_Y_RANGE[1]
        else:
            self.J_MIN = 0
            self.J_MAX = self.REFERENCE_MAP_NY

        if STR_TYPE == "netCDF4":
            NCDF4_DATA = NC.Dataset("{0:s}/{1:s}".format(STR_DIR, STR_FILENAME), 'r', format="NETCDF4")
            # READ DIMENSIONs:
            self.REFERENCE_MAP_NY = NCDF4_DATA.dimensions["Y"].size
            self.REFERENCE_MAP_NX = NCDF4_DATA.dimensions["X"].size
            
            # CREATE_VARIABLES:
            INDEX          = NCDF4_DATA.variables["INDEX"          ]
            INDEX_J        = NCDF4_DATA.variables["INDEX_J"        ]
            INDEX_I        = NCDF4_DATA.variables["INDEX_I"        ]
            CENTER_LON     = NCDF4_DATA.variables["CENTER_LON"     ]
            CENTER_LAT     = NCDF4_DATA.variables["CENTER_LAT"     ]
            INDEX_REF      = NCDF4_DATA.variables["INDEX_REF"      ]
            INDEX_REF_J    = NCDF4_DATA.variables["INDEX_REF_J"    ]
            INDEX_REF_I    = NCDF4_DATA.variables["INDEX_REF_I"    ]
            CENTER_REF_LON = NCDF4_DATA.variables["CENTER_REF_LON" ]
            CENTER_REF_LAT = NCDF4_DATA.variables["CENTER_REF_LAT" ]
            
            for j in range(self.J_MIN, self.J_MAX):
                for i in range(self.I_MIN, self.I_MAX):
                    OBJ_ELEMENT = {"INDEX"       :                         0 ,\
                                   "INDEX_I"     :                         0 ,\
                                   "INDEX_J"     :                         0 ,\
                                   "CENTER"      :  {"LAT": 0.0, "LON": 0.0} ,\
                                   "INDEX_REF"   :                         0 ,\
                                   "INDEX_REF_I" :                         0 ,\
                                   "INDEX_REF_J" :                         0 ,\
                                   "CENTER_REF"  :  {"LAT": 0.0, "LON": 0.0} }
                    if INDEX         [j][i] != None:
                        OBJ_ELEMENT["INDEX"]             = INDEX         [j][i] 
                        OBJ_ELEMENT["INDEX_J"]           = INDEX_J       [j][i]
                        OBJ_ELEMENT["INDEX_I"]           = INDEX_I       [j][i]
                        OBJ_ELEMENT["INDEX_REF"]         = INDEX_REF     [j][i]
                        OBJ_ELEMENT["INDEX_REF_J"]       = INDEX_REF_J   [j][i]
                        OBJ_ELEMENT["INDEX_REF_I"]       = INDEX_REF_I   [j][i]
                        OBJ_ELEMENT["CENTER"]["LAT"]     = CENTER_LAT    [j][i]
                        OBJ_ELEMENT["CENTER"]["LON"]     = CENTER_LON    [j][i]
                        OBJ_ELEMENT["CENTER_REF"]["LAT"] = CENTER_REF_LAT[j][i]
                        OBJ_ELEMENT["CENTER_REF"]["LON"] = CENTER_REF_LON[j][i]
                    else:
                        OBJ_ELEMENT["INDEX"]             = INDEX         [j][i] 
                        OBJ_ELEMENT["INDEX_I"]           = INDEX_J       [j][i]
                        OBJ_ELEMENT["INDEX_J"]           = INDEX_I       [j][i]
                        OBJ_ELEMENT["INDEX_REF"]         = -999
                        OBJ_ELEMENT["INDEX_REF_J"]       = -999
                        OBJ_ELEMENT["INDEX_REF_I"]       = -999
                        OBJ_ELEMENT["CENTER"]["LAT"]     = CENTER_LAT    [j][i]
                        OBJ_ELEMENT["CENTER"]["LON"]     = CENTER_LON    [j][i]
                        OBJ_ELEMENT["CENTER_REF"]["LAT"] = -999
                        OBJ_ELEMENT["CENTER_REF"]["LON"] = -999
                    self.ARR_REFERENCE_MAP.append(OBJ_ELEMENT)
                    self.NUM_MIN_INDEX_RS = min(self.NUM_MIN_INDEX_RS, INDEX_REF[j][i])
                    self.NUM_MAX_INDEX_RS = max(self.NUM_MAX_INDEX_RS, INDEX_REF[j][i])

                if IF_PB: TOOLS.progress_bar((j - self.J_MIN + 1)/float(self.J_MAX - self.J_MIN), STR_DES="IMPORTING")
            if self.NUM_MIN_INDEX_RS == 0:
                self.NUM_MAX_RS = self.NUM_MAX_INDEX_RS + 1
            NCDF4_DATA.close()    

    def create_resample_map(self, ARR_REFERENCE_MAP=[], ARR_VARIABLES=["Value"], ARR_GRID_IN=[],\
                             IF_PB=False, NUM_NT=0, NUM_NX=0, NUM_NY=0, NUM_NULL=-9999.999):
        if NUM_NT == 0:
            NUM_NT = self.NUM_NT
        if NUM_NX == 0:
            NUM_NX = self.NUM_NX
        if NUM_NY == 0:
            NUM_NY = self.NUM_NY
        if len(ARR_REFERENCE_MAP) == 0:
            self.ARR_RESAMPLE_OUT = []
            self.ARR_RESAMPLE_OUT_PARA = {"EDGE": {"N": 0.0,"S": 0.0,"E": 0.0,"W": 0.0}}
            NUM_END_J = self.NUM_GRIDS_LAT - 1
            NUM_END_I = self.NUM_GRIDS_LON - 1
            ARR_EMPTY = [float("NaN") for n in range(self.NUM_NT)]
            for J in range(self.NUM_GRIDS_LAT):
                for I in range(self.NUM_GRIDS_LON):
                    NUM_IND = I + J * self.NUM_GRIDS_LON
                    self.add_an_geo_element(self.ARR_RESAMPLE_OUT, NUM_INDEX=NUM_IND, NUM_J=J, NUM_I=I, \
                               NUM_NX= self.NUM_GRIDS_LON, NUM_NY= self.NUM_GRIDS_LAT,\
                               ARR_VALUE_STR=ARR_VARIABLES, NUM_NT=NUM_NT) 
            self.ARR_RESAMPLE_MAP_PARA["EDGE"]["N"] = max( self.ARR_LAT[NUM_END_J][0], self.ARR_LAT[NUM_END_J][NUM_END_I] )
            self.ARR_RESAMPLE_MAP_PARA["EDGE"]["S"] = min( self.ARR_LAT[0][0], self.ARR_LAT[0][NUM_END_I] )
            self.ARR_RESAMPLE_MAP_PARA["EDGE"]["W"] = min( self.ARR_LAT[0][0], self.ARR_LAT[NUM_END_J][0] )
            self.ARR_RESAMPLE_MAP_PARA["EDGE"]["E"] = max( self.ARR_LAT[0][NUM_END_I], self.ARR_LAT[NUM_END_J][NUM_END_I] )
            self.NUM_MAX_INDEX_RS = NUM_IND
            
        else:
            if ARR_GRID_IN == []: ARR_GRID_IN = self.ARR_GRID
            self.ARR_RESAMPLE_OUT = [ {} for n in range(NUM_NX * NUM_NY)]
            for IND in range(len(self.ARR_RESAMPLE_OUT)):
                for VAR in ARR_VARIABLES:
                    self.ARR_RESAMPLE_OUT[IND][VAR] = [{"VALUE" : []} for T in range(NUM_NT) ]
            #for IND in range(len(ARR_REFERENCE_MAP)):
            for IND in range(len(ARR_GRID_IN)):
                R_IND = ARR_REFERENCE_MAP[IND]["INDEX_REF"] 
                R_J   = ARR_REFERENCE_MAP[IND]["INDEX_REF_J"]
                R_I   = ARR_REFERENCE_MAP[IND]["INDEX_REF_I"]
                R_IND_FIX = TOOLS.fix_ind(R_IND, R_J, R_I, ARR_XRANGE=self.ARR_RESAMPLE_LIM_X, ARR_YRANGE=self.ARR_RESAMPLE_LIM_Y, NX=NUM_NX, NY=NUM_NY)
                if R_IND != None:
                    for VAR in ARR_VARIABLES:
                        for T in range(NUM_NT):
                            #print("R_IND:{0:d}, T:{1:d}, IND:{2:d} ".format(R_IND, T, IND))
                            NUM_VAL_IN = ARR_GRID_IN[IND][VAR][T]["VALUE"]
                            self.ARR_RESAMPLE_OUT[R_IND][VAR][T]["VALUE"].append(NUM_VAL_IN)
                    self.ARR_RESAMPLE_OUT[R_IND]["INDEX"]         = ARR_REFERENCE_MAP[IND]["INDEX_REF"] 
                    self.ARR_RESAMPLE_OUT[R_IND]["INDEX_J"]       = ARR_REFERENCE_MAP[IND]["INDEX_REF_J"]
                    self.ARR_RESAMPLE_OUT[R_IND]["INDEX_I"]       = ARR_REFERENCE_MAP[IND]["INDEX_REF_I"]
                    self.ARR_RESAMPLE_OUT[R_IND]["CENTER"] = {"LAT": 0.0, "LON": 0.0 }
                    self.ARR_RESAMPLE_OUT[R_IND]["CENTER"]["LAT"] = ARR_REFERENCE_MAP[IND]["CENTER"]["LAT"]
                    self.ARR_RESAMPLE_OUT[R_IND]["CENTER"]["LON"] = ARR_REFERENCE_MAP[IND]["CENTER"]["LON"]
                if IF_PB: TOOLS.progress_bar(TOOLS.cal_loop_progress([IND], [len(ARR_GRID_IN)]), STR_DES="RESAMPLING PROGRESS")
                
    def cal_resample_map(self, ARR_VARIABLES, ARR_GRID_IN=[], NUM_NT=0, IF_PB=False, \
                         DIC_PERCENTILE={ "P05": 0.05, "P10": 0.1, "P25": 0.25, "P75": 0.75, "P90": 0.90, "P95": 0.95}, NUM_NULL=-9999.999):
        if NUM_NT == 0:
            NUM_NT = self.NUM_NT
        NUM_RS_OUT_LEN = len(self.ARR_RESAMPLE_OUT)
        for IND in range(NUM_RS_OUT_LEN):
            for VAR in ARR_VARIABLES:
                for T in range(NUM_NT):
                    ARR_IN          = self.ARR_RESAMPLE_OUT[IND][VAR][T]["VALUE"]
                    if len(ARR_IN) > 0:
                        ARR_IN.sort()
                        NUM_ARR_LEN     = len(ARR_IN)
                        NUM_ARR_MEAN    = sum(ARR_IN) / float(NUM_ARR_LEN)
                        NUM_ARR_S2SUM   = 0
                        if math.fmod(NUM_ARR_LEN,2) == 1:
                            NUM_MPOS = [int((NUM_ARR_LEN-1)/2.0), int((NUM_ARR_LEN-1)/2.0)]
                        else: 
                            NUM_MPOS = [int(NUM_ARR_LEN/2.0)    , int(NUM_ARR_LEN/2.0 -1) ]
                        
                        self.ARR_RESAMPLE_OUT[IND][VAR][T]["MIN"]     = min(ARR_IN)
                        self.ARR_RESAMPLE_OUT[IND][VAR][T]["MAX"]     = max(ARR_IN)
                        self.ARR_RESAMPLE_OUT[IND][VAR][T]["MEAN"]    = NUM_ARR_MEAN 
                        self.ARR_RESAMPLE_OUT[IND][VAR][T]["MEDIAN"]  = ARR_IN[NUM_MPOS[0]] *0.5 + ARR_IN[NUM_MPOS[1]] *0.5
                        
                        for STVA in DIC_PERCENTILE:
                            self.ARR_RESAMPLE_OUT[IND][VAR][T][STVA]  = ARR_IN[ round(NUM_ARR_LEN * DIC_PERCENTILE[STVA])-1]
                        for VAL in ARR_IN:
                            NUM_ARR_S2SUM += (VAL - NUM_ARR_MEAN)**2
                        self.ARR_RESAMPLE_OUT[IND][VAR][T]["STD"]     =  (NUM_ARR_S2SUM / max(1, NUM_ARR_LEN-1))**0.5
            if IF_PB: TOOLS.progress_bar(TOOLS.cal_loop_progress([IND], [NUM_RS_OUT_LEN]), STR_DES="RESAMPLING CALCULATION")            

    def convert_grid2map(self, ARR_GRID_IN, STR_VAR, STR_VAR_TYPE="", NX=0, NY=0, NT=0, IF_PB=False, NC_TYPE=""):

        if NC_TYPE == "INT":
            if NT == 0:
                ARR_OUT = NP.empty([NY, NX], dtype=NP.int8)
            else:
                ARR_OUT = NP.empty([NT, NY, NX], dtype=NP.int8)
        elif NC_TYPE == "FLOAT":
            if NT == 0:
                ARR_OUT = NP.empty([NY, NX], dtype=NP.float64)
            else:
                ARR_OUT = NP.empty([NT, NY, NX], dtype=NP.float64)
        else:
            if NT == 0:
                ARR_OUT = [[ self.NUM_NULL for i in range(NX)] for j in range(NY) ]
            else:
                ARR_OUT = [[[ self.NUM_NULL for i in range(NX)] for j in range(NY) ] for t in range(NT)]

        if STR_VAR_TYPE == "":
            for I, GRID in enumerate(ARR_GRID_IN):
                if GRID["INDEX"] != -999:
                    if NT == 0:
                        #print(GRID["INDEX_J"], GRID["INDEX_I"], GRID[STR_VAR])
                        ARR_OUT[ GRID["INDEX_J"] ][ GRID["INDEX_I"] ] = GRID[STR_VAR]
                    else:
                        for T in range(NT):
                            ARR_OUT[T][ GRID["INDEX_J"] ][ GRID["INDEX_I"] ] = GRID[STR_VAR][T]
                if IF_PB==True: TOOLS.progress_bar(((I+1)/(len(ARR_GRID_IN))))
        else:
            for I, GRID in enumerate(ARR_GRID_IN):
                if GRID["INDEX"] != -999:
                    if NT == 0:
                        ARR_OUT[ GRID["INDEX_J"] ][ GRID["INDEX_I"] ] = GRID[STR_VAR][STR_VAR_TYPE]
                    else:
                        for T in range(NT):
                            ARR_OUT[T][ GRID["INDEX_J"] ][ GRID["INDEX_I"] ] = GRID[STR_VAR][T][STR_VAR_TYPE]
                if IF_PB==True: TOOLS.progress_bar(((I+1)/(len(ARR_GRID_IN))))
        return ARR_OUT

    def mask_grid(self, ARR_GRID_IN, STR_VAR, STR_VAR_TYPE, NUM_NT=0, STR_MASK="MASK",\
              ARR_NUM_DTM=[0,1,2], ARR_NUM_DTM_RANGE=[0,1]):
        if NUM_NT == 0:
            NUM_NT= self.NUM_NT
        for IND, GRID in enumerate(ARR_GRID_IN):
            for T in range(NUM_NT):
                NUM_DTM = GEO_TOOLS.mask_dtm(GRID[STR_VAR][T][STR_VAR_TYPE], ARR_NUM_DTM=ARR_NUM_DTM, ARR_NUM_DTM_RANGE=ARR_NUM_DTM_RANGE)
                ARR_GRID_IN[IND][STR_VAR][T][STR_MASK] = NUM_DTM

class MATH_TOOLS:
    """ Some math tools that help us to calculate. 
        gau_kde: kernel density estimator by Gaussian Function
        standard_dev: The Standard deviation
    """
    def GaussJordanEli(arr_in):
        num_ydim = len(arr_in)
        num_xdim = len(arr_in[0])
        arr_out = arr_in
        
        if num_ydim -num_xdim == 0 or num_xdim - num_ydim == 1:
            arr_i   = NP.array([[0.0 for j in range(num_ydim)] for i in range(num_ydim)])
            for ny in range(num_ydim):
                arr_i[ny][ny] = 1.0
            #print(arr_i)
            for nx in range(num_xdim):                        
                for ny in range(nx+1, num_ydim):
                    arr_i  [ny] = arr_i  [ny] - arr_i  [nx] * arr_out[ny][nx] / float(arr_out[nx][nx])
                    arr_out[ny] = arr_out[ny] - arr_out[nx] * arr_out[ny][nx] / float(arr_out[nx][nx])
            if num_xdim - num_ydim == 1:
                for nx in range(num_xdim-1,-1,-1):                        
                    for ny in range(num_ydim-1,nx, -1):
                        print(nx,ny)
                        arr_i  [nx] = arr_i  [nx] - arr_i  [ny] * arr_out[nx][ny] / float(arr_out[ny][ny])
                        arr_out[nx] = arr_out[nx] - arr_out[ny] * arr_out[nx][ny] / float(arr_out[ny][ny])
            else:
                for nx in range(num_xdim,-1,-1):                        
                    for ny in range(num_ydim-1, nx, -1):
                        print(nx,ny)
                        arr_i  [nx] = arr_i  [nx] - arr_i  [ny] * arr_out[nx][ny] / float(arr_out[ny][ny])
                        arr_out[nx] = arr_out[nx] - arr_out[ny] * arr_out[nx][ny] / float(arr_out[ny][ny])
                    
            if num_xdim - num_ydim == 1:
                arr_sol = [0.0 for n in range(num_ydim)]
                for ny in range(num_ydim):
                    arr_sol[ny] = arr_out[ny][num_xdim-1]/arr_out[ny][ny]
                return arr_out, arr_i, arr_sol
            else:
                return arr_out, arr_i
        else:
            print("Y dim: {0:d}, X dim: {1:d}: can not apply Gaussian-Jordan".format(num_ydim, num_xdim))
            return [0]
        
    def finding_XM_LSM(arr_in1, arr_in2, m=2):
        # Finding the by least square method
        arr_out=[[0.0 for i in range(m+2)] for j in range(m+1)]
        arr_x_power_m  = [0.0 for i in range(m+m+1)]
        arr_xy_power_m = [0.0 for i in range(m+1)]
        for n in range(len(arr_x_power_m)):
            for x in range(len(arr_in1)):
                arr_x_power_m[n] += arr_in1[x] ** n
        
        for n in range(len(arr_xy_power_m)):
            for x in range(len(arr_in1)):
                arr_xy_power_m[n] += arr_in1[x] ** n * arr_in2[x]
        
        for j in range(m+1):
            for i in range(j,j+m+1):
                arr_out[j][i-j] = arr_x_power_m[i]
            arr_out[j][m+1]   = arr_xy_power_m[j]
        return arr_out

    def cal_modelperform (arr_obs , arr_sim , num_empty=-999.999):
        # Based on Vazquez et al. 2002 (Hydrol. Process.)
        num_arr = len(arr_obs)
        num_n_total = num_arr
        num_sum = 0
        num_obs_sum = 0
    
        for n in range( num_arr ):
            if  math.isnan(arr_obs[n]) or arr_obs[n] == num_empty:
                num_n_total += -1
            else: 
                num_sum = num_sum + ( arr_sim[n] - arr_obs[n] )    ** 2
                num_obs_sum = num_obs_sum +    arr_obs[n]
        if num_n_total == 0 or num_obs_sum == 0:
            RRMSE   = -999.999
            RMSE    = -999.999
            obs_avg = -999.999
        else:
            RRMSE = ( num_sum / num_n_total ) ** 0.5 *    ( num_n_total / num_obs_sum )
            RMSE  = ( num_sum / num_n_total ) ** 0.5 
            obs_avg = num_obs_sum / num_n_total
    
        num_n_total = num_arr
        oo_sum = 0
        po_sum = 0
        for nn in range( num_arr ):
            if  math.isnan(arr_obs[nn]) or arr_obs[nn] == num_empty:
                num_n_total = num_n_total - 1
            else:
                oo_sum = oo_sum + ( arr_obs[nn] - obs_avg )     ** 2
                po_sum = po_sum + ( arr_sim[nn] - arr_obs[nn] ) ** 2
    
        if num_n_total == 0 or oo_sum * po_sum == 0:
            EF = -999.999
            CD = -999.999
        else:
            EF = ( oo_sum - po_sum ) / oo_sum
            CD = oo_sum / po_sum
        return RRMSE,EF,CD,RMSE, num_arr

    def cal_kappa(ARR_IN, NUM_n=0, NUM_N=0, NUM_k=0):
        """ Fleiss' kappa
            Mustt input with ARR_IN in the following format:
            ARR_IN = [ [ NUM for k in range(catalogue)] for N in range(Subjects)]
            Additional parameters: NUM_n is the number of raters (e.g. sim and obs results)
            Additional parameters: NUM_N is the number of subjects (e.g the outputs
            Additional parameters: NUM_k is the number of catalogue (e.g. results )
        """
        if NUM_N == 0:
            NUM_N = len(ARR_IN)
        if NUM_n == 0:
            NUM_n = sum(ARR_IN[0])
        if NUM_k == 0:
            NUM_k = len(ARR_IN[0])
        ARR_p_out = [ 0 for n in range(NUM_k)]
        ARR_P_OUT = [ 0 for n in range(NUM_N)]

        for N in range(NUM_N):
            for k in range(NUM_k):
                ARR_p_out[k] += ARR_IN[N][k]
                ARR_P_OUT[N] += ARR_IN[N][k] ** 2        
            ARR_P_OUT[N] -= NUM_n
            ARR_P_OUT[N]  = ARR_P_OUT[N] * (1./(NUM_n *(NUM_n - 1)))
        for k in range(NUM_k):
            ARR_p_out[k] = ARR_p_out[k] / (NUM_N * NUM_n)
        NUM_P_BAR = 0
        for N in range(NUM_N):
            NUM_P_BAR += ARR_P_OUT[N]
        NUM_P_BAR = NUM_P_BAR / float(NUM_N)

        NUM_p_bar = 0
        for k in ARR_p_out:
            NUM_p_bar += k **2

        return (NUM_P_BAR - NUM_p_bar) / (1 - NUM_p_bar)

    def gau_kde(ARR_IN_X, ARR_IN_I, NUM_BW=0.1  ):
        NUM_SUM = 0.
        NUM_LENG = len(ARR_IN_X)
        ARR_OUT  = [ 0. for n in range(NUM_LENG)]
        for IND_J, J in enumerate(ARR_IN_X):
            NUM_SUM = 0.0
            for I in ARR_IN_I:
                NUM_SUM += 1 / (2 * math.pi)**0.5 * math.e ** (-0.5 * ((J-I)/NUM_BW) ** 2 ) 
            ARR_OUT[IND_J] = NUM_SUM / len(ARR_IN_I) / NUM_BW
        return ARR_OUT 
        
    def standard_dev(ARR_IN):
        NUM_SUM = sum(ARR_IN)
        NUM_N   = len(ARR_IN)
        NUM_MEAN = 1.0*NUM_SUM/NUM_N
        NUM_SUM2 = 0.0
        for N in ARR_IN:
            if not math.isnan(N):
                NUM_SUM2 += (N-NUM_MEAN)**2
            else:
                NUM_N += -1
        return (NUM_SUM2 / (NUM_N-1)) ** 0.5
    
    def h_esti(ARR_IN):
        #A rule-of-thumb bandwidth estimator
        NUM_SIGMA = standard_dev(ARR_IN)
        NUM_N     = len(ARR_IN)
        return ((4 * NUM_SIGMA ** 5) / (3*NUM_N)  ) ** 0.2
    
    def data2array(ARR_IN, STR_IN="MEAN"):
        NUM_J = len(ARR_IN)
        NUM_I = len(ARR_IN[0])
        ARR_OUT = [[ 0.0 for i in range(NUM_I)] for j in range(NUM_J) ]
        for j in range(NUM_J):
            for i in range(NUM_I):
                ARR_OUT[j][i] = ARR_IN[j][i][STR_IN]
        return ARR_OUT

    def reshape2d(ARR_IN, NUM_NULL = 0.0):
        ARR_OUT=[]
        for A in ARR_IN:
            for B in A:
                if not math.isnan(B):
                    ARR_OUT.append(B)
        return ARR_OUT

    def NormalVector( V1, V2):
        return [(V1[1]*V2[2] - V1[2]*V2[1]), (V1[2]*V2[0] - V1[0]*V2[2]),(V1[0]*V2[1] - V1[1]*V2[0])]  
    
    def NVtoPlane( P0, P1, P2):
        """Input of P should be 3-dimensionals"""
        V1 = [(P1[0]-P0[0]),(P1[1]-P0[1]),(P1[2]-P0[2])]
        V2 = [(P2[0]-P0[0]),(P2[1]-P0[1]),(P2[2]-P0[2])]
        ARR_NV = MATH_TOOLS.NormalVector(V1, V2)
        D = ARR_NV[0] * P0[0] +  ARR_NV[1] * P0[1] + ARR_NV[2] * P0[2]
        return ARR_NV[0],ARR_NV[1],ARR_NV[2],D
        
    def FindZatP3( P0, P1, P2, P3):
        """ input of P: (X,Y,Z); but P3 is (X,Y) only """
        A,B,C,D = MATH_TOOLS.NVtoPlane(P0, P1, P2)
        return (D-A*P3[0] - B*P3[1])/float(C)

    def c2p(x,y):
        if y != 0.0:
            return math.fmod(math.atan(float(x/y)) * 180./math.pi + (x < 0) * 360  + (y < 0)*180, 360)
        else:
            return (x < 0) * 270 + (x > 0) * 90

    def GrtCirDist(lon1,lat1,lon2,lat2,R=6.37122E6):
        a2r = lambda x: x/180. * math.pi

        term1      = math.sin(a2r(lat1)) * math.sin(a2r(lat2))
        term2      = math.cos(a2r(lat1)) * math.cos(a2r(lat2))
        term3      = math.cos(a2r( lon2 - lon1  ))

        deltaSigma = math.acos( term1 + term2 * term3 ) * 180./math.pi
        numDist    = deltaSigma/360. * 2 * R * math.pi

        return { "deltaSigma" : deltaSigma, 
                 "dist"       : numDist     }

    def NearestNeighbor(arr_x, arr_y, target_x, target_y):
        """ This algorithm is used to find the nearest point
            as 1-NN method
            input: ARR_X/ARR_Y, the [j,i] array for x/y coordinates
                   target_x, target_y, the [j,i] array for target [x,y]
        """
        lenNY       = len(arr_x)
        lenNX       = len(arr_x[0])
        numDist_chk = 9.99E20

        for j in range(lenNY - 1):
            for i in range(lenNX - 1):
                chk1    = ( arr_x[j][i] - target_x ) ** 2
                chk3    = ( arr_y[j][i] - target_y ) ** 2
                numDist = (chk1 + chk3 ) ** 0.5
                if min(numDist_chk, numDist) == numDist:
                    numDist_chk = numDist
                    found_i = i
                    found_j = j
         
        return found_i, found_j    

    def NearestNeighbor1D(arr_x, target_v ):
        """ This algorithm is used to find the nearest point
            as 1-NN method
            input: ARR_X/ARR_Y, the [j,i] array for x/y coordinates
                   target_x, target_y, the [j,i] array for target [x,y]
                   K is the number of K
            Temperally solution for WRF and XX YY mesh. 
        """
        numOut = 9E12
        numIndOut = 0
        for ind, v_tmp in enumerate(arr_x):
            chkTmp = (( v_tmp - target_v) ** 2 ) ** 0.5
            if min( chkTmp, numOut) == chkTmp:
                numIndOut  = ind
                numOut     = chkTmp
        return numIndOut, numOut

    def KNearestNeighbor(arr_x, arr_y, target_x, target_y, numK, ):
        """ This algorithm is used to find the nearest point
            as 1-NN method
            input: ARR_X/ARR_Y, the [j,i] array for x/y coordinates
                   target_x, target_y, the [j,i] array for target [x,y]
                   K is the number of K
            Temperally solution for WRF and XX YY mesh. 
        """
        lenNY       = len(arr_x)
        lenNX       = len(arr_x[0])
        arrDistchk = [ [ 9.99E20, 0, 0 ] for n in range(numK) ]
        arrNumChk  = [ 9.99E20 for n in range(numK) ]
        for j in range(lenNY - 1):
            for i in range(lenNX - 1):
                chk1    = ( arr_x[j][i] - target_x ) ** 2
                chk3    = ( arr_y[j][i] - target_y ) ** 2
                numDist = (chk1 + chk3 ) ** 0.5
                for chk in arrDistChk:
                    if min( chk[0], numDist) == numDist:
                        arrDistChk = []
                        found_i = i
                        found_j = j
        return arrDistChk 

    def Find_Xest(arrIn, dim_x, dim_y, num_crop=20, numMax = 9E20):
        numMinOut    =  numMax
        numMaxOut    = -1 * numMax
        Min_i       = 0
        Min_j       = 0
        Max_i       = 0
        Max_j       = 0
        MinV        = 0
        MaxV        = 0
        for j in range(num_crop, dim_y-num_crop):
            for i in range(num_crop, dim_x-num_crop):
                numIn      =  arrIn[j][i]
                numMinOut =  min( numMinOut, numIn )
                numMaxOut =  max( numMaxOut, numIn )
                if numIn == numMinOut:
                    MinV   = numIn
                    Min_i = i
                    Min_j = j
                if numIn == numMaxOut:
                    MaxV   = numIn
                    Max_i = i
                    Max_j = j
        dicOut = {"MinV": MinV, "MinJ": Min_j, "MinI": Min_i,\
                  "MaxV": MaxV, "MaxJ": Max_j, "MaxI": Max_i}
        return dicOut



class TOOLS:
    """ TOOLS is contains:
        timestamp
        fix_ind
        progress_bar
        cal_progrss
    """
    ARR_HOY      = [0, 744, 1416, 2160, 2880, 3624, 4344, 5088, 5832, 6552, 7296, 8016, 8760] 
    ARR_HOY_LEAP = [0, 744, 1440, 2184, 2904, 3648, 4368, 5112, 5856, 6576, 7320, 8040, 8784]

    def NNARR(ARR_IN, IF_PAIRING=False):
        "Clean the NaN value in the array"

        if IF_PAIRING:
            ARR_SIZE = len(ARR_IN)
            ARR_OUT  = [ [] for N in range(ARR_SIZE)]
            for ind_n, N in enumerate(ARR_IN[0]):
                IF_NAN = False
                for ind_a in range(ARR_SIZE):
                    if math.isnan(ARR_IN[ind_a][ind_n]):
                        IF_NAN = True
                        break
                if not IF_NAN:
                    for ind_a in range(ARR_SIZE):
                        ARR_OUT[ind_a].append(ARR_IN[ind_a][ind_n])
        else: 
            ARR_OUT  = [ ]
            for N in ARR_IN:
                if not math.isnan(N):
                    ARR_OUT.append(N)
        return ARR_OUT

    def DATETIME2HOY(ARR_TIME, ARR_HOY_IN=[]):
        if math.fmod(ARR_TIME[0], 4) == 0 and len(ARR_HOY_IN) == 0:
            ARR_HOY_IN = [0, 744, 1440, 2184, 2904, 3648, 4368, 5112, 5856, 6576, 7320, 8040, 8784]
        elif math.fmod(ARR_TIME[0], 4) != 0 and len(ARR_HOY_IN) == 0:
            ARR_HOY_IN = [0, 744, 1416, 2160, 2880, 3624, 4344, 5088, 5832, 6552, 7296, 8016, 8760]
        else:
            ARR_HOY_IN = ARR_HOY_IN 
        return ARR_HOY_IN[ARR_TIME[1]-1] + (ARR_TIME[2]-1)*24 + ARR_TIME[3]

    def timestamp(STR_IN=""):
        str_out = "{0:04d}-{1:02d}-{2:02d}_{3:02d}:{4:02d}:{5:02d} {6:s}".format(time.gmtime().tm_year, time.gmtime().tm_mon, time.gmtime().tm_mday,\
        time.gmtime().tm_hour, time.gmtime().tm_min, time.gmtime().tm_sec, STR_IN) 
        print(str_out)
        return str_out

    def fix_ind(IND_IN, IND_J, IND_I, ARR_XRANGE=[], ARR_YRANGE=[], NX=0, NY=0):
        NUM_DY   = ARR_YRANGE[0]
        NUM_NX_F = ARR_XRANGE[0]
        NUM_NX_R = NX - (ARR_XRANGE[1]+1)
        if IND_J == ARR_YRANGE[0]:
            IND_OUT  = IND_IN - NUM_DY * NX - NUM_NX_F
        else:
            IND_OUT  = IND_IN - NUM_DY * NX - NUM_NX_F * (IND_J - NUM_DY +1) - NUM_NX_R * (IND_J - NUM_DY)
        return IND_OUT   

    def progress_bar(NUM_PROGRESS, NUM_PROGRESS_BIN=0.05, STR_SYS_SYMBOL="=", STR_DES="Progress"):
        NUM_SYM = int(NUM_PROGRESS / NUM_PROGRESS_BIN)
        sys.stdout.write('\r')
        sys.stdout.write('[{0:20s}] {1:4.2f}% {2:s}'.format(STR_SYS_SYMBOL*NUM_SYM, NUM_PROGRESS*100, STR_DES))
        sys.stdout.flush()

    def clean_arr(ARR_IN, CRITERIA=1):
        ARR_OUT=[]
        for i,n in enumerate(ARR_IN):
            if len(n)> CRITERIA:
                ARR_OUT.append(n)
        return ARR_OUT

    def cal_loop_progress(ARR_INDEX, ARR_INDEX_MAX, NUM_CUM_MAX=1, NUM_CUM_IND=1, NUM_TOTAL_MAX=1):
        """ Please list from smallest to largest, i.e.: x->y->z """
        if len(ARR_INDEX) == len(ARR_INDEX_MAX):
            for i, i_index in enumerate(ARR_INDEX):
                NUM_IND_PER = (i_index+1)/float(ARR_INDEX_MAX[i])
                NUM_TOTAL_MAX = NUM_TOTAL_MAX * ARR_INDEX_MAX[i]
                if i >0: NUM_CUM_MAX = NUM_CUM_MAX * ARR_INDEX_MAX[i-1]
                NUM_CUM_IND = NUM_CUM_IND + NUM_CUM_MAX * i_index
            return NUM_CUM_IND / float(NUM_TOTAL_MAX)
        else:
            print("Wrong dimenstion for in put ARR_INDEX ({0:d}) and ARR_INDEX_MAX ({1:d})".format(len(ARR_INDEX), len(ARR_INDEX_MAX)))

    def run_time_cal(ARR_TIME_IN, ARR_BASE_IN= [0,0,0,0,0,0 ], IF_LEAP=False, IF_UNI_MON=False):
        if IF_LEAP == True:
            ARR_DAY_LIM = [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]
        else:
            ARR_DAY_LIM = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]
        if IF_UNI_MON:
            DIFF_DAYS    = (ARR_TIME_IN[0] - ARR_BASE_IN[0] ) * 30 * 12 \
                         + (ARR_TIME_IN[1] - ARR_BASE_IN[1] ) * 30      \
                         + (ARR_TIME_IN[2] - ARR_BASE_IN[2] )
        else:
            DIFF_DAYS    = (ARR_TIME_IN[0] - ARR_BASE_IN[0] ) * sum(ARR_DAY_LIM)            \
                         + (ARR_DAY_LIM[ ARR_TIME_IN[1]] - ARR_DAY_LIM [ ARR_BASE_IN[1] ] ) \
                         + (ARR_TIME_IN[2] - ARR_BASE_IN[2] )
        DIFF_HOURS   = ARR_TIME_IN[3] - ARR_BASE_IN[3]
        DIFF_MINUTES = ARR_TIME_IN[4] - ARR_BASE_IN[4]
        DIFF_SECONDS = ARR_TIME_IN[5] - ARR_BASE_IN[5]
        DAYS         = DIFF_DAYS + DIFF_HOURS/24. + DIFF_MINUTES/(60.*24.) + DIFF_SECONDS/(60.*60.*24.)
        HOURS        = DIFF_DAYS * 24 + DIFF_HOURS + DIFF_MINUTES/60. + DIFF_SECONDS/(60.*60.)
        MINUTES      = DIFF_DAYS * 24 * 60 + DIFF_HOURS * 60 + DIFF_MINUTES + DIFF_SECONDS/(60.)
        SECONDS      = DIFF_DAYS * 24 * 3600 + DIFF_HOURS * 3600 + DIFF_MINUTES * 60 + DIFF_SECONDS
        return {"DAYS": DAYS, "HOURS": HOURS, "MINUTES": MINUTES, "SECONDS": SECONDS }


    def calendar_cal(ARR_START_TIME, ARR_INTERVAL, ARR_END_TIME_IN=[0, 0, 0, 0, 0, 0.0], IF_LEAP=False):
        ARR_END_TIME  = [ 0,0,0,0,0,0.0]
        ARR_DATETIME  = ["SECOND", "MINUTE", "HOUR","DAY", "MON", "YEAR"]
        NUM_ARR_DATETIME = len(ARR_DATETIME)
        IF_FERTIG = False
        ARR_FERTIG = [0,0,0,0,0,0]
        DIC_TIME_LIM = \
        {"YEAR"  : {"START": 0 , "LIMIT": 9999 },\
         "MON"   : {"START": 1 , "LIMIT": 12 },\
         "DAY"   : {"START": 1 , "LIMIT": 31 },\
         "HOUR"  : {"START": 0 , "LIMIT": 23 },\
         "MINUTE": {"START": 0 , "LIMIT": 59 },\
         "SECOND": {"START": 0 , "LIMIT": 59 },\
        }
        for I, T in enumerate(ARR_START_TIME):
            ARR_END_TIME[I] = T + ARR_INTERVAL[I]
        while IF_FERTIG == False:
            if math.fmod(ARR_END_TIME[0],4) == 0: IF_LEAP=True
            if IF_LEAP:
                ARR_DAY_LIM = [0,31,29,31,30,31,30,31,31,30,31,30,31]    
            else:
                ARR_DAY_LIM = [0,31,28,31,30,31,30,31,31,30,31,30,31]
            for I, ITEM in enumerate(ARR_DATETIME):
                NUM_ARR_POS = NUM_ARR_DATETIME-I-1
                if ITEM == "DAY":
                    if ARR_END_TIME[NUM_ARR_POS] > ARR_DAY_LIM[ARR_END_TIME[1]]:
                        ARR_END_TIME[NUM_ARR_POS] = ARR_END_TIME[NUM_ARR_POS] - ARR_DAY_LIM[ARR_END_TIME[1]] 
                        ARR_END_TIME[NUM_ARR_POS - 1] += 1
                else:
                    if ARR_END_TIME[NUM_ARR_POS] > DIC_TIME_LIM[ITEM]["LIMIT"]:
                        ARR_END_TIME[NUM_ARR_POS - 1] += 1
                        ARR_END_TIME[NUM_ARR_POS] = ARR_END_TIME[NUM_ARR_POS] - DIC_TIME_LIM[ITEM]["LIMIT"] - 1
            for I, ITEM in enumerate(ARR_DATETIME):
                NUM_ARR_POS = NUM_ARR_DATETIME-I-1
                if ITEM == "DAY":
                    if ARR_END_TIME[NUM_ARR_POS] <= ARR_DAY_LIM[ARR_END_TIME[1]]: ARR_FERTIG[NUM_ARR_POS] = 1
                else: 
                    if ARR_END_TIME[NUM_ARR_POS] <= DIC_TIME_LIM[ITEM]["LIMIT"]:  ARR_FERTIG[NUM_ARR_POS] = 1
                if sum(ARR_FERTIG) == 6: IF_FERTIG = True
        return ARR_END_TIME

    def DataIntersec(arrIn_3D, x = None, y = None, z = None):
        if x == None and y == None and z == None:
            print("You forgot to set the plane for intersection")
        else:
            dim_x = len(arrIn_3D[0][0])
            dim_y = len(arrIn_3D[0])
            dim_z = len(arrIn_3D)
            if x != None:
                print("Dealing X plane intersection")
                arrOut = [ [ 0.0 for i in range(dim_y) ] for j in range(dim_z) ]
                for k in range(dim_z):
                    for j in range(dim_y):
                        arrOut[k][j] = arrIn_3D[k][j][x]
         
            elif y != None :
                print("Dealing Y plane intersection")
                arrOut = [ [ 0.0 for i in range(dim_x) ] for j in range(dim_z) ]
                for k in range(dim_z):
                    for i in range(dim_x):
                        arrOut[k][i] = arrIn_3D[k][y][i]
            elif z != None:
                print("Dealing Z plane intersection")
                arrOut = [ [ 0.0 for i in range(dim_x) ] for j in range(dim_y) ]
                for j in range(dim_y):
                    for i in range(dim_x):
                        arrOut[j][i] = arrIn_3D[z][j][i]
            return arrOut, dim_x, dim_y, dim_z
    
    def crop_array_2d(arrIn, x_start, x_end, y_start, y_end):
        arrOut = [ [ 0.0 for i in range(x_end - x_start) ] for j in range(y_end - y_start) ]
        for j in range(y_start, y_end):
            for i in range(x_start, x_end):
                arrOut[j-y_start][i-x_start] = arrIn[j][i]
        return arrOut


class MPI_TOOLS:
    def __init__(self, \
                 NUM_NX_END=1, NUM_NY_END=1, NUM_NX_START=0, NUM_NY_START=0, NUM_NX_CORES=0 ,\
                 NUM_NX_TOTAL=1, NUM_NY_TOTAL=1 ):

        """ END number follow the python philisophy: End number is not included in the list """
        try:
            from mpi4py import MPI
        except:
            print("Can not find the mpi4py")
            sys.exit()
 
        self.COMM         = MPI.COMM_WORLD
        self.NUM_MPI_SIZE = self.COMM.Get_size()
        self.NUM_MPI_RANK = self.COMM.Get_rank()

        self.NUM_NX_START = NUM_NX_START
        self.NUM_NY_START = NUM_NY_START
        self.NUM_NX_SIZE  = NUM_NX_END - NUM_NX_START
        self.NUM_NY_SIZE  = NUM_NY_END - NUM_NY_START
        if NUM_NX_CORES == 0:
            self.NUM_NX_CORES = self.NUM_MPI_SIZE
        else:
            self.NUM_NX_CORES = int(NUM_NX_CORES)
        self.NUM_NY_CORES = max(1, int(self.NUM_MPI_SIZE / self.NUM_NX_CORES))
        self.ARR_RANK_DESIGN = [ {} for n in range(self.NUM_MPI_SIZE)]

    def UPDATE_NX_CORES(self, NUM_NX_CORES):
        self.NUM_NX_CORES = int(NUM_NX_CORES)
        self.NUM_NY_CORES = max(1, int(self.NUM_MPI_SIZE / self.NUM_NX_CORES))

    def UPDATE_NY_CORES(self):
        self.NUM_NY_CORES = max(1, int(self.NUM_MPI_SIZE / self.NUM_NX_CORES))

    def CPU_GEOMETRY_2D(self):

        NUM_NX_REMAIN = self.NUM_NX_SIZE  %  self.NUM_NX_CORES
        NUM_NY_REMAIN = self.NUM_NY_SIZE  %  self.NUM_NY_CORES
    
        NUM_NX_DIFF   = int((self.NUM_NX_SIZE - NUM_NX_REMAIN) / self.NUM_NX_CORES  )
        NUM_NY_DIFF   = int((self.NUM_NY_SIZE - NUM_NY_REMAIN) / self.NUM_NY_CORES  )
   
        NUM_NY_DIFF_P1 = NUM_NY_DIFF + 1 
        NUM_NX_DIFF_P1 = NUM_NX_DIFF + 1

        IND_RANK = 0
        ARR_RANK_DESIGN = [ 0 for n in range(self.NUM_MPI_SIZE)]
        for ny in range(self.NUM_NY_CORES):
            for nx in range(self.NUM_NX_CORES):
                NUM_RANK = ny * self.NUM_NX_CORES + nx
                DIC_IN   = {"INDEX_IN": NUM_RANK, "NX_START": 0, "NY_START": 0, "NX_END": 0, "NY_END": 0  }
                if ny < NUM_NY_REMAIN: 
                    DIC_IN["NY_START"] = (ny + 0) * NUM_NY_DIFF_P1 + self.NUM_NY_START
                    DIC_IN["NY_END"  ] = (ny + 1) * NUM_NY_DIFF_P1 + self.NUM_NY_START
                else:
                    DIC_IN["NY_START"] = (ny - NUM_NY_REMAIN + 0) * NUM_NY_DIFF + NUM_NY_REMAIN * NUM_NY_DIFF_P1 + self.NUM_NY_START
                    DIC_IN["NY_END"  ] = (ny - NUM_NY_REMAIN + 1) * NUM_NY_DIFF + NUM_NY_REMAIN * NUM_NY_DIFF_P1 + self.NUM_NY_START

                if nx < NUM_NX_REMAIN:
                    DIC_IN["NX_START"] = (nx + 0) * NUM_NX_DIFF_P1 + self.NUM_NX_START
                    DIC_IN["NX_END"  ] = (nx + 1) * NUM_NX_DIFF_P1 + self.NUM_NX_START
                else:
                    DIC_IN["NX_START"] = (nx - NUM_NX_REMAIN + 0) * NUM_NX_DIFF + NUM_NX_REMAIN * NUM_NX_DIFF_P1 + self.NUM_NX_START
                    DIC_IN["NX_END"  ] = (nx - NUM_NX_REMAIN + 1) * NUM_NX_DIFF + NUM_NX_REMAIN * NUM_NX_DIFF_P1 + self.NUM_NX_START

                ARR_RANK_DESIGN[NUM_RANK] = DIC_IN
        self.ARR_RANK_DESIGN = ARR_RANK_DESIGN
        return ARR_RANK_DESIGN

    def CPU_MAP(self ):
        ARR_CPU_MAP = [ [ NP.nan for i in range(self.NUM_NX_TOTAL)] for j in range(self.NUM_NY_TOTAL) ]
        for RANK in range(len(ARR_RANK_DESIGN)):
            print("DEAL WITH {0:d} {1:d}".format(RANK, ARR_RANK_DESIGN[RANK]["INDEX_IN"] ))
            for jj in range(ARR_RANK_DESIGN[RANK]["NY_START"], ARR_RANK_DESIGN[RANK]["NY_END"]):
                for ii in range(ARR_RANK_DESIGN[RANK]["NX_START"], ARR_RANK_DESIGN[RANK]["NX_END"]):
                    ARR_CPU_MAP[jj][ii] = ARR_RANK_DESIGN[RANK]["INDEX_IN"]  
        return MAP_CPU

    def GATHER_ARR(self, ARR_IN, NUM_GATHERING_RANK=0, ARR_RANK_DESIGN=[]):
        if ARR_RANK_DESIGN == []:
            ARR_RANK_DESIGN = self.ARR_RANK_DESIGN
        if self.NUM_MPI_RANK == NUM_GATHERING_RANK:
            ARR_TMP = [0.0 for size in range(self.NUM_MPI_SIZE)]
        else: 
            ARR_TMP = None
        ARR_TMP = self.COMM.gather(ARR_IN, root=NUM_GATHERING_RANK)
        if self.NUM_MPI_RANK == NUM_GATHERING_RANK:
            ARR_OUT = [[0.0 for i in range(self.NUM_NX_SIZE)] for j in range(self.NUM_NY_SIZE)  ]
            for ind, rank in enumerate(ARR_RANK_DESIGN):
                for j in range(rank["NY_START"], rank["NY_END"]):
                    for i in range(rank["NX_START"], rank["NX_END"]):
                        ARR_OUT[j][i] = ARR_TMP[ind][j][i]
            return ARR_OUT
        else:
            return 0

    def MPI_MESSAGE(self, STR_TEXT=""):
        TIME_NOW = time.gmtime()
        print("MPI RANK: {0:5d} @ {1:02d}:{2:02d}:{3:02d} # {4:s}"\
              .format(self.NUM_MPI_RANK, TIME_NOW.tm_hour, TIME_NOW.tm_min, TIME_NOW.tm_sec, STR_TEXT ))
    
class GEO_TOOLS:
    def __init__(self):
        STR_NCDF4PY = NC.__version__
        print("Using netCDF4 for Python, Version: {0:s}".format(STR_NCDF4PY))

    def mask_dtm(self, NUM, ARR_DTM=[0,1,2], ARR_DTM_RANGE=[0,1], ARR_DTM_STR=["OUT","IN","OUT"]):
        """ The determination algorithm is : x-1 < NUM <= x   """

        for i, n in enumerate(ARR_DTM):
            if i == 0:
                if NUM <= ARR_DTM_RANGE[i]: NUM_OUT = n 
            elif i == len(ARR_DTM_RANGE):
                if NUM > ARR_DTM_RANGE[i-1]: NUM_OUT = n
            else: 
                if NUM > ARR_DTM_RANGE[i-1]  and NUM <= ARR_DTM_RANGE[i]: NUM_OUT = n
        return NUM_OUT 
    
    def mask_array(self, ARR_IN, ARR_MASK_OUT=[], ARR_DTM=[0,1,2], ARR_DTM_RANGE=[0,1], ARR_DTM_STR=["OUT","IN","OUT"], IF_2D=False):
        if IF_2D:
            NUM_NX = len(ARR_IN[0])
            NUM_NY = len(ARR_IN)
            ARR_OUT    = [ [ self.NUM_NULL for i in range(NUM_NX)] for j in range(NUM_NY) ]
            for J in range(NUM_NY):
                for I in range(NUM_NY):
                    ARR_OUT[J][I] = self.mask_dtm(ARR_IN[J][I], ARR_NUM_DTM=ARR_NUM_DTM, ARR_NUM_DTM_RANGE=ARR_NUM_DTM_RANGE, ARR_STR_DTM=ARR_STR_DTM)
        else:
            NUM_NX = len(ARR_IN)
            ARR_OUT = [0 for n in range(NUM_NX)]
            for N in range(NUM_NX):
                ARR_OUT[N] = self.mask_dtm(ARR_IN[N], ARR_NUM_DTM=ARR_NUM_DTM, ARR_NUM_DTM_RANGE=ARR_NUM_DTM_RANGE, ARR_STR_DTM=ARR_STR_DTM)
        return ARR_OUT

    def MAKE_LAT_LON_ARR(self, FILE_NC_IN, STR_LAT="lat", STR_LON="lon", source="CFC"):
        """ Reading LAT and LON from a NC file """

        NC_DATA_IN = NC.Dataset(FILE_NC_IN, "r", format="NETCDF4")
        if source == "CFC":
            arr_lat_in = NC_DATA_IN.variables[STR_LAT]
            arr_lon_in = NC_DATA_IN.variables[STR_LON]
            num_nlat = len(arr_lat_in)
            num_nlon = len(arr_lon_in)
            arr_lon_out = [[0.0 for i in range(num_nlon)] for j in range(num_nlat)]
            arr_lat_out = [[0.0 for i in range(num_nlon)] for j in range(num_nlat)]
            for j in range(num_nlat):
                for i in range(num_nlon):
                    arr_lon_out[j][i] = arr_lat_in[j]
                    arr_lat_out[j][i] = arr_lon_in[i]
        return arr_lat_out, arr_lon_out

        
class NETCDF4_HELPER:
    def __init__(self):
        STR_NCDF4PY = NC.__version__
        print("Using netCDF4 for Python, Version: {0:s}".format(STR_NCDF4PY))

    def create_wrf_ensemble(self, STR_FILE_IN, STR_FILE_OUT, ARR_VAR=[], STR_DIR="./", NUM_ENSEMBLE_SIZE=1 ):
        FILE_OUT = NC.Dataset("{0:s}/{1:s}".format(STR_DIR, STR_FILE_OUT), "w",format="NETCDF4")
        FILE_IN  = NC.Dataset("{1:s}/{1:s}".format(STR_DIR, STR_FILE_IN ), "r",format="NETCDF4")
        
        # CREATE DIMENSIONS:
        for DIM in FILE_IN.dimensions:
            FILE_OUT.createDimension(DIM, FILE_IN.dimensions[DIM].size )
        FILE_OUT.createDimension("Ensembles", NUM_ENSEMBLE_SIZE )
        
        # CREATE ATTRIBUTES:
        FILE_OUT.TITLE                               = FILE_IN.TITLE                      
        FILE_OUT.START_DATE                          = FILE_IN.START_DATE    
        FILE_OUT.SIMULATION_START_DATE               = FILE_IN.SIMULATION_START_DATE               
        FILE_OUT.DX                                  = FILE_IN.DX
        FILE_OUT.DY                                  = FILE_IN.DY
        FILE_OUT.SKEBS_ON                            = FILE_IN.SKEBS_ON  
        FILE_OUT.SPEC_BDY_FINAL_MU                   = FILE_IN.SPEC_BDY_FINAL_MU           
        FILE_OUT.USE_Q_DIABATIC                      = FILE_IN.USE_Q_DIABATIC        
        FILE_OUT.GRIDTYPE                            = FILE_IN.GRIDTYPE  
        FILE_OUT.DIFF_OPT                            = FILE_IN.DIFF_OPT  
        FILE_OUT.KM_OPT                              = FILE_IN.KM_OPT
        
        if len(ARR_VAR) >0:
            for V in ARR_VAR:
                if V[1] == "2D":
                    FILE_OUT.createVariable(V[0],          "f8", ("Ensembles", "Time", "south_north", "west_east" ))
                elif V[1] == "3D":
                    FILE_OUT.createVariable(V[0],          "f8", ("Ensembles", "Time", "bottom_top", "south_north", "west_east" ))
        FILE_OUT.close() 
        FILE_IN.close()
        
    def add_ensemble(self, FILE_IN, FILE_OUT, STR_VAR, STR_DIM="2D", STR_DIR="./", IND_ENSEMBLE=0):
        FILE_OUT = NC.Dataset("{0:s}/{1:s}".format(STR_DIR, FILE_OUT), "a",format="NETCDF4")
        FILE_IN  = NC.Dataset("{0:s}/{1:s}".format(STR_DIR, FILE_IN ), "r",format="NETCDF4")
        
        ARR_VAR_IN = FILE_IN.variables[STR_VAR]
        NUM_NT = len(ARR_VAR_IN)
        NUM_NK = FILE_IN.dimensions["bottom_top"].size
        NUM_NJ = FILE_IN.dimensions["south_north"].size
        NUM_NI = FILE_IN.dimensions["west_east"].size
        for time in range(NUM_NT): 
            if STR_DIM   == "2D":
                FILE_OUT.variables[STR_VAR][IND_ENSEMBLE, time] = FILE_IN.variables[STR_VAR][time]
            elif STR_DIM == "3D":
                for k in range(NUM_NK):
                    FILE_OUT.variables[STR_VAR][IND_ENSEMBLE, time] = FILE_IN.variables[STR_VAR][time]
        FILE_OUT.close() 
        FILE_IN.close()

class WRF_TOOLS:
    STR_DIR_ROOT  = "./"
    NUM_TIME_INIT = 0
    NUM_SHIFT     = 0.001

    def __init__(self):
        """ 
        Remember: most array should be follow the rule of [j,i] instead of [x,y]. 
        """
        STR_NCDF4PY = NC.__version__
        print("Using netCDF4 for Python, Version: {0:s}".format(STR_NCDF4PY))


    def GEO_INFORMATER(self, STR_FILE="geo_em.d01.nc", STR_DIR=""):
        print("INPUT GEO FILE: {0:s}".format(STR_FILE))
        if STR_DIR == "":
            STR_DIR == self.STR_DIR_ROOT
        self.FILE_IN  = NC.Dataset("{0:s}/{1:s}".format(STR_DIR, STR_FILE ), "r",format="NETCDF4")
        self.MAP_LAT  = self.FILE_IN.variables["CLAT"] [self.NUM_TIME_INIT]
        self.MAP_LON  = self.FILE_IN.variables["CLONG"][self.NUM_TIME_INIT]
        ARR_TMP_IN     = self.FILE_IN.variables["CLONG"][0]
        # Since NetCDF4 for python does not support the hyphen in attributes, I 
        # am forced to calculate the NX and NY based on a map in the NC file. 
        self.NUM_NX    = len(ARR_TMP_IN[0])
        self.NUM_NY    = len(ARR_TMP_IN)

        self.NUM_DX   = self.FILE_IN.DX 
        self.NUM_DY   = self.FILE_IN.DX 
        
    def GEO_HELPER(self, ARR_LL_SW, ARR_LL_NE):
        self.MAP_CROP_MASK  = [[ 0 for i in range(self.NUM_NX)] for j in range(self.NUM_NY)]
        self.DIC_CROP_INFO = {"NE": {"LAT":0, "LON":0, "I":0, "J":0},\
                               "SW": {"LAT":0, "LON":0, "I":0, "J":0}}
        ARR_TMP_I   = []
        ARR_TMP_J   = []
        for j in range(self.NUM_NY):
            for i in range(self.NUM_NX):
                NUM_CHK_SW_J = self.MAP_LAT[j][i] - ARR_LL_SW[0]
                if NUM_CHK_SW_J == 0:
                    NUM_CHK_SW_J = self.MAP_LAT[j][i] - ARR_LL_SW[0] + self.NUM_SHIFT
                NUM_CHK_SW_I = self.MAP_LON[j][i] - ARR_LL_SW[1]
                if NUM_CHK_SW_I == 0:
                    NUM_CHK_SW_I = self.MAP_LAT[j][i] - ARR_LL_SW[1] - self.NUM_SHIFT
                NUM_CHK_NE_J = self.MAP_LAT[j][i] - ARR_LL_NE[0]
                if NUM_CHK_NE_J == 0:
                    NUM_CHK_NE_J = self.MAP_LAT[j][i] - ARR_LL_NE[0] + self.NUM_SHIFT
                NUM_CHK_NE_I = self.MAP_LON[j][i] - ARR_LL_NE[1]
                if NUM_CHK_NE_I == 0:
                    NUM_CHK_NE_I = self.MAP_LON[j][i] - ARR_LL_NE[1] - self.NUM_SHIFT

                NUM_CHK_NS_IN = NUM_CHK_SW_J * NUM_CHK_NE_J
                NUM_CHK_WE_IN = NUM_CHK_SW_I * NUM_CHK_NE_I
                if NUM_CHK_NS_IN < 0 and NUM_CHK_WE_IN < 0:
                    self.MAP_CROP_MASK[j][i] = 1
                    ARR_TMP_J.append(j) 
                    ARR_TMP_I.append(i)
        NUM_SW_J = min( ARR_TMP_J )
        NUM_SW_I = min( ARR_TMP_I )
        NUM_NE_J = max( ARR_TMP_J )
        NUM_NE_I = max( ARR_TMP_I )
        self.DIC_CROP_INFO["NE"]["J"]    = NUM_NE_J
        self.DIC_CROP_INFO["NE"]["I"]    = NUM_NE_I
        self.DIC_CROP_INFO["NE"]["LAT"]  = self.MAP_LAT[NUM_NE_J][NUM_NE_I]
        self.DIC_CROP_INFO["NE"]["LON"]  = self.MAP_LON[NUM_NE_J][NUM_NE_I]
        self.DIC_CROP_INFO["SW"]["J"]    = NUM_SW_J
        self.DIC_CROP_INFO["SW"]["I"]    = NUM_SW_I
        self.DIC_CROP_INFO["SW"]["LAT"]  = self.MAP_LAT[NUM_SW_J][NUM_SW_I]
        self.DIC_CROP_INFO["SW"]["LON"]  = self.MAP_LON[NUM_SW_J][NUM_SW_I]


    def PROFILE_HELPER(STR_FILE_IN, ARR_DATE_START, NUM_DOMS=3, NUM_TIMESTEPS=24, IF_PB=False):
        """
        This functions reads the filename, array of starting date, 
        and simulation hours and numbers of domains
        to profiling the time it takes for WRF.
        """
        FILE_READ_IN   = open("{0:s}".format(STR_FILE_IN))
        ARR_READ_IN    = FILE_READ_IN.readlines()
        NUM_TIME       = NUM_TIMESTEPS
        NUM_DOMAIN     = NUM_DOMS
        NUM_DATE_START = ARR_DATE_START
        NUM_LEN_IN     = len(ARR_READ_IN)

        ARR_TIME_PROFILE = [[0 for T in range(NUM_TIME)] for D in range(NUM_DOMS)]
        for I, TEXT_IN in enumerate(ARR_READ_IN):
            ARR_TEXT = re.split("\s",TEXT_IN.strip())
            if ARR_TEXT[0] == "Timing":
                if ARR_TEXT[2] == "main:" or ARR_TEXT[2] == "main":
                    for ind, T in enumerate(ARR_TEXT):
                        if T == "time"   : ind_time_text    = ind + 1
                        if T == "elapsed": ind_elapsed_text = ind - 1
                        if T == "domain" : ind_domain_text  = ind + 3
                    arr_time_in = re.split("_", ARR_TEXT[ind_time_text])
                    arr_date    = re.split("-", arr_time_in[0])
                    arr_time    = re.split(":", arr_time_in[1])
                    num_domain  = int(re.split(":", ARR_TEXT[ind_domain_text])[0])
                    num_elapsed = float(ARR_TEXT[ind_elapsed_text])
                    NUM_HOUR_FIX = (int(arr_date[2]) - NUM_DATE_START[2]) * 24
                    NUM_HOUR     = NUM_HOUR_FIX + int(arr_time[0])
                    ARR_TIME_PROFILE[num_domain-1][NUM_HOUR] += num_elapsed
            if IF_PB: TOOLS.progress_bar(I/float(NUM_LEN_IN))
        #self.ARR_TIME_PROFILE = ARR_TIME_PROFILE
        return ARR_TIME_PROFILE

    def wrf_reading_time(strIn):
        strDate, strTime = re.split("_", strIn)
        arrDate = re.split("-", strDate)
        arrTime = re.split(":", strTime)
        return [  int(arrDate[0]), int(arrDate[1]), int(arrDate[2]),\
                int(arrTime[0]), int(arrTime[1]), int(arrTime[2]),]


class DATA_READER:
    """
    The DATA_READER is based on my old work: gridtrans.py. 
    """
    def __init__(self, STR_NULL="noData", NUM_NULL=-999.999, STR_endian='b'):
        self.STR_NULL=STR_NULL
        self.NUM_NULL=NUM_NULL
        self.endian   = STR_endian
        if self.endian=="b":
           sym_end=">"
        elif self.endian=="l":
           sym_end="<"

    def stripblnk(arr,*num_typ):
        new_arr=[]
        for i in arr:
            if i == "":
                pass
            else:
                if num_typ[0] == 'int':
                    new_arr.append(int(i))
                elif num_typ[0] == 'float':
                    new_arr.append(float(i))
                elif num_typ[0] == '':
                    new_arr.append(i)
                else:
                    print("WRONG num_typ!")
        return new_arr
    
    def tryopen(self, sourcefile, ag):
        try:
            opf=open(sourcefile,ag)
            return opf
        except :
            print("No such file.")
            return "error"

    def READPFB (self, sourcefile):
         opf = tryopen(sourcefile,'rb')
         print("reading source file {0:s}".format(sourcefile))
         t1=struct.unpack('{0:s}ddd'.format(sym_end),opf.read(24))
         tn=struct.unpack('{0:s}iii'.format(sym_end),opf.read(12))
         td=struct.unpack('{0:s}ddd'.format(sym_end),opf.read(24))
         tns=struct.unpack('{0:s}i'.format(sym_end),opf.read(4))
         x1,y1,z1=t1
         nx,ny,nz=tn
         dx,dy,dz=td
         ns=tns[0]
         result_arr=list([[[0.0 for i in range(nx)] for j in range(ny)] for k in range(nz)])
         for isub in range(0,ns):
             ix,iy,iz,nnx,nny,nnz,rx,ry,rz=struct.unpack('{0:s}9i'.format(sym_end),opf.read(36))
             tmp_total = nnx * nny * nnz
             tvalue = struct.unpack('{1:s}{0:d}d'.format(tmp_total,sym_end), opf.read(8*tmp_total))
             for k in range(nnz):
                 for j in range(nny):
                     for i in range(nnx):
                                 result_arr[k+iz][j+iy][i+ix]=tvalue[ k*(nny*nnx) + j*nnx + i  ]
 
         opf.close()
         print("Completed reading pfb format from {0}".format(sourcefile))
         return result_arr,nx,ny,nz,dx,dy,dz
 
    def WRITEPFB(self, write_file_name,input_arr,nx,ny,nz,dx=0,dy=0,dz=0,x=0,y=0,z=0,ns=1):
        # Can choose for big endian type or little endian type
        sym_end = self.sym_end
        wtf=open(write_file_name,"w")
        wtf.write(struct.pack('{0:s}3d'.format(sym_end),x,y,z))
        wtf.write(struct.pack('{0:s}3i'.format(sym_end),nx,ny,nz))
        wtf.write(struct.pack('{0:s}3d'.format(sym_end),dx,dy,dz))
        wtf.write(struct.pack('{0:s}1i'.format(sym_end),ns))

        for isub in range(0,ns):
            iz,iy,ix=int(z),int(y),int(x)
            nnz,nny,nnx=int(nz),int(ny),int(nx)
            wtf.write(struct.pack('{0:s}3i'.format(sym_end),0,0,0))
            wtf.write(struct.pack('{0:s}3i'.format(sym_end),nx,ny,nz))
            wtf.write(struct.pack('{0:s}3i'.format(sym_end),dx,dy,dz))
            for i in range(iz,iz+nnz):
                for j in range(iy,iy+nny):
                    for k in range(ix,ix+nnx):
                        wtf.write(struct.pack('{0:s}d'.format(sym_end),input_arr[i][j][k]))
        wtf.close()

    def READCSV(self, sourcefile):
        opf    = self.tryopen(sourcefile,'r')
        opfchk = self.tryopen(sourcefile,'r')
        print("reading source file {0:s}".format(sourcefile))
        chk_lines = opfchk.readlines()
        num_totallines = len(chk_lines)
        ncols = 0
        num_notnum = 0
        for n in range(num_totallines):
            line_in = chk_lines[n]
            c_first = re.findall(".",line_in.strip())
            if c_first[0] == "#":
                num_notnum += 1
            else:
                ncols = len( re.split(",",line_in.strip()) )
                break
        if ncols == 0:
            print("something wrong with the input file! (all comments?)")
        else:
            del opfchk
            nrows=num_totallines - num_notnum
            result_arr=[[self.NUM_NULL for j in range(ncols)] for i in range(nrows)]
            result_arr_text=[]
            num_pass = 0
            for j in range(0,num_totallines):
                # chk if comment
                #print (j,i,chk_val)
                line_in = opf.readline()
                c_first = re.findall(".",line_in.strip())[0]
                if c_first == "#":
                    result_arr_text.append(line_in)
                    num_pass += 1
                else:
                    arr_in = re.split(",",line_in.strip())
                    for i in range(ncols):
                        chk_val = arr_in[i]
                        if chk_val == self.STR_NULL: 
                            result_arr[j-num_pass][i] = self.NUM_NULL
                        else:
                            result_arr[j-num_pass][i] = float(chk_val)
        return result_arr,result_arr_text

    def IbtracsReader(self, strFileIn):
        fileIn       = open(strFileIn, 'r')
        arrLines     = fileIn.readlines()
        numTimesteps = len(arrLines)

        arrITEM15  = ["BASIN", "TIME",   "NATURE", 
                      "LAT"  ,  "LON", "WMO_WIND", "WMO_PRES",
                      "USA_WIND", "USA_PRES", "TOKYO_WIND", "TOKYO_PRES",
                      "CMA_WIND", "CMA_PRES", "HKO_WIND"  ,   "HKO_PRES"]

        arrITEM16  = ["BASIN", "DATE", "TIME"    , "NATURE"  , 
                      "LAT"  ,  "LON", "WMO_WIND", "WMO_PRES",
                      "USA_WIND", "USA_PRES", "TOKYO_WIND", "TOKYO_PRES",
                      "CMA_WIND", "CMA_PRES",   "HKO_WIND",   "HKO_PRES"]

        arrFMT15   = ["str", "arr", "str", "float", "float",
                      "float", "float", "float", "float", "float", 
                      "float", "float", "float", "float", "float", "float"] 

        arrFMT16   = ["str", "arr", "arr", "str", "float", "float",
                      "float", "float", "float", "float", "float", 
                      "float", "float", "float", "float", "float", "float"] 

        # Made the output dic
        dicOut    = {}
        dicOut["TIMESTEPS"] = numTimesteps

        for item in arrITEM16:
            dicOut[item] = [ "" for n in range(numTimesteps) ]

        arrDATE = []
        # Analysing the dat and into the dictionary output
        for ind,line in enumerate(arrLines):
            arr_tmp = re.split("\s", line)
            numChk  = len(arr_tmp)
            if numChk == 17:
                for ind2 in range(16):
                    if arrFMT16[ind2] == "str":
                        dicOut[ arrITEM16[ind2] ][ind] = str(arr_tmp[ind2]) 
                    elif arrFMT16[ind2] == "arr":
                        dicOut[ arrITEM16[ind2] ][ind] = \
                            InputTool.readDateTime(arr_tmp[ind2])[arrITEM16[ind2]]
                    elif arrFMT16[ind2] == "float":
                        if arr_tmp[ind2] == '':
                            dicOut[ arrITEM16[ind2] ][ ind ] = 0.0
                        else:
                            dicOut[ arrITEM16[ind2] ][ ind ] = float(arr_tmp[ ind2 ]) 
                    if arrITEM16[ind2] == "DATE":
                        arrDATE = InputTool.readDateTime(arr_tmp[ind2])[arrITEM16[ind2]]
            elif numChk == 16:
                for ind2 in range(15):
                    if arrFMT15[ind2] == "str":
                        dicOut[ arrITEM15[ind2] ][ind] = str(arr_tmp[ind2]) 
                    elif arrFMT15[ind2] == "arr":
                        dicOut[ arrITEM15[ind2] ][ind] = \
                            InputTool.readDateTime(arr_tmp[ind2])[arrITEM15[ind2]] 
                    elif arrFMT15[ind2] == "float":
                        if arr_tmp[ind2] == '':
                            dicOut[ arrITEM15[ind2] ][ ind ] = -999.999
                        else:
                            dicOut[ arrITEM15[ind2] ][ ind ] = float(arr_tmp[ ind2 ]) 
                    if arrITEM16[ind2] == "DATE":
                        dicOut[ "DATE" ][ ind ] = arrDATE             
            else:
                print("something wrong with the dat file")
        return dicOut  

class InputTool:
    def readDateTime(strIn):
        arrTmp1 = re.split("-", strIn) 
        arrTmp2 = re.split(":", strIn) 
        if len(arrTmp1) == 3:
            return {"DATE" : [int(tmp) for tmp in arrTmp1]  }
        elif len(arrTmp2) == 3:
            return {"TIME" : [int(tmp) for tmp in arrTmp2]  }
        else:
            return {"NONE" : None }
class AGORITHM:
    def ParticleFilter(ARR_IN):
        NUM_LEN_ARR_IN    = len(ARR_IN)
        NPARR_IN          = NP.array(ARR_IN)
        NPARR_WEIGHT_OUT  = NPARR_IN / NPARR_IN.sum()
        NUM_WEIGHT_CUMSUM = NPARR_WEIGHT_OUT.cumsum()
        NUM_THRE          = 1./len(NUM_WEIGHT_CUMSUM)
        NUM_START         = RD.random() * NUM_THRE
        ARR_THRE          = NP.array([ NUM_THRE for n in range(len(NUM_WEIGHT_CUMSUM))])
        ARR_THRE [0]      = NUM_START
        ARR_CUM_THRE      = ARR_THRE.cumsum()
        ARR_WEIGHT2       = NP.insert(NPARR_WEIGHT_OUT, 0, 0)
        ARR_IND_NUMBERS = NP.array([0 for n in range(NUM_LEN_ARR_IN)])
        for ind, num in enumerate(ARR_WEIGHT2):
            for num_in in ARR_CUM_THRE:
                if ind < NUM_LEN_ARR_IN:
                    num_chk = ( ARR_WEIGHT2[ind] - num_in ) * (ARR_WEIGHT2[ind+1] - num_in)
                    if num_chk < 0:
                        ARR_IND_NUMBERS[ind] += 1
        ARR_IND_CONFIRM = NP.array([ int(n >= 1) for n in ARR_IND_NUMBERS])

        return ARR_IND_NUMBERS, ARR_IND_CONFIRM

    def ParticleFilterOrigin(ARR_IN):
        NUM_ENSEMBLE   = len(ARR_IN)
        NPARR_PROB_OUT = NP.array(ARR_IN )
        NUM_MAX_PROB      = NPARR_PROB_OUT.max()
        NUM_MIN_PROB      = NPARR_PROB_OUT.min()
        NPARR_WEIGHT_OUT  = NPARR_PROB_OUT/NPARR_PROB_OUT.sum()

        cumsum_weight = NPARR_WEIGHT_OUT.cumsum()
        num_thre      = 1./len(cumsum_weight)
        num_start     = (RD.random() ) * num_thre
        arr_thre      = NP.array([num_thre for n in range(len(cumsum_weight))])
        arr_thre[0]   = num_start
        arr_cum_thre  = arr_thre.cumsum()

        cumsum_weight2 = NP.insert(cumsum_weight,0,0)

        ARR_IND_NUMBERS = NP.array([0 for n in range(NUM_ENSEMBLE)])
        for ind, num in enumerate(cumsum_weight2):
            for num_in in arr_cum_thre:
                if ind < NUM_ENSEMBLE:
                    num_chk = ( cumsum_weight2[ind] - num_in ) * (cumsum_weight2[ind+1] - num_in)
                    if num_chk < 0:
                        ARR_IND_NUMBERS[ind] += 1

        ARR_IND_CONFIRM = NP.array([ int(n >= 1) for n in ARR_IND_NUMBERS])

        return ARR_IND_NUMBERS, ARR_IND_CONFIRM


