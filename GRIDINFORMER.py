import math,re,sys
import random as RD
import time
try:
    import netCDF4 as NC
except:
    print("Do not import netCDF4")
try:
    import numpy as NP
except:
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
    # FROM WRF: module_cam_shr_const_mod.f90
    NUM_CONST_EARTH_R = 6.37122E6 
    NUM_CONST_PI      = 3.14159265358979323846
    
    
    def __init__(self, name="GRID", ARR_LAT=[], ARR_LON=[], DIMENSIONS=2 ):
        self.STR_NAME = name
        self.NUM_DIMENSIONS = DIMENSIONS
        self.NUM_LAST_INDEX = -1
        self.ARR_GRID = []
        
        self.ARR_LAT  = ARR_LAT
        self.ARR_LON  = ARR_LON
       
        if len(ARR_LAT) != 0 and len(ARR_LON) != 0:
            NUM_ARR_NY_T1 = len(ARR_LAT)
            NUM_ARR_NY_T2 = len(ARR_LON)
            NUM_ARR_NX_T1 = len(ARR_LAT[0])
            NUM_ARR_NX_T2 = len(ARR_LON[0])
            self.NUM_ARR_NX = NUM_ARR_NX_T1
            self.NUM_ARR_NY = NUM_ARR_NY_T1
            if NUM_ARR_NY_T1 - NUM_ARR_NY_T2 + NUM_ARR_NX_T1 - NUM_ARR_NX_T2 != 0:
                print("The gridcell of LAT is {0:d}&{1:d}, and LON is {2:d}&{3:d} are not match"\
                      .format(NUM_ARR_NY_T1,NUM_ARR_NY_T2,NUM_ARR_NX_T1,NUM_ARR_NX_T2))
    def index_map(self):
        self.INDEX_MAP = [[0 for i in range(self.NUM_ARR_NX)] for j in range(self.NUM_ARR_NY)]
        NUM_ALL_INDEX = len(self.ARR_GRID)
        for n in range(NUM_ALL_INDEX):
            self.INDEX_MAP[self.ARR_GRID[n]["INDEX_J"]][self.ARR_GRID[n]["INDEX_I"]] =\
            self.ARR_GRID[n]["INDEX"]
                                                     
    def add_an_element(self, ARR_GRID, NUM_INDEX=0, STR_VALUE=STR_VALUE_INIT, NUM_VALUE=NUM_VALUE_INIT ):
        OBJ_ELEMENT = {"INDEX" : NUM_INDEX, \
                       STR_VALUE : NUM_VALUE}
        ARR_GRID.append(OBJ_ELEMENT)
        
    def add_an_geo_element(self, ARR_GRID, NUM_INDEX=-999, NUM_J=0, NUM_I=0, \
                           NUM_MAX_NX = 0, NUM_MAX_NY = 0, \
                           STR_VALUE=STR_VALUE_INIT, NUM_VALUE=NUM_VALUE_INIT ):
        if NUM_MAX_NX == 0 or NUM_MAX_NY == 0:
            NUM_MAX_NX = self.NUM_ARR_NX 
            NUM_MAX_NY = self.NUM_ARR_NY
        NUM_CENTER_LON = self.ARR_LON[NUM_J][NUM_I] 
        NUM_CENTER_LAT = self.ARR_LAT[NUM_J][NUM_I]
        if NUM_I == 0:
            NUM_WE_LON =      ( self.ARR_LON[NUM_J][NUM_I] - self.ARR_LON[NUM_J][NUM_I + 1] ) * 0.5
            NUM_EW_LON = -1 * ( self.ARR_LON[NUM_J][NUM_I] - self.ARR_LON[NUM_J][NUM_I + 1] ) * 0.5
        elif NUM_I == NUM_MAX_NX - 1:
            NUM_WE_LON = -1 * ( self.ARR_LON[NUM_J][NUM_I] - self.ARR_LON[NUM_J][NUM_I - 1] ) * 0.5
            NUM_EW_LON =      ( self.ARR_LON[NUM_J][NUM_I] - self.ARR_LON[NUM_J][NUM_I - 1] ) * 0.5
        else:
            NUM_WE_LON = ( self.ARR_LON[NUM_J][NUM_I] - self.ARR_LON[NUM_J][NUM_I + 1] ) * 0.5
            NUM_EW_LON = ( self.ARR_LON[NUM_J][NUM_I] - self.ARR_LON[NUM_J][NUM_I - 1] ) * 0.5
        if NUM_J == 0:
            NUM_SN_LAT = -1 * ( self.ARR_LAT[NUM_J][NUM_I] - self.ARR_LAT[NUM_J + 1][NUM_I ] ) * 0.5
            NUM_NS_LAT =      ( self.ARR_LAT[NUM_J][NUM_I] - self.ARR_LAT[NUM_J + 1][NUM_I ] ) * 0.5
        elif NUM_J == NUM_MAX_NY - 1:
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
                                   "E": NUM_CENTER_LON + NUM_EW_LON,"W": NUM_CENTER_LON + NUM_WE_LON},\
                       STR_VALUE : NUM_VALUE}
        ARR_GRID.append(OBJ_ELEMENT)
        
    def add_an_geo_variable(self, ARR_GRID, NUM_INDEX=-999, NUM_J=0, NUM_I=0,\
                           STR_VALUE=STR_VALUE_INIT, NUM_VALUE=NUM_VALUE_INIT ):
        if NUM_INDEX == -999:
            NUM_INDEX = self.INDEX_MAP[NUM_J][NUM_I]
        ARR_GRID[NUM_INDEX][STR_VALUE] = NUM_VALUE
        
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
                
    def create_resample_map(self):
        NUM_IND = 0
        self.ARR_RESAMPLE_MAP = []
        for j in range(self.NUM_GRIDS_LAT):
            for i in range(self.NUM_GRIDS_LAT):
                self.add_an_geo_element(self.ARR_RESAMPLE_MAP, NUM_INDEX=NUM_IND, NUM_J=j, NUM_I=i, \
                           NUM_MAX_NX= self.NUM_GRIDS_LON, NUM_MAX_NY= self.NUM_GRIDS_LAT)
                NUM_IND += 1
        self.NUM_MAX_INDEX_RS = NUM_IND
                
    def create_reference_map(self, MAP_TARGET, MAP_RESAMPLE, STR_TYPE="FIX", NUM_SHIFT=0.001, IF_PB=False):
        """Must input with OBJ_REFERENCE"""
        self.ARR_REFERENCE_MAP = []
        if STR_TYPE=="FIX":
            NUM_OBJ_G_LEN = len(MAP_TARGET)
            for OBJ_G in MAP_TARGET:
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
                if IF_PB: TOOLS.progress_bar(TOOLS.loop_progress_cal([OBJ_G["INDEX"]], [NUM_OBJ_G_LEN]), STR_DES="CREATING REFERENCE MAP")
                
    def export_reference_map(self, STR_DIR, STR_FILENAME, STR_TYPE="netCDF4", IF_PB=False ):
        TIME_NOW = time.gmtime()
        self.STR_DATE_NOW = "{0:04d}-{1:02d}-{2:02d}".format(TIME_NOW.tm_year, TIME_NOW.tm_mon, TIME_NOW.tm_mday) 
        self.STR_TIME_NOW = "{0:04d}:{1:02d}:{2:02d}".format(TIME_NOW.tm_hour, TIME_NOW.tm_min, TIME_NOW.tm_sec)

        if STR_TYPE == "netCDF4":
            NCDF4_DATA = NC.Dataset("{0:s}/{1:s}".format(STR_DIR, STR_FILENAME), 'w', format="NETCDF4")
            # CREATE ATTRIBUTEs:
            NCDF4_DATA.description = \
            "The netCDF4 version of reference map which contains grid information for resampling"
            NCDF4_DATA.history = "Create on {0:s} at {1:s}".format(self.STR_DATE_NOW, self.STR_TIME_NOW)
            
            # CREATE DIMENSIONs:
            NCDF4_DATA.createDimension("Y",self.NUM_ARR_NY)
            NCDF4_DATA.createDimension("X",self.NUM_ARR_NX)
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
            NUM_MAX_I     = self.NUM_ARR_NX
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
                if IF_PB: TOOLS.progress_bar((i+j*NUM_MAX_I)/float(NUM_TOTAL_OBJ), STR_DES="CREATING REFERENCE MAP")
            NCDF4_DATA.close()
        
    def import_reference_map(self, STR_DIR, STR_FILENAME, STR_TYPE="netCDF4", IF_PB=False):
        self.ARR_REFERENCE_MAP = []
        self.NUM_MAX_INDEX_RS = 0
        self.NUM_MIN_INDEX_RS = 999
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
            
            for j in range(self.REFERENCE_MAP_NY):
                for i in range(self.REFERENCE_MAP_NX):
                    OBJ_ELEMENT = {"INDEX"       :                         0 ,\
                                   "INDEX_I"     :                         0 ,\
                                   "INDEX_J"     :                         0 ,\
                                   "CENTER"      :  {"LAT": 0.0, "LON": 0.0} ,\
                                   "INDEX_REF"   :                         0 ,\
                                   "INDEX_REF_I" :                         0 ,\
                                   "INDEX_REF_J" :                         0 ,\
                                   "CENTER_REF"  :  {"LAT": 0.0, "LON": 0.0} }
                    OBJ_ELEMENT["INDEX"]       = INDEX      [j][i] 
                    OBJ_ELEMENT["INDEX_I"]     = INDEX_J    [j][i]
                    OBJ_ELEMENT["INDEX_J"]     = INDEX_I    [j][i]
                    OBJ_ELEMENT["INDEX_REF"]   = INDEX_REF  [j][i]
                    OBJ_ELEMENT["INDEX_REF_J"] = INDEX_REF_J[j][i]
                    OBJ_ELEMENT["INDEX_REF_I"] = INDEX_REF_I[j][i]
                    OBJ_ELEMENT["CENTER"]["LAT"]     = CENTER_LAT    [j][i]
                    OBJ_ELEMENT["CENTER"]["LON"]     = CENTER_LON    [j][i]
                    OBJ_ELEMENT["CENTER_REF"]["LAT"] = CENTER_REF_LAT[j][i]
                    OBJ_ELEMENT["CENTER_REF"]["LON"] = CENTER_REF_LON[j][i]
                    self.ARR_REFERENCE_MAP.append(OBJ_ELEMENT)
                    self.NUM_MIN_INDEX_RS = min(self.NUM_MIN_INDEX_RS, INDEX_REF[j][i])
                    self.NUM_MAX_INDEX_RS = max(self.NUM_MAX_INDEX_RS, INDEX_REF[j][i])
                if IF_PB: TOOLS.progress_bar((j+1)/float(self.REFERENCE_MAP_NY), STR_DES="IMPORTING ... ")
            if self.NUM_MIN_INDEX_RS == 0:
                self.NUM_MAX_RS = self.NUM_MAX_INDEX_RS + 1
            NCDF4_DATA.close()    
        
    def resample_calculation(self, ARR_REFERENCE_MAP, ARR_VARIABLES, ARR_GRID_IN=[], STR_TYPE="ALL",\
                             DIC_THRESHOLD={"0":0.1,"1":0.9, "2":1.0},\
                             DIC_PERCENTILE={ "P05": 0.05, "P10": 0.1, "P25": 0.25, "P75": 0.75, "P90": 0.90, "P95": 0.95},\
                             IF_PB=False):
        """
        STR_TYPE includes:
           ALL : MEAN, MEDIAN, 5,10,25,75,90,95 Percentile
           MASK: (This is done by identified by DIC_THRESHOLD)     
        """
        if ARR_GRID_IN == []:
            ARR_GRID_IN = self.ARR_GRID
        self.ARR_RESAMPLE_OUT = [ {} for n in range(self.NUM_MAX_RS)]
        for IND in range(len(ARR_GRID_IN)):
            IND_REF2RS = ARR_REFERENCE_MAP[IND]["INDEX_REF"]
            for VAR in ARR_VARIABLES:
                OBJ_ELEMENT = {"VALUE"   :  []}
                self.ARR_RESAMPLE_OUT[IND_REF2RS][VAR] = OBJ_ELEMENT
        for IND in range(len(ARR_GRID_IN)):
            IND_REF2RS = ARR_REFERENCE_MAP[IND]["INDEX_REF"]
            for VAR in ARR_VARIABLES:
                self.ARR_RESAMPLE_OUT[IND_REF2RS][VAR]["VALUE"].append(ARR_GRID_IN[IND][VAR])
            if IF_PB: TOOLS.progress_bar(loop_progress_cal([IND], [len(ARR_GRID_IN)]), STR_DES="RESAMPLING PROGRESS")
        NUM_RS_OUT_LEN = len(self.ARR_RESAMPLE_OUT)
        for IND in range(NUM_RS_OUT_LEN):
            for VAR in ARR_VARIABLES:
                ARR_IN          = self.ARR_RESAMPLE_OUT[IND][VAR]["VALUE"]
                ARR_IN.sort()
                NUM_ARR_LEN     = len(ARR_IN)
                NUM_ARR_MEAN    = sum(ARR_IN) / float(NUM_ARR_LEN)
                NUM_ARR_S2SUM   = 0
                if math.fmod(len(ARR_IN),2) == 1:
                    NUM_MPOS = [int((NUM_ARR_LEN-1)/2), int((NUM_ARR_LEN-1)/2)]
                else:
                    NUM_MPOS = [int(NUM_ARR_LEN/2), int(NUM_ARR_LEN/2 -1) ]
                self.ARR_RESAMPLE_OUT[IND][VAR]["MIN"]     = min(ARR_IN)
                self.ARR_RESAMPLE_OUT[IND][VAR]["MAX"]     = max(ARR_IN)
                self.ARR_RESAMPLE_OUT[IND][VAR]["MEAN"]    = NUM_ARR_MEAN 
                self.ARR_RESAMPLE_OUT[IND][VAR]["MEDIAN"]  = (ARR_IN[NUM_MPOS[0]] + ARR_IN[NUM_MPOS[1]]) / 2.0
                for STVA in DIC_PERCENTILE:
                    self.ARR_RESAMPLE_OUT[IND][VAR][STVA]  = ARR_IN[ round(NUM_ARR_LEN * DIC_PERCENTILE[STVA])-1]
                for VAL in ARR_IN:
                    NUM_ARR_S2SUM += (VAL - NUM_ARR_MEAN)**2
                self.ARR_RESAMPLE_OUT[IND][VAR]["STD"]     =  (NUM_ARR_S2SUM / (NUM_ARR_LEN-1))**0.5
            if IF_PB: TOOLS.progress_bar(TOOLS.loop_progress_cal([IND], [NUM_RS_OUT_LEN]), STR_DES="CALCULATION PROGRESS")


class TOOLS:
    def progress_bar(NUM_PROGRESS, NUM_PROGRESS_BIN=0.05, STR_SYS_SYMBOL="=", STR_DES="Progress"):
        NUM_SYM = int(NUM_PROGRESS / NUM_PROGRESS_BIN)
        sys.stdout.write('\r')
        sys.stdout.write('[{0:20s}]{1:4.1f}% {2:s}'.format(STR_SYS_SYMBOL*NUM_SYM, NUM_PROGRESS*100, STR_DES))
        sys.stdout.flush()

    def loop_progress_cal(ARR_INDEX, ARR_INDEX_MAX, NUM_CUM_MAX=1, NUM_CUM_IND=1, NUM_TOTAL_MAX=1):
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

# MAKING Calandar calculator
    def calendar_cal(ARR_START_TIME, ARR_INTERVAL, ARR_END_TIME_IN=[0, 0, 0, 0, 0, 0.0], IF_LEAP=False):
        # DICTIONARY VERSION
        ARR_END_TIME  = [ 0,0,0,0,0,0.0]
        ARR_DATETIME  = ["SECOND", "MINUTE", "HOUR","DAY", "MON", "YEAR"]
        NUM_ARR_DATETIME = len(ARR_DATETIME)
        if IF_LEAP:
            ARR_DAY_LIM = [0,31,29,31,30,31,30,31,31,30,31,30,31]    
        else:
            ARR_DAY_LIM = [0,31,28,31,30,31,30,31,31,30,31,30,31]
        
        DIC_TIME_LIM = \
        {"YEAR": {"START": 0 , "LIMIT": 999999 },\
         "MON": {"START": 1 , "LIMIT": 12 },\
         "DAY": {"START": 1 , "LIMIT": 31 },\
         "HOUR": {"START": 0 , "LIMIT": 23 },\
         "MINUTE": {"START": 0 , "LIMIT": 59 },\
         "SECOND": {"START": 0 , "LIMIT": 59 },\
        }
        for I, T in enumerate(ARR_START_TIME):
            ARR_END_TIME[I] += T + ARR_INTERVAL[I]
        if math.fmod(ARR_START_TIME[0], 4) == 0: IF_LEAP = True 
        for I, ITEM in enumerate(ARR_DATETIME):
            NUM_ARR_POS = NUM_ARR_DATETIME-I-1
            if ITEM == "DAY":
                if ARR_END_TIME[NUM_ARR_POS] > ARR_DAY_LIM[ARR_END_TIME[1]]:
                    ARR_END_TIME[NUM_ARR_POS - 1] += 1
                    ARR_END_TIME[NUM_ARR_POS] = DIC_TIME_LIM[ITEM]["START"]
            else:
                if ARR_END_TIME[NUM_ARR_POS] > DIC_TIME_LIM[ITEM]["LIMIT"]:
                    ARR_END_TIME[NUM_ARR_POS - 1] += 1
                    ARR_END_TIME[NUM_ARR_POS] = DIC_TIME_LIM[ITEM]["START"]
            
        return ARR_END_TIME
    


class GEO_TOOLS:
    def MAKE_LAT_LON_ARR(FILE_NC_IN, source="CFC"):
        NC_DATA_IN = NC.Dataset(FILE_NC_IN, "r", format="NETCDF4")
        if source == "CFC":
            arr_lat_in = NC_DATA_IN.variables["lat"]
            arr_lon_in = NC_DATA_IN.variables["lon"]
            num_nlat = len(arr_lat_in)
            num_nlon = len(arr_lon_in)
            arr_lon_out = [[0.0 for i in range(num_nlon)] for j in range(num_nlat)]
            arr_lat_out = [[0.0 for i in range(num_nlon)] for j in range(num_nlat)]
            for j in range(num_nlat):
                for i in range(num_nlon):
                    arr_lon_out[j][i] = arr_lat_in[j]
                    arr_lat_out[j][i] = arr_lon_in[i]
        return arr_lat_out, arr_lon_out



