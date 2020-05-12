class GRIDINFORMATER:
    """
    This object is the information of the input gridcells/array/map. 
    Using .an_element to add an element/gridcell    
    """
    STR_VALUE_INIT = "None"
    NUM_VALUE_INIT = -9999.9
    def __init__(self, name="GRID", dimensions=2, arr_lat=[], arr_lon=[] ):
        self.STR_NAME = name
        self.NUM_DIMENSIONS = dimensions
        self.NUM_LAST_INDEX = 0
        self.ARR_GRID = []
        self.ARR_LAT  = arr_lat
        self.ARR_LON  = arr_lon
        
        if  len(arr_lat) != 0 and len(arr_lon) != 0:
            NUM_ARR_NY_T1 = len(arr_lat)
            NUM_ARR_NY_T2 = len(arr_lon)
            NUM_ARR_NX_T1 = len(arr_lat[0])
            NUM_ARR_NX_T2 = len(arr_lon[0])
            self.NUM_ARR_NX = NUM_ARR_NX_T1
            self.NUM_ARR_NY = NUM_ARR_NY_T1
            if NUM_ARR_NY_T1 - NUM_ARR_NY_T2 + NUM_ARR_NX_T1 - NUM_ARR_NX_T2 != 0:
                print("The gridcell of LAT is {0:d}&{1:d}, and LON is {2:d}&{3:d} and wrong"\
                      .format(NUM_ARR_NY_T1,NUM_ARR_NY_T2,NUM_ARR_NX_T1,NUM_ARR_NX_T2))
    def add_an_element(self, NUM_INDEX=0, STR_VALUE=STR_VALUE_INIT, NUM_VALUE=NUM_VALUE_INIT ):
        OBJ_ELEMENT = {"INDEX" : NUM_INDEX, \
                       STR_VALUE : NUM_VALUE}
        self.ARR_GRID.append(OBJ_ELEMENT)
        
    def add_an_geo_element(self, NUM_INDEX=-999, num_j=0, num_i=0, STR_VALUE=STR_VALUE_INIT, NUM_VALUE=NUM_VALUE_INIT ):
        NUM_CENTER_LON = self.ARR_LON[num_j][num_i] 
        NUM_CENTER_LAT = self.ARR_LAT[num_j][num_i]
        if num_i == 0:
            NUM_WE_LON =      ( self.ARR_LON[num_j][num_i] - self.ARR_LON[num_j][num_i + 1] ) * 0.5
            NUM_EW_LON = -1 * ( self.ARR_LON[num_j][num_i] - self.ARR_LON[num_j][num_i + 1] ) * 0.5
        elif num_i == self.NUM_ARR_NX - 1:
            NUM_WE_LON = -1 * ( self.ARR_LON[num_j][num_i] - self.ARR_LON[num_j][num_i - 1] ) * 0.5
            NUM_EW_LON =      ( self.ARR_LON[num_j][num_i] - self.ARR_LON[num_j][num_i - 1] ) * 0.5
        else:
            NUM_WE_LON = ( self.ARR_LON[num_j][num_i] - self.ARR_LON[num_j][num_i + 1] ) * 0.5
            NUM_EW_LON = ( self.ARR_LON[num_j][num_i] - self.ARR_LON[num_j][num_i - 1] ) * 0.5
        if num_j == 0:
            NUM_SN_LAT = -1 * ( self.ARR_LAT[num_j][num_i] - self.ARR_LAT[num_j + 1][num_i ] ) * 0.5
            NUM_NS_LAT =      ( self.ARR_LAT[num_j][num_i] - self.ARR_LAT[num_j + 1][num_i ] ) * 0.5
        elif num_j == self.NUM_ARR_NY - 1:
            NUM_SN_LAT =      ( self.ARR_LAT[num_j][num_i] - self.ARR_LAT[num_j - 1][num_i ] ) * 0.5 
            NUM_NS_LAT = -1 * ( self.ARR_LAT[num_j][num_i] - self.ARR_LAT[num_j - 1][num_i ] ) * 0.5
        else:       
            NUM_SN_LAT = ( self.ARR_LAT[num_j][num_i] - self.ARR_LAT[num_j - 1][num_i ] ) * 0.5
            NUM_NS_LAT = ( self.ARR_LAT[num_j][num_i] - self.ARR_LAT[num_j + 1][num_i ] ) * 0.5
        ARR_NE = [ NUM_CENTER_LON + NUM_EW_LON , NUM_CENTER_LAT + NUM_NS_LAT ]
        ARR_NW = [ NUM_CENTER_LON + NUM_WE_LON , NUM_CENTER_LAT + NUM_NS_LAT ] 
        ARR_SE = [ NUM_CENTER_LON + NUM_EW_LON , NUM_CENTER_LAT + NUM_SN_LAT ] 
        ARR_SW = [ NUM_CENTER_LON + NUM_WE_LON , NUM_CENTER_LAT + NUM_SN_LAT ] 
        if NUM_INDEX == -999: 
            NUM_INDEX = self.NUM_LAST_INDEX +1
            self.NUM_LAST_INDEX += 1
        OBJ_ELEMENT = {"INDEX"      : NUM_INDEX,\
                       "INDEX_I"    : num_i,\
                       "INDEX_J"    : num_j,\
                       "CENTER" : {"LAT" : NUM_CENTER_LAT, "LON" : NUM_CENTER_LON},\
                       "VERTEX" : {"NE": ARR_NE, "SE": ARR_SE, "SW": ARR_SW, "NW": ARR_NW},\
                       STR_VALUE : NUM_VALUE}
        self.ARR_GRID.append(OBJ_ELEMENT)

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



