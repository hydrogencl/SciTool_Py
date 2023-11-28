import netCDF4 as NC
import sys, re, math

def make_new_nc(STR_FILE_IN, STR_FILE_OUT, ARR_VAR=[], STR_DIR="./", NUM_ENSEMBLE_SIZE=1,\
                NUM_MPI_RANK=0, NUM_MPI_SIZE=0):
    FILE_OUT = NC.Dataset("{0:s}/{1:s}".format(STR_DIR, STR_FILE_OUT), "w",format="NETCDF4")
    FILE_IN  = NC.Dataset("{0:s}/{1:s}".format(STR_DIR, STR_FILE_IN ), "r",format="NETCDF4")
    
    # CREATE DIMENSIONS:
    
    for DIM in FILE_IN.dimensions:
        FILE_OUT.createDimension(DIM, FILE_IN.dimensions[DIM].size )
    FILE_OUT.createDimension("Ensembles", NUM_ENSEMBLE_SIZE )
    
    # CREATE ATTRIBUTES:
    if NUM_MPI_SIZE > 0:
        FILE_OUT.MPI_RANK = NUM_MPI_RANK
        FILE_OUT.MPI_SIZE = NUM_MPI_SIZE

    FILE_OUT.TITLE                               = FILE_IN.TITLE                      
    FILE_OUT.START_DATE                          = FILE_IN.START_DATE    
    FILE_OUT.SIMULATION_START_DATE               = FILE_IN.SIMULATION_START_DATE               
    #FILE_OUT.WEST-EAST_GRID_DIMENSION            = FILE_IN.WEST-EAST_GRID_DIMENSION                  
    #FILE_OUT.SOUTH-NORTH_GRID_DIMENSION          = FILE_IN.SOUTH-NORTH_GRID_DIMENSION                    
    #FILE_OUT.BOTTOM-TOP_GRID_DIMENSION           = FILE_IN.BOTTOM-TOP_GRID_DIMENSION                   
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

def read_to_new_nc(FILE_IN, FILE_OUT, STR_VAR, STR_DIM="2D", STR_DIR="./", IND_ENSEMBLE=0):
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
                    



