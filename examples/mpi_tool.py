### 
#  This is the example for using MPI_TOOLS
#     1) Creating the geometry for mpi rank for the input array
#     2) Creating the netCDF file
###

import math, sys,os, time, re
import netCDF4 as NC
from mpi4py import MPI
from GRIDINFORMER import MPI_TOOLS as MT

IF_PARALLEL = True

# MPI configuration
# Initialing the MPI Commworld
comm     = MPI.COMM_WORLD
# Get the rank number (aka index of the CPU)
NUM_MPI_RANK = comm.Get_rank()
# Get the total size of CPUs
NUM_MPI_SIZE = comm.Get_size()

numNX = 200
numNY = 100
numNT =   6

# Design Cores, can be the same as NUM_MPI_SIZE
# numCores   = 48
numCores   = NUM_MPI_SIZE

numCoresNX =  6

MPI_SET = MT(MPI_SIZE=NUM_MPI_SIZE, MPI_RANK = NUM_MPI_RANK,\
             NUM_NX_END=numNX ,\
             NUM_NY_END=numNY ,\
             NUM_NX_CORES=numCoresNX)

MPI_SET.CPU_GEOMETRY_2D()
ARR_RANK_DESIGN = MPI_SET.ARR_RANK_DESIGN

# Preparing the NC file
STR_OUTPUT_FILE_TMP = "./test_mpi.nc" 

if not os.path.exists(STR_OUTPUT_FILE_TMP):
    if NUM_MPI_RANK == 0:

        MPI_SET.MPI_MESSAGE("Starting to create File ")
        NC_OUT_TMP         = NC.Dataset(STR_OUTPUT_FILE_TMP, 'w')
        # Copying Dimensions:
        NC_OUT_TMP.createDimension( "NX", numNX            )
        NC_OUT_TMP.createDimension( "NY", numNY            )
        NC_OUT_TMP.createDimension( "NT", numNT            )

        # Adding statistics name:
        NC_OUT_TMP.createVariable ('MAP2D', "f4", ('NY', 'NX'                ) )
        NC_OUT_TMP.createVariable ('MAP3D', "i8", ('NT', 'NY', 'NX'          ) )
        MPI_SET.MPI_MESSAGE("Finishing Creating the NC output{}".format(STR_OUTPUT_FILE_TMP))
        NC_OUT_TMP.close()
    else: 
        time.sleep(3)
        MPI_SET.MPI_MESSAGE("Done Waiting")

    
    NC_OUT     = NC.Dataset(STR_OUTPUT_FILE_TMP, 'a', parallel=IF_PARALLEL)
else:
    NC_OUT     = NC.Dataset(STR_OUTPUT_FILE_TMP, 'a', parallel=IF_PARALLEL)
  

MPI_SET.MPI_MESSAGE("Starting to work!")

NUM_NX_START = ARR_RANK_DESIGN[NUM_MPI_RANK]["NX_START"]
NUM_NX_END   = ARR_RANK_DESIGN[NUM_MPI_RANK]["NX_END"]
NUM_NY_START = ARR_RANK_DESIGN[NUM_MPI_RANK]["NY_START"]
NUM_NY_END   = ARR_RANK_DESIGN[NUM_MPI_RANK]["NY_END"]
print(NUM_MPI_RANK, NUM_NX_END, NUM_NY_END)
VAR_IN2D = NC_OUT.variables["MAP2D"]
for j in range(NUM_NY_START, NUM_NY_END):
    for i in range(NUM_NX_START, NUM_NX_END):
        VAR_IN2D[j,i] = NUM_MPI_RANK

MPI_SET.MPI_MESSAGE("End of MAP2D     ")
VAR_IN3D = NC_OUT.variables["MAP3D"]
for j in range(NUM_NY_START, NUM_NY_END):
    for i in range(NUM_NX_START, NUM_NX_END):
        for z in range(numNT):
            VAR_IN3D[z,j,i] = NUM_MPI_RANK * z 

MPI_SET.MPI_MESSAGE("End of MAP3D     ")


NC_OUT.close()
