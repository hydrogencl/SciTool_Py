### 
#  This is the example for using MPI_TOOLS
#     1) Creating the geometry for mpi rank for the input array
#     2) Calculating some number based on the mpi_rank.
###

import math, sys,os, time, re
import numpy as NP
from GRIDINFORMER import MPI_TOOLS as MT

numNX = 16
numNY = 5 

# Start to setup the MPI from the MPI_TOOLS
# Firstly, set setup the NX and NY
MPI_SET = MT(NUM_NX_END=numNX, NUM_NY_END=numNY)

# Update the cores topology, NX * NY = MPI_SIZE
# NY will be calculated automatically
MPI_SET.UPDATE_NX_CORES(MPI_SET.NUM_MPI_SIZE/2)

# Designing the CPU Geometry for each Size. 
# My script will try to distribute the array index evenly. 

MPI_SET.CPU_GEOMETRY_2D()
ARR_RANK_DESIGN = MPI_SET.ARR_RANK_DESIGN

# This is for the example for the 3D array. 
numNT = MPI_SET.NUM_MPI_SIZE

# The message from each Rank
MPI_SET.MPI_MESSAGE("Starting to work!")

# The final design of the array index (for x and y)
# Usually using only x is fine if NX_CORES = NUM_MPI_SIZE
NUM_NX_START = ARR_RANK_DESIGN[MPI_SET.NUM_MPI_RANK]["NX_START"]
NUM_NX_END   = ARR_RANK_DESIGN[MPI_SET.NUM_MPI_RANK]["NX_END"]
NUM_NY_START = ARR_RANK_DESIGN[MPI_SET.NUM_MPI_RANK]["NY_START"]
NUM_NY_END   = ARR_RANK_DESIGN[MPI_SET.NUM_MPI_RANK]["NY_END"]

VAR_IN2D = [[ 0 for x in range(numNX)    ] for y in range(numNY) ]
for j in range(NUM_NY_START, NUM_NY_END):
    for i in range(NUM_NX_START, NUM_NX_END):
        VAR_IN2D[j][i] = MPI_SET.NUM_MPI_RANK + 1

time.sleep(MPI_SET.NUM_MPI_RANK*0.1)
MPI_SET.MPI_MESSAGE("Printing out:  Rank ( +1 to identify the first rank (index=0) )    ")
print(NP.array(VAR_IN2D))
MPI_SET.MPI_MESSAGE("End of MAP2D     ")

#VAR_IN3D = [[[ 0 for x in range(numNX)    ] for y in range(numNY) ] for t in range(numNT)   ] 
#for j in range(NUM_NY_START, NUM_NY_END):
#    for i in range(NUM_NX_START, NUM_NX_END):
#        for z in range(numNT):
#            VAR_IN3D[z][j][i] = MPI_SET.NUM_MPI_RANK + 1 + z 

#time.sleep(MPI_SET.NUM_MPI_RANK*0.1)
#MPI_SET.MPI_MESSAGE("Printing out: Z + Rank + 1    ")
#print(NP.array(VAR_IN3D))
#MPI_SET.MPI_MESSAGE("End of MAP3D     ")

