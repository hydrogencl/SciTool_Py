### 
#  This is the example for using MPI_TOOLS
#     1) Creating the geometry for mpi rank for the input array
#     2) Calculating some number based on the mpi_rank.
###

import math, sys,os, time, re
from mpi4py import MPI
import numpy as NP
from GRIDINFORMER import MPI_TOOLS as MT

IF_PARALLEL = True

# MPI configuration
# Initialing the MPI Commworld
comm     = MPI.COMM_WORLD
# Get the rank number (aka index of the CPU)
NUM_MPI_RANK = comm.Get_rank()
# Get the total size of CPUs
NUM_MPI_SIZE = comm.Get_size()

numNX = 12
numNY = 12
numNT = NUM_MPI_SIZE

# Design Cores, can be the same as NUM_MPI_SIZE
# numCores   = 48
numCores   = NUM_MPI_SIZE

numCoresNX =  int(NUM_MPI_SIZE / 2)

MPI_SET = MT(MPI_SIZE=NUM_MPI_SIZE, MPI_RANK = NUM_MPI_RANK,\
             NUM_NX_END=numNX ,\
             NUM_NY_END=numNY ,\
             NUM_NX_CORES=numCoresNX)

MPI_SET.CPU_GEOMETRY_2D()
ARR_RANK_DESIGN = MPI_SET.ARR_RANK_DESIGN

MPI_SET.MPI_MESSAGE("Starting to work!")

NUM_NX_START = ARR_RANK_DESIGN[NUM_MPI_RANK]["NX_START"]
NUM_NX_END   = ARR_RANK_DESIGN[NUM_MPI_RANK]["NX_END"]
NUM_NY_START = ARR_RANK_DESIGN[NUM_MPI_RANK]["NY_START"]
NUM_NY_END   = ARR_RANK_DESIGN[NUM_MPI_RANK]["NY_END"]

VAR_IN2D = [[ 0 for x in range(numNX)    ] for y in range(numNY) ]
for j in range(NUM_NY_START, NUM_NY_END):
    for i in range(NUM_NX_START, NUM_NX_END):
        VAR_IN2D[j][i] = NUM_MPI_RANK + 1

time.sleep(z*0.1)
MPI_SET.MPI_MESSAGE("Printing out:  Rank (start from 1)    ")
print(NP.array(VAR_IN2D))
MPI_SET.MPI_MESSAGE("End of MAP2D     ")

VAR_IN3D = [[[ 0 for x in range(numNX)    ] for y in range(numNY) ] for t in range(numNT)   ] 
for j in range(NUM_NY_START, NUM_NY_END):
    for i in range(NUM_NX_START, NUM_NX_END):
        for z in range(numNT):
            VAR_IN3D[z][j][i] = NUM_MPI_RANK + 1 + z 

time.sleep(z*0.1)
MPI_SET.MPI_MESSAGE("Printing out: Z + Rank + 1    ")
print(NP.array(VAR_IN3D))
MPI_SET.MPI_MESSAGE("End of MAP3D     ")

