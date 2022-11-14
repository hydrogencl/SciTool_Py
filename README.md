SciTool_Py
==========

All the files are protect by MIT Lincen, Use it, distribute it, and NEVER copyright them. 

The main works are here: 

+ The algorithm to read, write, and deal with grids and arrays in geoscience. 
+ The handy tools for programing.
+ Some handy tools for plotting. 
+ Handy tools for WRF


## LIGHTWEIGHT-ENSEMBLE-FRAMEWORK (LEF)

+ The framwork that using Slurm CLI to redistribute the simulation for 
  ensemble simulations. 
+ Based on OO, this framework can be working with only few commandlines:
    + SlurmController
        `strProject`      : The name of project
        `strPartition`    : The partition name matching the setup for Slurm
        `numCoresPerNode` : The Cores per Node matching the setup for Slurm
        `strJobname`      : The Jobname that used as the stdout/stderr matching Slurm
        `strRootdir`      : The root of application folder, e.g. `WRF`
        `ifServerLog`     : ON/OFF for storing the server log of the LEF

        + `InitEnsemble`  : To init the ensemble simulation, to indicate that:
            + `members`     : How many members for the simulation
            + `UsingNodes`  : How many nodes for one ensemble member ( can be a fraction)

        + `CreateMembers` : To create the members with some prefix
            + 

        + `FileControl`   : To manipulate the files for all the ensemble members

        + `RunMembers`    : To execute the simulation of the ensemble members, 
            + The executor need to be indicated

        + `CheckMembersWRF` : A simple tool to check if the WRF of each member is still alive

## GRIDINFORMATER 

Will deal with gridcells. 

### progress bar ###
Based on the work from stackflow, the progress bar can be used in `Jupyter-Notebook` and `Terminal` and Python IDE. 

## Taylor's Diagram - YS

This is a simple example of Taylor's diagram, which requires only `matplotlib`. 
The spacing of the correlation axis is depending on a sin fixing. 

## Gridtrans

This is used to create and read the output/input files for these format:

- `.pfb`
- `.csv`
- `.sa`
- `.dat`  (seperater is `,`)
- `.xyz`

With `matplotlib` user can use it to quickly plot the array files. 
One day this will be integrated into GRIDINFORMER

## Bib_Corrected

A small tool to fix the .bib file, to remove all the round bracket that can not be read by `Gummi`. 
Before running the bibliothek function of Gummi, run this script to let `Gummi` read the .bib file. 


