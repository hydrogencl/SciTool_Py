SciTool_Py
==========

All the files are protect by LGPL. Use it, distribute it, and NEVER copyright them. 

The main works are here: 

+ The algorithm to read, write, and deal with grids and arrays in geoscience. 
+ The handy tools for programing.
+ Some handy tools for plotting. 
+ Handy tools for WRF

### GRIDINFORMATER 

Will deal with gridcells. 

### progress bar ###
Based on the work from stackflow, the progress bar can be used in `Jupyter` and `Terminal` and Python IDE. 

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


