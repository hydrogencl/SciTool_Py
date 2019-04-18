#!/usr/bin/python
#Purpose: 1. To convert the grid data of following format to each other:
#            .sa .xyz .asc .pfb
#            (simple ascii, xyz, ascii grid, parflow binary)
#         2. To simply plot the figure from the single output file.
#ChangeLog: 20150429: Changing the reading method of read pfb. Make it read faster. 
#

import argparse
import re,math
import struct
plot_funk=1

try:
    import numpy as numpy
except ImportError: 
    print("You no install numpy? ")
    print("You could not use the plot function of this script until you install Numpy")
    plot_funk=0
try:
    from mpl_toolkits.mplot3d import axes3d
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from pylab import *
    import matplotlib.animation as animation
except ImportError: 
    print("You could not use the plot function of this script until you install Matplotlib")
    plot_funk=0

# -----
# Strip the blank while readline
# -----
#Read and Strip the Blanks
#Default re_pattern : \s
#    if len(re_pattern) == 0:
#        str_arr=re.split('\s',fdata.readline())
#    elif len(re_pattern) == 1:
#        str_arr=re.split(re_pattern,fdata.readline())
#    else:
#        print("Wrong for re_pattern, Bug#001")
#    new_arr=[]
#    for s in str_arr:
#        if s =="":
#            pass
#        else:
 #            new_arr.append(s)
#    return new_arr
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


def tryopen(sourcefile,ag):
    try:
        opf=open(sourcefile,ag)
        return opf
    except :
        print("No such file.")
        return "error"

def checkformat(sourcefile):
    fmtt=re.split('\.',sourcefile)
    fmt=fmtt[len(fmtt)-1]
    return fmt

# -----
# Read .pfb
# -----

# -----
# Read .csv 
# -----

def readcsv(sourcefile, str_null="noData", num_null=-999.999):
    opf    = tryopen(sourcefile,'r')
    opfchk = tryopen(sourcefile,'r')
    print("reading source file {0:s}".format(sourcefile))
    chk_lines = opfchk.readlines()
    num_totallines = len(chk_lines)
    ncols = 0
    num_notnum = 0
    for n in range(num_totallines):
        line_in = chk_lines[n]
    #print line_in
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
        result_arr=[[num_null for j in range(nrows)] for i in range(ncols)]
        result_arr_text=[]
        num_pass = 0
        for j in range(0,num_totallines):
            # chk if comment
            line_in = opf.readline()
            c_first = re.findall(".",line_in.strip())[0]
            if c_first == "#":
               result_arr_text.append(line_in)
                num_pass += 1
            else:
                arr_in = re.split(",",line_in.strip())
                for i in range(ncols):
                    chk_val = arr_in[i]
                    if chk_val == str_null: 
                        result_arr[i][j-num_pass] = num_null
                        print j,i,chk_val
                    else:
                        result_arr[i][j-num_pass] = float(chk_val)
    return result_arr,result_arr_text

def readpfb(sourcefile):
    opf = tryopen(sourcefile,'rb')
    print("reading source file {0:s}".format(sourcefile))
    
    t1=struct.unpack('>ddd',opf.read(24))
    tn=struct.unpack('>iii',opf.read(12))
    td=struct.unpack('>ddd',opf.read(24))
    tns=struct.unpack('>i',opf.read(4))
    x1,y1,z1=t1
    nx,ny,nz=tn
    dx,dy,dz=td
    ns=tns[0]
    result_arr=list([[[0.0 for i in range(nx)] for j in range(ny)] for k in range(nz)])
    for isub in range(0,ns):
        ix,iy,iz,nnx,nny,nnz,rx,ry,rz=struct.unpack('>9i',opf.read(36))
        tmp_total = nnx * nny * nnz
        tvalue = struct.unpack('>{0:d}d'.format(tmp_total), opf.read(8*tmp_total))
        for k in range(nnz):
            for j in range(nny):
                for i in range(nnx):
                    result_arr[k+iz][j+iy][i+ix]=tvalue[ k*(nny*nnx) + j*nnx + i  ]
             
    opf.close()
    print("Completed reading pfb format from {0}".format(sourcefile))
    return result_arr,nx,ny,nz,dx,dy,dz


# -----
# Read .sa
# -----
def readsa(sourcefile):
    print("reading source file {0:s}".format(sourcefile))
    result_arr=[]
        opf = tryopen(sourcefile,'r')
    headt=re.split('\s',opf.readline().strip())
    head=stripblnk(headt,'int')
    nx=int(head[0])
    ny=int(head[1])
    nz=int(head[2])
    for j in range(0,ny):
        tmp=[]
        for i in range(0,nx):
                ans=re.split('\s',opf.readline().strip())
            tmp.append(float(ans[0]))
            result_arr.append(tmp)
        print("Completed reading sa format from {0}".format(sourcefile))
    return result_arr,nx,ny,nz

# -----
# Read .asc (ascii grid from ESRI)
# -----

def readascgrid(sourcefile):
    opf = tryopen(sourcefile,'r')
    print("reading source file {0:s}".format(sourcefile))
    ncols=int(re.split('\s',opf.readline().strip())[1])
    nrows=int(re.split('\s',opf.readline().strip())[1])
    xllcorner=float(re.split('\s',opf.readline().strip())[1])
    yllcorner=float(re.split('\s',opf.readline().strip())[1])
    cellsize=float(re.split('\s',opf.readline().strip())[1])
    nodata_v=float(re.split('\s',opf.readline().strip())[1])

    result_arr=[]
    for j in range(0,nrows):
        valuet=[]
        arr_valuet=re.split('\s',opf.readline().strip())
        for i in range(0,ncols):
            valuet.append(float(arr_valuet[i]))    
        result_arr.append(valuet)

        print("Completed reading ascii-grid format from {0}".format(sourcefile))
    return result_arr,ncols,nrows,xllcorner,yllcorner,cellsize,nodata_v

# -----
# Read .xyz
# -----

def readxyz(sourcefile,*nxny):
    chk = tryopen(sourcefile,'r')
    print("reading source file {0:s}".format(sourcefile))
    opf = tryopen(sourcefile,'r')
    result_arr=[]
    if len(nxny) == 2:
        print("Specific nx and ny is indicated.")
        ncol=nxny[0]
        nrow=nxny[1]

    elif len(nxny) == 0:
        
            print("Checking the nx and ny")
            count=1
            k=chk.readlines()
            start_col=float(re.split('\s',k[0].strip())[0])
            check = -9999
            while check != 0.0:
                    check_col=float(re.split('\s',k[count].strip())[0])
                    check=check_col-start_col
                    if check ==0.0 :
                            ncol=count
                        
                    else:
                           count=count+1
            nrow=len(k)/ncol

    for j in range(0,nrow):
        valuex=[]
        for i in range(0,ncol):
            tline=re.split('\s',opf.readline().strip())
            line=stripblnk(tline,'float')
            valuex.append(line[2])
            if i == 0 and j == 0:
                xll=float(line[0])
                print('xll: {0}'.format(xll))
            if i == 1 and j == 0:
                refx=float(line[0])
                print('refx: {0}'.format(refx))
            if i == 0 and j == 1:
                refy=float(line[1])
                print('refy: {0}'.format(refy))
            if i == ncol-1 and j == nrow-1:
                yll=float(line[1])
                print('yll: {0}'.format(yll))
        result_arr.append(valuex)

    #dx=xll-refx
    #dy=refy-yll
        print("Completed reading ascii-grid format from {0}".format(sourcefile))
    #print result_arr
    return result_arr,ncol,nrow,xll,yll,#dx,dy

# ----- 
# Read Custom 1 col raster
# -----

def read1col(sourcefile,nx,ny,nz,*skiplines):
    print("reading source file {0:s}".format(sourcefile))
    opf = tryopen(sourcefile,'r')
    result_arr=[]
    if len(skiplines) == 0:
        pass
    else:
        for m in range(0,skipelines[0]):
            opf.readline()

    for j in range(0,ny):
        tmp=[]
        for i in range(0,nx):
                ans=re.split('\s',opf.readline().strip())
            tmp.append(float(ans[0]))
            result_arr.append(tmp)
        print("Completed reading sa format from {0}".format(sourcefile))
    return result_arr,nx,ny,nz

# -----
# Read Custom 2d grid raster
# -----
def read2d(sourcefile,nx,ny,num_typ,*skiplines):
    print("reading source file {0:s}".format(sourcefile))
    opf = tryopen(sourcefile,'r')
    result_arr=[]
    if len(skiplines) == 0:
        pass
    else:
        for m in range(0,skiplines[0]):
            opf.readline()
    
    for j in range(0,ny):
        ans=re.split('\s',opf.readline().strip())
        t_arr=stripblnk(ans,num_typ)
        result_arr.append(t_arr)

    return result_arr,nx,ny


# -----
# Write from 2D array to 1D array
# -----

# -----
# Write from 2D array to 2D array
# -----

# -----
# Write .xyz
# -----

def writexyz(write_file_name,input_arr,ncols,nrows,xllco,yllco,cellsize,nodata_value):
    wtf=open(write_file_name)
    for j in range(ncows):
        for i in range(ncols):
            wtf.write("{0} {1} {2}\n".format(xllco+i*cellsize,yllco+nrows*cellsize-i*cellsize,value[j][i]))
    wtf.close()


# -----
# Write .pfb
# -----
def writepfb(write_file_name,input_arr,nx,ny,nz,dx=0,dy=0,dz=0,x=0,y=0,z=0,ns=1,nodata_value=-999.999):
    wtf=open(write_file_name,"w")
    wtf.write(struct.pack('>3d',x,y,z))
    wtf.write(struct.pack('>3i',nx,ny,nz))
    wtf.write(struct.pack('>3d',dx,dy,dz))
    wtf.write(struct.pack('>1i',ns))

    for isub in range(0,ns):
        iz,iy,ix=int(z),int(y),int(x)
        nnz,nny,nnx=int(nz),int(ny),int(nx)
        wtf.write(struct.pack('>3i',0,0,0))
        wtf.write(struct.pack('>3i',nx,ny,nz))
        wtf.write(struct.pack('>3i',dx,dy,dz))
        for i in range(iz,iz+nnz):
            for j in range(iy,iy+nny):
                for k in range(ix,ix+nnx):
                    wtf.write(struct.pack('>d',input_arr[i][j][k]))
    wtf.close()


# -----
# Write .sa
# -----

def writesa(write_file_name,input_arr,nz,ny,nx):
    wtf=open(write_file_name,"w")
    wtf.write('{0} {1} {2}\n'.format(nx,ny,nz))
    for k in range(0,nz):
        for j in range(0,ny):
            for i in range(0,nx):
                wtf.write(str(input_arr[j][i]) + '\n' )
    wtf.close()
# -----
# Write .asc 
# -----
def writeasc(write_file_name,input_arr,ncols,nrows,xllco,yllco,cellsize,nodata_v):
    wtf=open(write_file_name,'w')
    wtf.write("ncols         {0}\n".format(ncols))
    wtf.write("nrows         {0}\n".format(nrows))
    wtf.write("xllcorner     {0}\n".format(xllco))
    wtf.write("yllcorner     {0}\n".format(yllco))
    wtf.write("cellsize      {0}\n".format(cellsize))
    wtf.write("NODATA_value  {0}\n".format(nodata_v))
    for j in range(0,nrows):
            for i in range(0,ncols):
                for x in input_arr[j]:
                        wtf.write("{0} ".format(x))
                wtf.write("\n")
    wtf.close()
    
# -----
# Plot in default_im 
# -----
def im_subplots(title,ax,figx,ncols,nrows,array,*inver):
    figx=plt.figure()
    ax=figx.add_subplot(111)
    cax=ax.imshow(array)
    #cax=ax.imshow('array[0][layer],vmin=%f,vmax=%f' % (float(vmin),float(vmax)))
    ax.set_title(title)
    ax.set_xlim((-0.5,ncols-0.5))
    if len(inver) == 0:
        ax.set_ylim((-0.5,nrows-0.5))
    else: 
        ax.set_ylim((nrows-0.5,-0.5))
        cbar= figx.colorbar(cax,format='%.5f')


# -----
# Plot .pfb in im
# -----
def pfb_im_subplots(title,ax,figx,array,ncols,nrows,layer,*inver):
    figx=plt.figure()
    #ax=figx.add_subplot('%d%d%d'%(col,row,pos))
    ax=figx.add_subplot(111)
    v_max = numpy.amax(array[layer])
    v_min = max(numpy.amin(array[layer]),-15)
    if v_max == v_min: v_max=v_min+0.5
    cax=ax.imshow(array[layer],vmin=v_min,vmax=v_max,interpolation='nearest')
    ax.set_title('File:{0:s}, Layer:{1:d}'.format(title,layer))
    ax.set_xlim((-0.5,ncols-0.5))
    if inver == True:
        ax.set_ylim((nrows-0.5,-0.5))
    else:
        ax.set_ylim((-0.5,nrows-0.5))
    cbar= figx.colorbar(cax,format='%.5f')

# -----
# Chunk the value from pfb file
# -----

def pfb_chunk3d(pfba,fixplane,layer,ncols,nrows,nlayers):
    #Must call the array[0] from readpfb
    fp=fixplane
    nx=ncols
    ny=nrows
    nz=nlayers
    if fp == "x":
        print("Chunk from X plane, layer: {0}".format(layer))
        value=zeros((nz,ny))
        for j in range(nz):
            for i in range(ny):
                value[j][i]=pfba[j][i][layer]
        return value,ny,nz

    elif fp =="y":
        print("Chunk from Y plane, layer: {0}".format(layer))
        value=zeros((nz,nx))
        for j in range(nz):
            for i in range(nx):
                value[j][i]=pfba[j][layer][i]
        return value,nx,nz

    elif fp =="z":
        print("Chunk from Z plane, layer: {0}".format(layer))
        value=zeros((ny,nx))
        for j in range(ny):
            for i in range(nx):
                value[j][i]=pfba[layer][j][i]
        return value,nx,ny

    else:
        print("Wrong fix plane, \"x\", \"y\" or \"z\" only")
        print("Get chunk at {0}-plane,{1} layer".format(fixplane,layer))

  
  # -----
  # Persing arguments 
  # -----

if __name__=="__main__":
  
      parser = argparse.ArgumentParser(description="""Converting grid from format to format and plot it.
      
      
      ***Following package should be install for plot: matplotlib and numpy  ***
      ***NO WARRENTY FOR ANY PURPOSE.                                        ***
      *** Free to use, please report bug to yslu@uni-bonn.de                 ***""")
      parser.add_argument('-s','--source-file',nargs='*', dest='sourcefile', help='Input ONE source FILE NAME and/or FORMAT. Will automatically detect the sub-filename. Support format: sa(simple ascii), pfb(Parflow Binary), xyz (xyz grid),asc (ascii grid),1D or 2D array')
      
      parser.add_argument('-o','--output-file',nargs='*', dest='outputfile', help='Input ONE output FILE NAME and/or FORMAT. Will automatically detect the sub-filename')
      parser.add_argument('-p','--parameters',nargs='*', dest='parameter',default=[0], help='Support Parameters (in order) for output file: NX,NY,NZ,XLLcorner,YLLcorner,Cell Size. The unnecessary parameter could be empty. The default parameter is: NX=1, NY=1, NZ= 1, XLLcorner=YLLcorner=0, Cell Size=1, No data value=-9999')
      parser.add_argument('-ufp','--uf-parameters',action='store_true', help='User-Friendly parameter input interface')
      
      parser.add_argument('-im','--im_plot', action='store_true', help='Plot the grid file. The pfb will be plot on the surface.')
      parser.add_argument('-imch3d','--im_plot_chunk_3d', nargs='*', dest='para_ch3d' , help='Plot the pfb in the speciic plane (x,y,z). Indicate the parameters of chunk:[plane] [layer] (all must set).')
      parser.add_argument('-ims','--im_plot_slider', action='store_true', help='Plot the grid file. The pfb will be plot on the surface.')
      
      args = parser.parse_args()
        
      #-----
      # READ THE INPUT
      #-----
      
      
      # -----
      # Default parameters if not indicate:
      # -----
      s_nlays=1
      s_xllcorner = 0.00
      s_xllcorner = 0.00
      s_cellsize =  1.00
      s_nodata_v = -9999
      
      o_nlays=1
      o_xllcorner = 0.00
      o_xllcorner = 0.00
      o_cellsize =  1.00
      o_nodata_v = -9999
      
      if args.sourcefile == None:
      
          print("IMPORTING GRIDTRANS, NO WARRENTY FOR ANY PURPOSE AND USAGE")
      
      elif len(args.sourcefile) > 0:
      
      #---
      # Read the source file name
      #---
          s=args.sourcefile[0]
          if len(args.sourcefile)==1:
              sf=checkformat(s)
          elif len(args.sourcefile) == 2:
              sf=args.sourcefile[1]
          else:
              print("TOO MUCH INPUT, please check!")
              sf=""
      
      #---
      # Read the source and store the values
      #---
      
          if sf == 'sa':
              value=readsa(s)
              nlays=value[3]
          elif sf == 'pfb':
              value=readpfb(s)
              nlays=value[3]
              dz=value[6]
          elif sf == 'asc':
              value=readascgrid(s)
              s_xllcorner = value[3]
              s_xllcorner = value[4]
              s_cellsize =  value[5]
          elif sf == 'xyz':
      #        print len(args.parameter)
              if len(args.parameter) >= 2:
                  value=readxyz(s,args.parameter[0],args.parameter[1])
              else:
                  value=readxyz(s)
              s_xllcorner=value[3]
              s_yllcorner=value[4]
              dn = 0 #value[5]-value[6]
              if dn == 0:
                  s_cellsize=dn
              else:
                  s_cellsize=(value[5]+value[6])/2.0
                  print("Warning! DX not equal to DY.")
          elif sf == '2d':
              value=read2d(s,ncol,nrow,num_typ)
          elif sf == '1d':
              value=read1col(s)
          else:
              print("No recognized file or format")
              value=[0]
                  
          result_arr=value[0]
          ncols=value[1]
          nrows=value[2]
      
      #-----
      # READ PARAMETERS
      #-----
      
      if len(args.parameter) > 3:
          nrows=int(args.parameter[0])
          ncols=int(args.parameter[1])
          nlays=int(args.parameter[2])
          o_xllcorner=args.parameter[3]
          o_yllcorner=args.parameter[4]
          o_cellsize=args.parameter[5]
          o_nodata_v=args.parameter[6]
      elif len(args.parameter) == 3:
          nrows=int(args.parameter[0])
          ncols=int(args.parameter[1])
          nlays=int(args.parameter[2])
          
      
      #-----
      # MAKE THE OUTPUT
      #-----
      
      #args.outputfile
      
      #-----
      # MAKE THE PLOT
      #-----
      if args.im_plot == True:
          if sf=='pfb':
              layer=nlays-1
              pfb_im_subplots("{0:s}".format(args.sourcefile[0]),"ax1","fig2",result_arr,ncols,nrows,layer)
      
          else:
              im_subplots("{0:s}".format(args.sourcefile[0]),"ax1","fig1",ncols,nrows,result_arr)
          plt.grid(True)
          plt.show()
      
      
      if args.para_ch3d != None:
          if len(args.para_ch3d) !=2:
              print("Please input the parameters: plane \"x\",\"y\",\"z\", and the layer")
          else:
              if sf=='pfb':
                  chunk_arr=pfb_chunk3d(result_arr,str(args.para_ch3d[0]),int(args.para_ch3d[1]),ncols,nrows,nlays)
                  im_subplots("{}".format(args.sourcefile[0]),"ax1","fig1",chunk_arr[1],chunk_arr[2],chunk_arr[0])
                  plt.grid(True)
                  plt.show()
              else:
                  print("Sorry, pfb only")
