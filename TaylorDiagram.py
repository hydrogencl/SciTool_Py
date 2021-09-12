import mpl_toolkits.axisartist.floating_axes as floating_axes
import numpy as NP
import mpl_toolkits.axisartist.angle_helper as angle_helper
from   matplotlib.projections import PolarAxes
from   matplotlib.transforms import Affine2D
from   mpl_toolkits.axisartist.grid_finder import (FixedLocator, MaxNLocator,
                                                   DictFormatter)
import matplotlib.pyplot as PLT
import math, re

"""
 REF: 
 Taylor, K. E. (2001), Summarizing multiple aspects of model performance in a single diagram,
 J. Geophys. Res., 106(D7), 7183–7192, http://dx.doi.org/10.1029/2000JD900719.1

"""


class TaylorsDiagram:
    Sc  = lambda X : math.cos(X*NP.pi*0.50)
    NUM_PI      = NP.pi
    angle_ticks = [(1.0 *.5*NUM_PI, r"0   "),
                       (0.90*.5*NUM_PI, r"0.2 "),
                       (0.80*.5*NUM_PI, r"0.3 "),
                       (0.70*.5*NUM_PI, r"0.4 "),
                       (0.60*.5*NUM_PI, r"0.5 "),
                       (0.50*.5*NUM_PI, r"0.6 "),
                       (0.40*.5*NUM_PI, r"0.7 "),
                       (0.30*.5*NUM_PI, r"0.8 "),
                       (0.20*.5*NUM_PI, r"0.9 "),
                       (0.10*.5*NUM_PI, r"0.95"),
                       (0.01*.5*NUM_PI, r"0.99"),
                       (0.00*.5*NUM_PI, r"    ")]

    angle_grids = [(Sc(0.0 )*.5*NUM_PI, r"0   "),
                   (Sc(0.20)*.5*NUM_PI, r"0.2 "),
                   (Sc(0.30)*.5*NUM_PI, r"0.3 "),
                   (Sc(0.40)*.5*NUM_PI, r"0.4 "),
                   (Sc(0.50)*.5*NUM_PI, r"0.5 "),
                   (Sc(0.60)*.5*NUM_PI, r"0.6 "),
                   (Sc(0.70)*.5*NUM_PI, r"0.7 "),
                   (Sc(0.80)*.5*NUM_PI, r"0.8 "),
                   (Sc(0.90)*.5*NUM_PI, r"0.9 "),
                   (Sc(0.95)*.5*NUM_PI, r"0.95"),
                   (Sc(1.00)*.5*NUM_PI, r"    ")]

    angle_ticks = [(Sc(0.0 )*.5*NUM_PI, r"0   "),
                   (Sc(0.20)*.5*NUM_PI, r"0.2 "),
                   (Sc(0.30)*.5*NUM_PI, r"0.3 "),
                   (Sc(0.40)*.5*NUM_PI, r"0.4 "),
                   (Sc(0.50)*.5*NUM_PI, r"0.5 "),
                   (Sc(0.60)*.5*NUM_PI, r"0.6 "),
                   (Sc(0.70)*.5*NUM_PI, r"0.7 "),
                   (Sc(0.80)*.5*NUM_PI, r"0.8 "),
                   (Sc(0.90)*.5*NUM_PI, r"0.9 "),
                   (Sc(0.91)*.5*NUM_PI, r"    "),
                   (Sc(0.92)*.5*NUM_PI, r"    "),
                   (Sc(0.93)*.5*NUM_PI, r"    "),
                   (Sc(0.94)*.5*NUM_PI, r"    "),
                   (Sc(0.95)*.5*NUM_PI, r"0.95"),
                   (Sc(0.96)*.5*NUM_PI, r"    "),
                   (Sc(0.97)*.5*NUM_PI, r"    "),
                   (Sc(0.98)*.5*NUM_PI, r"    "),
                   (Sc(0.99)*.5*NUM_PI, r"0.99"),
                   (Sc(1.00)*.5*NUM_PI, r"    ")]

    def __init__(self, FIG_IN=None, arr_pos=[1,1,1], num_refstd=1.0, num_mag=1.5, if_rel=False):
        """
        Adopting from the examples on Matplotlib.org
        """
        
        tr = PolarAxes.PolarTransform()
        grid_locator1   = FixedLocator([v for v, s in self.angle_grids])
        tick_formatter1 = DictFormatter(dict(self.angle_ticks))
        grid_locator2   = MaxNLocator(2)
        
        self.num_mag    = num_mag
        self.num_refstd = num_refstd
        self.num_max    = self.num_refstd * self.num_mag
    
        if if_rel:
            grid_locator2 = FixedLocator([0.0,0.5,1.0,1.5,2.0])
        else:
            grid_locator2 = FixedLocator(NP.arange(0.0, self.num_max+0.1, 0.5))
            
    
        grid_helper = floating_axes.GridHelperCurveLinear(
            tr, extremes=(0.5*self.NUM_PI,0, self.num_max, 0),
            grid_locator1=grid_locator1,
            grid_locator2=grid_locator2,
            tick_formatter1=tick_formatter1,
            tick_formatter2=None           )
        if FIG_IN == None:
            FIG_IN = PLT.figure()

        ax1 = floating_axes.FloatingSubplot(FIG_IN, arr_pos[0],arr_pos[1],arr_pos[2], grid_helper=grid_helper)
         
        FIG_IN.add_subplot( ax1 )
        ax1.axis["bottom"].set_axis_direction("bottom")  # "Angle axis"
        ax1.axis["bottom"].major_ticklabels.set_axis_direction("top")
        ax1.axis["bottom"].major_ticklabels.set_pad(10)
        ax1.axis["bottom"].minor_ticklabels.set_axis_direction("top")
        ax1.axis["bottom"].minor_ticklabels.set_pad(10)
        ax1.axis["bottom"].label.set_axis_direction("top")
        ax1.axis["bottom"].label.set_text("Correlation")
        ax1.axis["left"].set_axis_direction("top") # "X axis"
        if if_rel:
          ax1.axis["left"].label.set_text("Relative Std")
        else:
          ax1.axis["left"].label.set_text("Standard deviation")
    
        ax1.axis["right"].set_axis_direction("bottom")   # "Y axis"
        ax1.axis["right"].toggle(ticklabels=True)
        ax1.axis["right"].major_ticklabels.set_axis_direction("bottom")
        ax1.axis["right"].major_ticklabels.set_pad(10)
        if if_rel:
          ax1.axis["right"].label.set_text("Relative Std")
        else:
          ax1.axis["right"].label.set_text("Standard deviation")
    
        # create a parasite axes whose transData in RA, cz
        aux_ax = ax1.get_aux_axes(tr)
        ax1.grid(linestyle=":", zorder=1)
        aux_ax.patch = ax1.patch  # for aux_ax to have a clip path as in ax
        ax1.patch.zorder = 0.9    # but this has a side effect that the patch is
        aux_ax.tick_params(axis='both', which='both', pad=30)
        ax1.tick_params(axis='both', which='both', pad=30)
        self.ax  = ax1
        self.aux = aux_ax


    def add_sample(self, NUM_Cor, NUM_Std, if_rel=False, STR_TEXT="", *args, **kwargs ):
        Sc  = lambda X : math.cos(X*NP.pi*0.50)
        NUM_theta = Sc(NUM_Cor) * 0.5 * self.NUM_PI 
        #NUM_theta = (NUM_Cor) * 0.5 * self.NUM_PI 
        NUM_Std  =  NUM_Std/self.num_refstd 
        self.aux.scatter(NUM_theta, NUM_Std, zorder=3, *args, **kwargs )
        if not STR_TEXT == "": 
            self.aux.text (NUM_theta*1.05, NUM_Std*1.05, STR_TEXT, fontsize=10)
    
    def add_contour(self, num_intvl=0.25, num_cnts=0, alpha=0.5, color='#b4b4b4', if_rel=False):
        if num_cnts == 0:
            num_cnts = int(self.num_mag / num_intvl )
        arr_psi = NP.arange(0, self.NUM_PI+0.01, 0.05)
        arr_r   = [[] for n in range(num_cnts)]
        arr_theta   = [[] for n in range(num_cnts)]
        for n in range(num_cnts):
            for psi in arr_psi:
                r     = num_intvl * (n+1) * self.num_refstd
                R     = (self.num_refstd**2 - 2 * self.num_refstd * r * math.cos( psi) + r **2) ** 0.5
                if R <=0:
                    theta = 0
                else:
                    theta =  math.asin( r * math.sin(psi) / R )
                arr_r[n]    .append( R     )
                arr_theta[n].append( theta )
            self.aux.plot(arr_theta[n], arr_r[n], c=color, ls='--', alpha=alpha, zorder=2)


##########################################################

if __name__ == "main":
    print("This is an example")

    FIG = PLT.figure(figsize=(5,5))
    TDYY = TDY.TaylorsDiagram(FIG_IN=FIG, if_rel=True)

    TDYY.add_sample(0.9, 0.78,  marker="^", STR_TEXT="0.9, 0.78")
    TDYY.add_sample(0.95, 0.87, c="r", STR_TEXT="0.95, 0.87")
    TDYY.add_sample(0.7, 0.5, c="g", STR_TEXT="0.7, 0.5")
    TDYY.add_cntour()

    




