from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes
import numpy as np
import mpl_toolkits.axisartist.angle_helper as angle_helper
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist.grid_finder import (FixedLocator, MaxNLocator,
                                                 DictFormatter)
import matplotlib.pyplot as plt
import math, re
def Taylor_Diagram(fig, arr_pos, num_refstd, num_mag=1.5, if_rel=False):
    """
    Adopting from the examples on Matplotlib.org
    """
    tr = PolarAxes.PolarTransform()

    pi = np.pi
    angle_ticks = [(0.01*.5*pi, r"0.99"),
                   (0.1*.5*pi, r"0.9"),
                   (0.2*.5*pi, r"0.8"),
                   (0.3*.5*pi, r"0.7"),
                   (0.4*.5*pi, r"0.6"),
                   (0.5*.5*pi, r"0.5"),
                   (0.6*.5*pi, r"0.4"),
                   (0.7*.5*pi, r"0.3"),
                   (0.8*.5*pi, r"0.2"),
                   (0.9*.5*pi, r"0.1"),
                   (1.0*.5*pi, r"0.0")]
    grid_locator1 = FixedLocator([v for v, s in angle_ticks])
    tick_formatter1 = DictFormatter(dict(angle_ticks))

    grid_locator2 = MaxNLocator(2)
    if if_rel:
      grid_locator2 = FixedLocator([0.0,0.5,1.0,1.5,2.0])
      num_refstd = float(num_refstd/num_refstd)
    num_max = num_refstd * num_mag

    grid_helper = floating_axes.GridHelperCurveLinear(
        tr, extremes=(0.5*pi,0, num_max, 0),
        grid_locator1=grid_locator1,
        grid_locator2=grid_locator2,
        tick_formatter1=tick_formatter1,
        tick_formatter2=None)

    ax1 = floating_axes.FloatingSubplot(fig, arr_pos[0],arr_pos[1],arr_pos[2], grid_helper=grid_helper)
    fig.add_subplot(ax1)
    ax1.axis["bottom"].set_axis_direction("bottom")  # "Angle axis"
    #ax1.axis["top"].toggle(ticklabels=True, label=True)
    ax1.axis["bottom"].major_ticklabels.set_axis_direction("top")
    ax1.axis["bottom"].major_ticklabels.set_pad(10)
    ax1.axis["bottom"].label.set_axis_direction("top")
    ax1.axis["bottom"].label.set_text("Correlation")
    ax1.axis["left"].set_axis_direction("top") # "X axis"
    #ax1.Labels.get_texts_widths_heights_descents()
    if if_rel:
      ax1.axis["left"].label.set_text("Relative Std")
    else:
      ax1.axis["left"].label.set_text("Standard deviation")

    ax1.axis["right"].set_axis_direction("bottom")   # "Y axis"
    ax1.axis["right"].toggle(ticklabels=True)
    ax1.axis["right"].major_ticklabels.set_axis_direction("bottom")
    ax1.axis["right"].major_ticklabels.set_pad(10)
    #ax1.TickLabels.get_texts_widths_heights_descents()
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
    
    return ax1, aux_ax

def add_sample(ax, arr_Cor=[], arr_Std=[], color=['r'], marker=['o'], size=100, alpha=1.0, if_rel=False, num_refstd=1.0, arr_text=[], lw=0.1, ec=[],fc=[] ):
  if len(fc)==0:
    fc=color
  if len(ec)==0:
    ec=color
  lenarr1 = len(arr_Cor)
  lenarr2 = len(arr_Std)
  if lenarr1 == lenarr2:
    arr_radius = [0 for n in range(lenarr1)]
    arr_theta  = [0 for n in range(lenarr1)]
    for n in range(lenarr1):
      arr_theta [n] = (1.0-arr_Cor[n]) * 0.5 * math.pi
      if if_rel:
        if num_refstd == 0:
          arr_Std [n] = 0.0
        else:
          arr_Std [n] = arr_Std[n]/num_refstd 
    #ax.scatter(arr_theta, arr_Std, c=color,marker=marker, s=marker_size, alpha=0.0)
      ax.scatter(arr_theta[n], arr_Std[n], s=size, c=color[n],marker=marker[n], alpha=alpha, linewidths=lw, edgecolors=ec[n], zorder=3,facecolor=fc[n] )
      if len(arr_text) > 0: 
        ax.text   (arr_theta[n], arr_Std[n], arr_text[n], fontsize=10)
    return True  
  else:
    print("!! length of arr_Cor is not equal to arr_Std")

def add_cntour(ax, num_refstd, num_intvl=0.25, num_cnts=10, alpha=0.5, color='#b4b4b4', marker='o', if_rel=False):
  if if_rel: num_refstd = 1.0
  arr_psi = np.arange(0, math.pi, 0.05)
  arr_r   = [[] for n in range(num_cnts)]
  arr_theta   = [[] for n in range(num_cnts)]
  for n in range(num_cnts):
    for psi in arr_psi:
      r     = num_intvl * (n+1) * num_refstd
      R     = (num_refstd**2 - 2 * num_refstd * r * math.cos( psi) + r **2) ** 0.5
      if R <=0:
        theta = 0
      else:
        theta =  math.asin( r * math.sin(psi) / R )
      arr_r[n]    .append( R     )
      arr_theta[n].append( theta )
    ax.plot(arr_theta[n], arr_r[n], c=color, ls='--', alpha=alpha, zorder=2)
  return ax


##########################################################
a = """
example
fig = plt.figure(1, figsize=(8, 4))
fig.subplots_adjust(wspace=0.3, left=0.05, right=0.95)
refstd = 3.0
r = 0.5
ax2,aux_ax2 = Taylor_Diagram( fig, 111, refstd )
add_sample(aux_ax2, arr_Cor=[0.5, 0.7, 0.99] ,arr_Std=[2.0, 2.3, 1.5]   )
add_cntour(aux_ax2, refstd)
plt.show()
"""
