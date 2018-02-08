import numpy as np
import matplotlib.pyplot as plt
from params import *

# revise matplotlib rc parameters to customize features of figure, like line width, axis width, label fontsize....
myparams = {
          'axislw': ('axes.linewidth', 4),                       # axes.linewidth
          'titlesz': ('axes.titlesize', 28),                     # axes.titlesize, title size of axes
          'linewd': ('lines.linewidth', 3),                      # lines.linewidth  //Line size of Line2D instance
          'ticksz': ('xtick.major.size', 'ytick.major.size', 6),  # xtick/ytick.major.size    //Tick size of x,y axis
          'tickwd': ('xtick.major.width', 'ytick.major.width', 2),  # xtick/ytick.major.width     // Tick width of x,y axis
          'labelsz': ('axes.labelsize', 36),                     # axes.labelsize             // label size of x,y axis
          'ticklbsz': ('xtick.labelsize', 'ytick.labelsize', 24), # xtick/ytick.labelsize     // Tick label size of x,y axis
          'legendsz': ('legend.fontsize', 20),                   # legend.fontsize          // Font size of legend
          'figsize': ('figure.figsize', (16, 10)),                # figure.figsize
          'figdpi' : ('figure.dpi', 100),                        # figure,dpi
          'tickpad': ('xtick.major.pad', 'ytick.major.pad', 20)}   # xtick/ytick.major.pad     // Distance of tick label to ticks
parm = Setparm(myparams)
parm.setrc()


def set_format(ax):
    ax.xaxis.labelpad = 25
    ax.yaxis.labelpad = 25
    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()
    yt = ax.yaxis.get_ticklabels()
    yt[0].set_visible(False)


def set_minortick_interval(ax, axis, divided_num):
    majortick_locs = eval('ax.%saxis.get_ticklocs()' % axis)
    majortick_interval = majortick_locs[1]-majortick_locs[0]
    minortick_interval = majortick_interval/divided_num
    


def set_majortick_interval(ax, axis, interval):
    ticklocs = eval('ax.%saxis.get_ticklocs()' % axis)
    maxvalue = ticklocs[-1]
    current = ticklocs[0]
    newtickloc = []
    newticklabel = []
    newtickloc.append(current)
    newticklabel.append(str(current))
    while current < maxvalue:
        current += interval
        newtickloc.append(current)
        newticklabel.append(str(current))
    eval('plt.%sticks(newtickloc,newticklabel)' % axis)


def set_ticklabel(ax, axis, label):
    ticklocs = eval('ax.%saxis.get_ticklocs()' % axis)
    ticklocs = list(ticklocs)
    newticklabels = []
    loc_set = label.keys()
    for key in label.keys():
        if key not in ticklocs:
            ticklocs.append(key)
    ticklocs.sort()
    newticklabels = [label[x] if x in label.keys() else str(x) for x in ticklocs]
    eval('plt.%sticks(ticklocs, newticklabels)' % axis)


def addRec(f, ax, dimension, color='k'):
    global AxesPercentage, FigPercentage
    left, bottom, width, height = dimension
    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()
    pos = ax.get_position()
    AxesWidth = xmax - xmin
    AxesHeight = ymax - ymin
    # AxesPercentage = [(left-xmin)/AxesWidth,(bottom-ymin)/AxesHeight,width/AxesWidth,height/AxesHeight]
    FigLeft = pos.x0+(left-xmin)/AxesWidth*pos.width
    FigBottom = pos.y0+(bottom-ymin)/AxesHeight*pos.height
    FigWidth = width/AxesWidth*pos.width
    FigHeight = height/AxesHeight*pos.height
    FigPercentage = [FigLeft, FigBottom, FigWidth, FigHeight]
    ax2 = f.add_axes(FigPercentage)
    rec = ax2.patch
    rec.set_facecolor(color)
    ax2.xaxis.set_visible(False)
    ax2.yaxis.set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    return ax2

### To set one side tick line, one can use ax.xaxis.tick_bottom(),ax.yaxis,tick_left() Only show tick line on bottom for x axis and on left for y axis
### Or one can use for loop to set properties of ticks. for tick in ax.yaxis.get_major_ticks(): tick.tick2On = False
### To change tick location and label. function xtick or ytick can be used. xtick([,,,,,],[,,,,]) the firt list includes position, and the second list includes string of tick label.
