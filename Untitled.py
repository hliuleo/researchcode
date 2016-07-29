# -*- coding: utf-8 -*-
"""
Created on Thu May 21 19:32:29 2015

@author: hliu
"""

def returnLoc(f,ax,dimension,color='k'):
    global AxesPercentage, FigPercentage
    left,bottom,width,height = dimension
    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()
    pos = ax.get_position()
    AxesWidth = xmax - xmin
    AxesHeight = ymax - ymin
    #AxesPercentage = [(left-xmin)/AxesWidth,(bottom-ymin)/AxesHeight,width/AxesWidth,height/AxesHeight]
    FigLeft = pos.x0+(left-xmin)/AxesWidth*pos.width
    FigBottom = pos.y0+(bottom-ymin)/AxesHeight*pos.height
    FigWidth = width/AxesWidth*pos.width
    FigHeight = height/AxesHeight*pos.height
    FigPercentage = [FigLeft, FigBottom, FigWidth, FigHeight]

    return FigPercentage
    
x = np.arange(3)
plt.plot(x,x)
ax = plt.gca()
f = plt.gcf()
ax2 = f.add_subplot(122)
dimension = [1,-1,0.5,0.5]
loc = returnLoc(f,ax,dimension)
#ax2.set_position(loc)
ax2.set_position(loc)
plt.savefig('/Users/hliu/Desktop/test.png')#,bbox_inches='tight')