'''
This module contains functions that draw user-defined matplotlib objects on a 
given axes.
'''

import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, Rectangle, Arc, Polygon
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker


def moon_phase(x,y,alpha_d,dx,c_br,c_dk,lw,zorder,ax=None):
    '''
    To draw moon/planet phases.

    Inputs
    ------
    x,y: center coordinates in axis units
    alpha: [deg], [-pi,pi], the phase angle (Sun-object-observer angle). 
        0=full, pi/2=east half (bright), +-pi=crescent, -pi/2=west half
    dx: diameter of the demo circle in x-axis units
    c_br: color of the bright part
    c_dk: color of the dark part
    lw: linewidth of the circle
    zorder: lowest zorder
    ax: the axes object to be drawn on

    Outputs
    -------
    e1,e2,rect,cir: matplotlib objects
    '''

    alpha = alpha_d*np.pi/180

    # get x-to-y conversion factor 
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    bbox_w = ax.bbox.width 
    bbox_h = ax.bbox.height 
    ratio = (ymax-ymin)/bbox_h*bbox_w/(xmax-xmin)

    # The bottom ellipse
    e1 = Ellipse((x,y),dx,ratio*dx,0,facecolor=c_dk,edgecolor='none',
        zorder=zorder)

    # The top ellipse 
    c = None 
    if np.abs(alpha)<np.pi/2:
        c = c_br 
    else:
        c = c_dk
    e2 = Ellipse((x,y),dx*np.abs(np.cos(alpha)),ratio*dx,facecolor=c,
        edgecolor='none',zorder=zorder+2)

    # the shading rectangle
    xy = None
    if alpha>0:
        xy = (x-dx/2, y-ratio*dx/2)
    else:
        xy = (x, y-ratio*dx/2)
    rect = Rectangle(xy,dx/2,ratio*dx,facecolor=c_br,edgecolor='none',
        zorder=zorder+1)

    # The circle 
    cir = Ellipse((x,y),dx,ratio*dx,facecolor='none',edgecolor=c_dk,
        linewidth=lw,zorder=zorder+3)

    # draw
    if not ax==None:
        ax.add_artist(e1)
        ax.add_artist(e2)
        ax.add_artist(rect)
        ax.add_artist(cir)

    return e1,e2,rect,cir


def multicolor_label(ax,strings,colors,axis='x',anchorpad=0,**kw):
    '''
    Source: 'https://stackoverflow.com/questions/33159134/
        matplotlib-y-axis-label-with-multiple-colors'
    This function creates axes labels with multiple colors
    ax: specifies the axes object where the labels should be drawn
    strings: list of all of the text items
    colors: list of colors for the strings
    axis: 'x', 'y', or 'both' and specifies which label(s) should be drawn
    '''

    # x-axis label
    if axis=='x' or axis=='both':
        boxes = [TextArea(text,textprops=dict(color=color,ha='left',va='bottom',
            **kw)) for text,color in zip(strings,colors) ]
        xbox = HPacker(children=boxes,align="center",pad=0, sep=5)
        anchored_xbox = AnchoredOffsetbox(loc=3, child=xbox, pad=anchorpad,
            frameon=False,bbox_to_anchor=(0.2, -0.09),
            bbox_transform=ax.transAxes, borderpad=0.)
        ax.add_artist(anchored_xbox)

    # y-axis label
    if axis=='y' or axis=='both':
        boxes = [TextArea(text, textprops=dict(color=color, ha='left',
            va='bottom',rotation=90,**kw)) 
            for text,color in zip(strings[::-1],colors) ]
        ybox = VPacker(children=boxes,align="center", pad=0, sep=5)
        anchored_ybox = AnchoredOffsetbox(loc=3, child=ybox, pad=anchorpad, 
            frameon=False, bbox_to_anchor=(-0.10, 0.2), 
            bbox_transform=ax.transAxes, borderpad=0.)
        ax.add_artist(anchored_ybox)



def test():
    ''' moon_phase

    # parameters 
    alpha_d = -23

    plt.figure(figsize=(5,5))
    ax = plt.axes([0,0,1,1])
    plt.xlim(-.6,.6)
    plt.ylim(-.6,.6)
    moon_phase(0,0,alpha_d,1,'w','k',lw=1,zorder=1,ax=ax)
    plt.show()
    #'''





