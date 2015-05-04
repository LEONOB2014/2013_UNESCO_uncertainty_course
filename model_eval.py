# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 14:00:00 2013

@author: VHOEYS
"""

import os
#import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors
import  matplotlib.axes as maxes 

from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid import make_axes_locatable

def Scatter_hist(data1, data2, data1b=False, data2b=False, xbinwidth = 0.5, 
                 ybinwidth=0.5, SSE=None, SSEb=None, vmax=1000., colormaps=cm.YlOrBr,
                 cleanstyle = False, roodlichter=0.5, *args,  **kwargs):
    '''
    Three-parts plot with the two parameter sitdirbutions plotted on the sides
    
    Parameters
    -----------
    data1: ndarray
        dataset 1 for x-axis
    data2: ndarray
        dataset 2 for y-axis
    data1b: ndarray
        dataset to plot along the data 1 set
    data2b: ndarray
        dataset to plot along the data 2 set
    binwidth: float
        defines the width of the bins relative to the data used
    cleanstyle: bool True|False
        if True, a more minimalistic version of the plot is given
    *args,  **kwargs: args
        arguments given toi the scatter plot
        example: s=15, marker='o', edgecolors= 'k',facecolor = 'white'

    Returns
    ---------
    fig: matplotlib.figure.Figure object
        the resulting figure
    axScatter: matplotlib.axes.AxesSubplot object
        the scatter plot with the datapoints, can be used to add labels or 
        change the current ticks settings
    axHistx: matplotlib.axes.AxesSubplot object
        the x-axis histogram
    axHisty: matplotlib.axes.AxesSubplot object
        the y-axis histogram
    
    Examples
    ----------
    >>> nMC = 1000
    >>> parval1 = np.random.gamma(5.5,size=nMC)   
    >>> parval2 = np.random.gamma(8.0,size=nMC)   
    >>> parnames = ['par1','par2']
    >>> fig,axScatter,axHistx,axHisty = Scatter_hist(parval1,parval2, 
                                                 cleanstyle = True, s=48, 
                                                 marker='o', edgecolors= 'k', 
                                                 facecolor = 'none',alpha=0.7)    
    >>> parval1b = np.random.uniform(low=0.0, high=30.0,size=nMC)   
    >>> parval2b = np.random.uniform(low=0.0, high=30.0,size=nMC)   
    >>> fig,axScatter,axHistx,axHisty = Scatter_hist(parval1,parval2,parval1b, 
                                                 parval2b, cleanstyle = True, 
                                                 s=48, marker='o', 
                                                 edgecolors= 'k', 
                                                 facecolor = 'none',
                                                 alpha=0.7)
    
    Notes
    ------
    Typical application is to check the dependency of two posterior 
    parameter distrbutions, eventually compared with their selected posteriors
    
    If a second dataset is added to compare, the style options of the scatter
    plot are fixed and the *args, **kwargs have no influence
    
    '''
    if not isinstance(data1, np.ndarray):
        raise Exception('dataset 1 need to be numpy ndarray')
    if not isinstance(data2, np.ndarray):
        raise Exception('dataset 2 need to be numpy ndarray')        
        
    if isinstance(data1b, np.ndarray):
        if not isinstance(data2b, np.ndarray):
            raise Exception('Always combine the data of both')
    if isinstance(data2b, np.ndarray):
        if not isinstance(data1b, np.ndarray):
            raise Exception('Always combine the data of both')

    fig = plt.figure(figsize=(10,10))
    axScatter = plt.subplot(111)
    divider = make_axes_locatable(axScatter)
    #axScatter.set_aspect('equal')
    axScatter.set_autoscale_on(True)

    # create a new axes with  above the axScatter
    axHistx = divider.new_vertical(1.5, pad=0.0001, sharex=axScatter)

    # create a new axes on the right side of the
    # axScatter
    axHisty = divider.new_horizontal(1.5, pad=0.0001, sharey=axScatter)

    fig.add_axes(axHistx)
    fig.add_axes(axHisty)

    # now determine nice limits by hand:
#    binwidth = binwidth
    xmin = np.min(data1)
    xmax = np.max(data1)
    ymin = np.min(data2)
    ymax = np.max(data2)
    #xymax = np.max( [np.max(np.fabs(data1)), np.max(np.fabs(data2))] )
    #lim = (int(xymax/binwidth) + 1) * binwidth

    binsx = np.arange(xmin, xmax + xbinwidth, xbinwidth)
    binsy = np.arange(ymin, ymax + ybinwidth, ybinwidth)
    #bins = np.arange(-lim, lim + binwidth, binwidth)
    
    # the scatter plot:  
    if isinstance(data1b, np.ndarray): #TWO DATA ENTRIES
        if SSE == None:
            print '*args, **kwargs do not have any influcence when using two\
            options'
            axScatter.scatter(data1, data2, facecolor = 'none', 
                              edgecolor='k',s=25)
            axScatter.scatter(data1b, data2b, facecolor='none', 
                              edgecolor='grey',s=25)  

            xminb = np.min(data1b)
            xmaxb = np.max(data1b)
            yminb = np.min(data2b)
            ymaxb = np.max(data2b)   
            binsxb = np.arange(xminb, xmaxb + xbinwidth, xbinwidth)
            binsyb = np.arange(yminb, ymaxb + ybinwidth, ybinwidth)
    
            axHistx.hist(data1b, bins=binsxb, edgecolor='None', 
                         color='grey',normed=True)
            axHisty.hist(data2b, bins=binsyb, orientation='horizontal',  
                         edgecolor='None', color='grey', normed=True) 
    
            axHistx.hist(data1, bins=binsx, edgecolor='None', 
                         color='k', normed=True)
            axHisty.hist(data2, bins=binsy, orientation='horizontal',  
                         edgecolor='None', color='k', normed=True)
        else:                     
            print '*args, **kwargs do not have any influcence when using two\
            options'
            sc1 = axScatter.scatter(data1b, data2b, c=SSEb, vmax=vmax,alpha=roodlichter,
                              edgecolors= 'none', cmap = colormaps, *args,  **kwargs) 
                              
            axScatter.scatter(data1, data2, c=SSE, vmax=vmax, 
                              edgecolors= 'none', cmap = colormaps, *args,  **kwargs)
 

            xminb = np.min(data1b)
            xmaxb = np.max(data1b)
            yminb = np.min(data2b)
            ymaxb = np.max(data2b)   
            binsxb = np.arange(xminb, xmaxb + xbinwidth, xbinwidth)
            binsyb = np.arange(yminb, ymaxb + ybinwidth, ybinwidth)
    
            axHistx.hist(data1b, bins=binsxb, edgecolor='None', 
                         color=colormaps(1.),normed=True)
            axHisty.hist(data2b, bins=binsyb, orientation='horizontal', color=colormaps(1.),
                         edgecolor='None', normed=True) 
    
            axHistx.hist(data1, bins=binsx, edgecolor='None', 
                         color=colormaps(0.), normed=True)
            axHisty.hist(data2, bins=binsy, orientation='horizontal',  
                         edgecolor='None', color=colormaps(0.), normed=True)  
                                     
    else:  #ONLY ONE DATA1 and DATA2
        if SSE == None:
            axScatter.scatter(data1, data2, c= 'black', *args,  **kwargs)                        
            axHistx.hist(data1, bins=binsx, edgecolor='None', color='k')
            axHisty.hist(data2, bins=binsy, orientation='horizontal',  
                         edgecolor='None', color='k')            
        else:
            axScatter.scatter(data1, data2, c=SSE, vmax=vmax, 
                              edgecolors= 'none', cmap = colormaps, *args,  **kwargs)            
            axHistx.hist(data1, bins=binsx, edgecolor='None', color='k')
            axHisty.hist(data2, bins=binsy, orientation='horizontal',  
                         edgecolor='None', color='k')

    # the xaxis of axHistx and yaxis of axHisty are shared with axScatter,
    # thus there is no need to manually adjust the xlim and ylim of these
    # axis.

    majloc1 = MaxNLocator(nbins=4, prune='lower')   
    axScatter.yaxis.set_major_locator(majloc1)   
    majloc2 = MaxNLocator(nbins=4)  
    axScatter.xaxis.set_major_locator(majloc2) 

    axScatter.grid(linestyle = 'dashed', color = '0.75',linewidth = 1.)    
    axScatter.set_axisbelow(True)
    axHisty.set_axisbelow(True)
    axHistx.set_axisbelow(True)
    
    # The 'clean' environment
    if cleanstyle == True:
        plt.setp(axHistx.get_xticklabels() + axHisty.get_yticklabels(),
                 visible=False)   
        plt.setp(axHistx.get_yticklabels() + axHisty.get_xticklabels(),
                 visible=False)             
        axHisty.set_xticks([])
        axHistx.set_yticks([])
        axHistx.xaxis.set_ticks_position('bottom')
        axHisty.yaxis.set_ticks_position('left')
    
        axHisty.spines['right'].set_color('none')
        axHisty.spines['top'].set_color('none')
        axHisty.spines['bottom'].set_color('none')    
        axHisty.spines['left'].set_color('none') 
        
        axHistx.spines['top'].set_color('none')
        axHistx.spines['right'].set_color('none')
        axHistx.spines['left'].set_color('none') 
        axHistx.spines['bottom'].set_color('none') 
        axScatter.spines['top'].set_color('none') 
        axScatter.spines['right'].set_color('none') 
    else:
        plt.setp(axHistx.get_xticklabels() + axHisty.get_yticklabels(),
                 visible=False)      
        for tl in axHisty.get_yticklabels():
            tl.set_visible(False)
        for tlp in axHisty.get_xticklabels():
            tlp.set_rotation(-90)
            
        majloc3 = MaxNLocator(nbins=4, prune='lower')   
        axHistx.yaxis.set_major_locator(majloc3)  
        axHistx.yaxis.tick_right()
        majloc4 = MaxNLocator(nbins=4, prune='lower')   
        axHisty.xaxis.set_major_locator(majloc4) 
        axHisty.xaxis.tick_top()

        axHisty.yaxis.grid(linestyle = 'dashed', color = '0.75',linewidth = 1.)
        axHistx.xaxis.grid(linestyle = 'dashed', color = '0.75',linewidth = 1.)
        
    return fig,axScatter,axHistx,axHisty,sc1

##colormaps: see http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps
#http://www.loria.fr/~rougier/teaching/matplotlib/#scatter-plots

colormapss=cm.RdYlGn_r
treshold=20.

###LOAD DATA IN ZOOM1
datapath="D:\Dropbox\PieterV_Stijn"
SSE = np.loadtxt(os.path.join(datapath,'SSE_zoom1.txt'))
pars = np.loadtxt(os.path.join(datapath,'Parameters_zoom1.txt'))
bwx=15
bwy=0.03
###LOAD DATA IN ZOOM2
#SSE = np.loadtxt(os.path.join(datapath,'SSE_zoom2.txt'))
#pars = np.loadtxt(os.path.join(datapath,'Parameters_zoom2.txt'))
#bwx=150
#bwy=0.01

#select subset
behavpars = pars[np.where(SSE<treshold)]
behavSSE = SSE[np.where(SSE<treshold)]


#PLOT POSTERIORS
fig,axScatter,axHistx,axHisty,sc1 = Scatter_hist(behavpars[:,0], behavpars[:,1], data1b=pars[:,0],data2b=pars[:,1], 
             xbinwidth = bwx, ybinwidth=bwy,
             cleanstyle = True, s=25, marker='o', SSE=behavSSE, SSEb=SSE, 
             vmax=treshold+treshold*1.5, colormaps= colormapss, roodlichter=1.0)
axScatter.set_ylabel(r'r$_H$',fontsize=16)
axScatter.set_xlabel(r'v$_0$',fontsize=16)
cbar = fig.colorbar(sc1, ax=axScatter, cmap=colormapss, orientation='vertical',ticks=[treshold,treshold+treshold*1.5],shrink=1.)                        
cbar.ax.set_yticklabels(['<'+str(treshold),'> '+str(treshold+treshold*1.5)])




    