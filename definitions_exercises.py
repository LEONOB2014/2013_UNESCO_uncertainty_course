# -*- coding: utf-8 -*-
"""
Created on Sat May 04 14:20:51 2013

@author: VHOEYS
"""

#Load packages for calculation and making plots
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

#IMPORT NORMAL DISTRIBUTION:
from scipy.stats import norm

from matplotlib import cm, colors
import  matplotlib.axes as maxes 

from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid import make_axes_locatable

#------------------------------------------------------------------------------
#DEFINITIONS TO RUN THE MODEL--------------------------------------------------
def deriv_works(u,t,Pars,Const):
    '''
    Differential equations of the respirometric model
    '''
    #Define the parameters
    mumax = np.float64(Pars[0])
    Y = np.float64(Pars[1])
    Ks = np.float64(Pars[2])
    tau = np.float64(Pars[3])
##    print 'mumax is %f ;Y is %f; Ks is %f and tau is %f' %(mumax,Y,Ks,tau)
    b = np.float64(Const[0])
    kla = np.float64(Const[1])
    SOeq = np.float64(Const[2])
##    print ' de kla is %f en b is %f en SOeq is %f' %(kla,b,SOeq)
    Monod=mumax*(u[1])/(u[1]+Ks)    #Monod Kinetic
    Expo=1.0-np.exp(-t/tau)         #Helpfunction

    dXdt = (Expo*Monod-b)*u[0]                    #Biomassa
    dSsdt = -(1.0/Y)*Expo*Monod*u[0]              #Substraat
    dOdt = kla*(SOeq-u[2])-((1-Y)/Y)*Expo*Monod*u[0]   #Oxygen

    return np.array([dXdt,dSsdt,dOdt]) #

def RespiroModel(Pars,Init_unc,time):
    '''
    Run the respirometric model
    '''
    #Define the constants
    b = 0.62
    kla = 369.7334962
    SOeq = 8.4
    Constt = np.array([b,kla,SOeq])

    #Define the initial conditions (Constants)Ss0
    Ss0 = 58.4899
    #Define the initial conditions (Uncertain) -> X0
    X0=Init_unc[0]
    yinit = np.array([X0,Ss0,SOeq])  

    #Define the necessary parameters
    mumax = np.float64(Pars[0])
    Y = np.float64(Pars[1])
    Ks = np.float64(Pars[2])
    tau = np.float64(Pars[3])

    #Solve with LSODA scheme
    y,infodic=odeint(deriv_works,yinit,time,full_output=True, printmessg=False, args=(Pars,Constt))

    #Get outputs
    X = y[:,0]
    Ss = y[:,1]
    O = y[:,2]

    OUR_ex=((1-np.exp(-time/tau))*mumax*(1-Y)/Y*Ss/(Ss+Ks)*X)/(24*60)
    #when using this deifnition, we get the time, biomass X, substrate S, oxygen O and Oxygen uptake rate OUR_ex back, 
    #together with some information about the model 
    return [time,X,Ss,O,OUR_ex,infodic]
#END OF DEFINITIONS TO RUN THE MODEL-------------------------------------------
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
#DEFINITIONS TO RUN THE SAMPLING-----------------------------------------------
def segments_part(nsegments):
    '''lhs help function
    Function used in our lhs-sampling procedures 
    
    Parameters
    -----------
    nsegments : int
        The number of samples to take
    
    Returns
    --------
    samples : narray
        An array with the different samples taken (in each segment between [0,1])
    '''
    #SAMPLING UNIFORM IN SEGMENTS----------------------------------------------
    #prepare the samples
    segments = np.linspace(0,1.,nsegments+1) 
    #we now sample as we would sample in a uniform distribution in each segment:
    samples=np.zeros(nsegments) #here we will store our sampled values
    for i in range(nsegments):
        samples[i]=np.random.uniform(segments[i],segments[i+1],1)
    #--------------------------------------------------------------------------    
    
    #we shuffle the values for higher order purposes---------------------------
    np.random.shuffle(samples)
    #--------------------------------------------------------------------------
    return samples

#function for Latin-hypercube sampling  of uniform distribution
def lhs_uniform(nsamples,pmin,pmax):
    '''lhs of uniform distribution
    
    '''
    #Run the help function to get the samples from the segments
    unif_samples = segments_part(nsamples)
    
    #transform into uniform samples of the requested interval [pmin, pmax]
    lhs_samples = pmin + unif_samples* (pmax-pmin)

    return lhs_samples

#function for Latin-hypercube sampling  of normal distribution
def lhs_normal(nsamples, ploc, pscale):   
    '''lhs of normal distribution
    
    '''
    #Run the help function to get the samples from the segments
    unif_samples = segments_part(nsamples)
    
    #transform into random samples with requested mean (ploc) and std (pscale)   
    norm_samples =  norm.ppf(unif_samples, ploc, pscale)
    
    return norm_samples
#END OF DEFINITIONS TO RUN THE SAMPLING----------------------------------------
#------------------------------------------------------------------------------



#RIVER MODEL
#------------------------------------------------------------------------------
#DEFINITIONS TO RUN THE MODEL--------------------------------------------------
def deriv_river(u,t,Pars, dummy):
    '''
    Differential equations of the river contamination model
    '''
    #Define the parameters
    k1 = np.float64(Pars[0])
    k2 = np.float64(Pars[1])
    BZVin = np.float64(Pars[2])
    DOsat = np.float64(Pars[3])

    dBZVdt = BZVin - k1*u[0]                        #BZV
    dDOdt = k2 *(DOsat - u[1])-k1*u[0]              #DO

    return np.array([dBZVdt,dDOdt])

def RiverModel(Pars,time):
    '''
    Run the river contamination model
    '''
    #Define the initial conditions (Constants)
    BZV0 = 7.33
    DO0 = 8.5
    yinit = np.array([BZV0, DO0])  

    #Solve with LSODA scheme
    y,infodic=odeint(deriv_river,yinit,time,full_output=True, printmessg=False, args=(Pars,[]))

    #Get outputs
    BZV = y[:,0]
    DO = y[:,1]

    return [time,BZV,DO,infodic]
#END OF DEFINITIONS TO RUN THE MODEL-------------------------------------------
#------------------------------------------------------------------------------


#GLUE VISUALISATION
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
        
    return fig,axScatter,axHistx,axHisty

def UncertOut(ObjFunc_behav,Output_behav, qval=0.05):
    '''
    Help function of the plotfunction UncPlot

    OF matrix is nx1 matrix
    Output_behav is nxkTimesteps

    '''

    Normed=ObjFunc_behav/ObjFunc_behav.sum()
    #print 'controle voor de som: %f' %Normal.sum()
    Likl=Normed

    qvalmin=qval
    qvalmax=1.0-qval
    lbm=[]
    ubm=[]

    for i in range(np.size(Output_behav,1)):  #! voor elke tijdstap
        Outp=Output_behav[:,i]

        tempi = np.vstack([Likl,Outp]) #omgezet naar rijvectoren...
        indices = np.lexsort(tempi)  #(last row is sortrow in default)
        cdf=tempi.take(indices, axis=-1)
        cdf[0,:]=np.cumsum(cdf[0,:])  #2xn matrix
        
        #Make first likelihood = 0.0 to get cdf
        cdf[0,0]=0.0
        
        lb=np.nonzero(cdf[0,:]>qvalmin)[0][0]
        ub=np.nonzero(cdf[0,:]<qvalmin)[0][-1]
        
        if (np.abs(qval-(cdf[0,:].take([lb]))))<(np.abs(qval-(cdf[0,:].take([ub])))):
            lbselect=lb
        else:
            lbselect=ub
            
        #upper boundary -indice
        lb=np.nonzero(cdf[0,:]>qvalmax)[0][0]
        ub=np.nonzero(cdf[0,:]<qvalmax)[0][-1]           
        if (np.abs(qval-(cdf[0,:].take([lb]))))<(np.abs(qval-(cdf[0,:].take([ub])))):
            ubselect=lb
        else:
            ubselect=ub

        lbm=np.append(lbm,cdf[1,:].take([lbselect]))
        ubm=np.append(ubm,cdf[1,:].take([ubselect]))

    return [lbm,ubm]