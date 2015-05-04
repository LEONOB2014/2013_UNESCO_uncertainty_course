# -*- coding: utf-8 -*-
"""
Created on Sat May 04 12:20:08 2013

@author: VHOEYS
"""

import numpy as np
import matplotlib.pyplot as plt

#IMPORT NORMAL DISTRIBUTION:
from scipy.stats import norm
#other distributions: http://docs.scipy.org/doc/scipy/reference/stats.html


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
