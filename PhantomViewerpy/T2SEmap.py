"""
Created on Fri Oct 11 16:30:54 2013
Each model is referred to using a modelname and must contain must contain three methods
  intializemodelname
  modelname
  fitmodelname
T2SE : T2 spin echo model 
last modification: 6-3-14
"""

import lmfit
import numpy as np


def init (TE=None, data=None):    
    """initialize parameters for T2SE model"""
    nT1IRparams =2      #max number of parameters, some may be fixed
    params = lmfit.Parameters()   #define parameter dictionary
    paramlist = []    # list of parameters used for this model
    T2guess=np.mean(TE)
    params.add('T2', value= T2guess, min=0, vary = True)
    paramlist.append('T2')
    params.add('Si', value= np.amax(data), vary = True)
    paramlist.append('Si')
    return [params,paramlist]

# define objective function: returns the array to be minimized
def objFunction(params, TE, data):
    """ T2-SE model abs(exponential); TE echo time array, T2 spin transverse relaxation time"""
    Si = params['Si'].value
    T2 = params['T2'].value
    model = Si*np.exp(-TE/T2)
    return (model - data)

def fit(params, TE, data):
    """fits signal vs TE data to T2SE model"""
    result = lmfit.minimize(objFunction, params, args=(TE, data))
    final = data + result.residual
    return final
