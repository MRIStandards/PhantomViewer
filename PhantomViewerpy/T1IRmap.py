"""
Created on Fri Oct 11 16:30:54 2013
Each model is referred to using a modelname and must contain must contain three methods
  intializemodelname
  modelname
  fitmodelname
T1IRabd : T1 inversion recovery absolute value model 
last modification: 6-3-14
"""

import lmfit
import numpy as np


def init (TI=None, data=None):    
    """initialize parameters for T1IR absolute value model"""
    nT1IRparams =3      #max number of parameters, some may be fixed
    params = lmfit.Parameters()   #define parameter dictionary
    paramlist = []    # list of parameters used for this model
    T1guess=TI[np.argmin(data)]/np.log(2) #minimum signal should occur at ln(2)T1
    params.add('T1', value= T1guess, min=0, vary = True)
    paramlist.append('T1')
    params.add('Si', value= np.amax(data), vary = True)
    paramlist.append('Si')
    params.add('delta',  value= 1,  min=0, max=2.0, vary = True)
    return [params,paramlist]

# define objective function: returns the array to be minimized
def objFunction(params, TI, data):
    """ T1-IR model abs(exponential); TI inversion time array, T1 recovery time"""
    delta = params['delta'].value
    Si = params['Si'].value
    T1 = params['T1'].value
    model = np.abs(Si*(1-(1+delta)*np.exp(-TI/T1)))
    return (model - data)

def fit(params, TI, data):
    """fits signal vs TI data to T1IRabs model"""
    result = lmfit.minimize(objFunction, params, args=(TI, data))
    final = data + result.residual
    return final
