'''
Each model is referred to using a modelname and must contain must contain three methods
  intializemodelname
  modelname
  fitmodelname
fits to a liquid crystal phase transition used to determine temperature. 
Assume the order parameter goes to zero as the temperature increase as (Tc-T)^gamma and the signal 
is proportional to the order parameter plus  constant offset 
last modification: 6-3-14
'''

import lmfit
import numpy as np


def initialize (Tc, data):    
    """initialize parameters for model"""
    nparams =4      #max number of parameters, some may be fixed

    params = lmfit.Parameters()   #define parameter dictionary
    paramlist = []    # list of parameters used for this model
    params.add('T', value= 20,  vary = True) #sample temperature
    paramlist.append('T')
    params.add('gamma', value= 0.2, min=0, max=1 ,vary = True)  #critical exponent
    paramlist.append('gamma')
    params.add('A',  value= 100, vary = True)   #~Signal max
    paramlist.append('A')
    params.add('B',  value= 0, vary = True)   #baseline
    paramlist.append('B')
    return [params,paramlist]

# define objective function: returns the array to be minimized
def model(params, Tc, data):
    """ Phase transition power law """
    T = params['T'].value
    gamma = params['gamma'].value
    A = params['A'].value
    B = params['B'].value
    modelValues = np.nan_to_num(A*(T-Tc)**gamma +B, nan=B)
    return (modelValues - data)

def fit(params, xi, data):
    """fits signal vs x data to power law model"""
    result = lmfit.minimize(model, params, args=(xi, data))
    final = data + result.residual
    return result

