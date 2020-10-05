'''
Each model is referred to using a modelname and must contain must contain three methods
  intializemodelname
  modelname
  fitmodelname
T2SE : Simple diffusion exponential decay model 
last modification: 6-3-14
'''

import lmfit
import numpy as np


def init (bvalue=None, data=None):    
    """initialize parameters for DifModel """
    nDifModelparams =2      #max number of parameters, some may be fixed
    params = lmfit.Parameters()   #define parameter dictionary
    paramlist = []    # list of parameters used for this model
    ADCguess=0.001
    params.add('ADC', value= ADCguess, min=0., vary = True)
    paramlist.append('ADC')
    params.add('Si', value= np.amax(data), vary = True)
    paramlist.append('Si')
#     params.add('B',  value= 0,  min=0., vary = False)  #baseline for Rician noise correction, not implemented
#     paramlist.append('B')
    return [params,paramlist]

# define objective function: returns the array to be minimized
def objFunction(params, bvalue, data):
    """ monoExponential Diffusion model; bvalue  array"""
    Si = params['Si'].value
    ADC = params['ADC'].value

    model = Si*np.exp(-bvalue*ADC)
    return (model - data)

def fit(params, bvalue, data):
    """fits signal vs bvalue data to DifModel model"""
    result = lmfit.minimize(DifModel, params, args=(bvalue, data))
    final = data + result.residual
    return final

