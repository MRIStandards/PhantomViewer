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


def initializeDifModel (nroi=None,bvalue=None, data=None, roi = None, useROIs = False, B=0):    
    """initialize parameters for DifModel """
    nDifModelparams =2      #max number of parameters, some may be fixed
    if nroi == None:    #if no parameters are passed return the number of fitting parameters for this model
      return nDifModelparams
    ADCparams = lmfit.Parameters()   #define parameter dictionary
    paramlist = []    # list of parameters used for this model
    if useROIs: #if true use ROI values else use best guess
      ADCguess = roi.ADC  
    else:
      ADCguess=0.001
    ADCparams.add('ADC', value= ADCguess, min=0., vary = True)
    paramlist.append('ADC')
    ADCparams.add('Si', value= np.amax(data), vary = True)
    paramlist.append('Si')
#     ADCparams.add('B',  value= 0,  min=0., vary = False)  #baseline for Rician noise correction, not implemented
#     paramlist.append('B')
    return [ADCparams,paramlist]

# define objective function: returns the array to be minimized
def DifModel(params, bvalue, data):
    """ Diffusion model; bvalue  array"""
    Si = params['Si'].value
    ADC = params['ADC'].value

    model = Si*np.exp(-bvalue*ADC)
    return (model - data)

def fitDifModel(params, bvalue, data):
    """fits signal vs bvalue data to DifModel model"""
    result = lmfit.minimize(DifModel, params, args=(bvalue, data))
    final = data + result.residual
    return final

