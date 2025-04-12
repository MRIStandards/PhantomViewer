'''
Each model is referred to using a modelname and must contain must contain three methods
  intializemodelname
  modelname
  fitmodelname
T2SE : Simple T2 spin echo exponential decay model, assumes fitting to magnitude data >0 
last modification: 6-3-14
'''

import lmfit
import numpy as np


def initializeT2SE (nroi=None,TE=None, data=None, roi = None, useROIs = False, B=0, varyB=False):    
    """initialize parameters for T2SE model"""
    nT2SEparams =3      #max number of parameters, some may be fixed
    if nroi == None:    #if no parameters are passed return the number of fitting parameters for this model
      return nT2SEparams
    T2params = lmfit.Parameters()   #define parameter dictionary
    paramlist = []    # list of parameters used for this model
    if useROIs: #if true use ROI values else use best guess
      T2guess = roi.T2  
    else:
      T2guess=200.0
    T2params.add('T2', value= T2guess, min=0.0, vary = True) #exponential time constant
    paramlist.append('T2')
    T2params.add('Si', value= np.amax(data), vary = True) #initial signal
    paramlist.append('Si')
    bmax=np.amax(data)/10
#     if varyB: #sets a maximum value for the baseline if it is allowed to vary
#       bmax=B
#     else:
#       bmax=np.inf
    T2params.add('B',  value= 0.0,min=0,max=bmax,vary = varyB)    #baseline must be >0
    paramlist.append('B')
    return [T2params,paramlist]

# define objective function: returns the array to be minimized
def T2SE(params, TE, data):
    """ T2-SE model; TE echo time array"""
    B = params['B'].value
    Si = params['Si'].value
    T2 = params['T2'].value

    model = Si*np.exp(-TE/T2)+ B
    return (model - data)

def fitT2SE(params, TE, data):
    """fits signal vs TE data to T2SE model"""
    result = lmfit.minimize(T2SE, params, args=(TE, data))
    final = data + result.residual
    return final

