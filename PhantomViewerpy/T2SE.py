'''
Each model is referred to using a modelname and must contain must contain three methods
  intializemodelname
  modelname
  fitmodelname
T2SE : Simple T2 spin echo exponential decay model 
last modification: 6-3-14
'''

import lmfit
import numpy as np


def initializeT2SE (nroi=None,TE=None, data=None, roi = None, useROIs = False, B=0):    
    """initialize parameters for T2SE model"""
    nT2SEparams =3      #max number of parameters, some may be fixed
    if nroi == None:    #if no parameters are passed return the number of fitting parameters for this model
      return nT2SEparams
    T2params = lmfit.Parameters()   #define parameter dictionary
    paramlist = []    # list of parameters used for this model
    if useROIs: #if true use ROI values else use best guess
      T2guess = roi.T2  
    else:
      T2guess=200
    T2params.add('T2', value= T2guess, min=0., vary = True)
    paramlist.append('T2')
    T2params.add('Si', value= np.amax(data), vary = True)
    paramlist.append('Si')
    T2params.add('B',  value= 0,  min=0., vary = False)
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

