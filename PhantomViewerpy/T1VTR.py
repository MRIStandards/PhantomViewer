
"""
Created on Jan 7, 2015
Each model is referred to using a modelname and must contain must contain three methods
  intializemodelname
  modelname
  fitmodelname
T1VTR : T1 variableflip angle  model 
last modification: 
"""

import lmfit
import numpy as np


def initializeT1VTR (nroi=None,TR=None, data=None, roi=None, ds = None, imageset = None):    
    """initialize parameters for T1 variable TR model"""
    nT1VTRparams =5      #max number of parameters, some may be fixed
    if nroi == None:    #if no parameters are passed return the number of fitting parameters for this model
      return nT1VTRparams
    T1params = lmfit.Parameters()   #define parameter dictionary
    paramlist = []    # list of parameters used for this model
    if roi == None:
      T1guess=100
    else:
      T1guess = roi.T1
    FA=ds.FA[imageset[0]]     # take fixed parameters from first image in image set. Cannot handle non-constant fixed parameters 
    TE=ds.TE[imageset[0]]
    B = 0.0
    T1params.add('T1', value= T1guess, min=1, vary = True)
    paramlist.append('T1')
    T1params.add('S0', value= np.amax(data), min=0, vary = True)
    paramlist.append('S0')
    T1params.add('FA', value= FA, min=0, vary = False)
    paramlist.append('FA')
    T1params.add('TE',  value= TE,  vary = False)
    paramlist.append('TE')
    T1params.add('B',  value= B,  vary = False)
    paramlist.append('B')
    return [T1params,paramlist]

# define objective function: returns the array to be minimized
def T1VTR(params, TR, data):  
    """ T1-VTR model ; TR repetition time"""
    FA = params['FA'].value*np.pi/180.
    TE = params['TE'].value
    T1 = params['T1'].value
    S0 = params['S0'].value
    B  = params['B'].value
    model = S0 * np.sin(FA) * (1 + np.exp(-TR / T1)-2*np.exp(-(TR-TE/2) / T1)) / (1 - np.exp(-TR / T1) * np.cos(FA)) + B
    return (model - data)

def fitT1VTR(params, xdata, ydata):
    """fits signal vs TI data to T1IRabs model"""
    result = lmfit.minimize(T1VTR, params, args=(xdata, ydata))
    final = ydata + result.residual
    return final
