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


def initializeT1IRabs (nroi=None,TI=None, data=None, roi = None, useROIs = False, delta=None):    
    """initialize parameters for T1IR absolute value model"""
    nT1IRparams =3      #max number of parameters, some may be fixed
    if nroi == None:    #if no parameters are passed return the number of fitting parameters for this model
      return nT1IRparams
    T1params = lmfit.Parameters()   #define parameter dictionary
    paramlist = []    # list of parameters used for this model
    if useROIs: #if true use ROI values else use best guess
      T1guess = roi.T1  
    else:
      T1guess=TI[np.argmin(data)]/np.log(2) #minimum signal should occur at ln(2)T1
    T1params.add('T1', value= T1guess, min=0, vary = True)
    paramlist.append('T1')
    T1params.add('Si', value= np.amax(data), vary = True)
    paramlist.append('Si')
    if delta==None:
        T1params.add('B',  value= 2,  min=1.5, max=3.0, vary = True)
    else:
        T1params.add('B',  value= 1+delta,  min=1.5, max=2.5, vary = False)
    paramlist.append('B')
#     T1params.add('A',  value= 0,  min=0,  vary = True)
#     paramlist.append('A')
    return [T1params,paramlist]

# define objective function: returns the array to be minimized
def T1IRabs(params, TI, data):
    """ T1-IR model abs(exponential); TI inversion time array, T1 recovery time"""
    #A = params['A'].value
    B = params['B'].value
    Si = params['Si'].value
    T1 = params['T1'].value
    
    model = np.abs(Si*(1-B * np.exp(-TI/T1)))
    obf=model-data
    return (obf)

def fitT1IRabs(params, TI, data):
    """fits signal vs TI data to T1IRabs model"""
    result = lmfit.minimize(T1IRabs, params, args=(TI, data))
    final = data + result.residual
    return final
