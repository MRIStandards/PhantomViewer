'''
Each model is referred to using a modelname and must contain must contain three methods
  intializemodelname
  modelname
  fitmodelname
Gaussian 
last modification: 6-3-14
'''

import lmfit
import numpy as np


def initialize (xi, data):    
    """initialize parameters for model"""
    nparams =4      #max number of parameters, some may be fixed

    params = lmfit.Parameters()   #define parameter dictionary
    paramlist = []    # list of parameters used for this model
    xmin=np.amin(xi)
    xmax=np.amax(xi)
    params.add('x0', value= np.average(xi),  vary = True) #center  min=xmin, max=xmax,
    paramlist.append('x0')
    params.add('sigma', value= (xmax-xmin)/6, min=0, max=xmax-xmin ,vary = True)
    paramlist.append('sigma')
    params.add('A',  value= np.amax(data), vary = True)
    paramlist.append('A')
    params.add('S0',  value= 0, vary = False)
    paramlist.append('S0')
#     params.update_constraints()
#     params.update()
    return [params,paramlist]

# define objective function: returns the array to be minimized
def model(params, xi, data):
    """ gaussian model; """
    x0 = params['x0'].value
    sigma = params['sigma'].value
    A = params['A'].value
    S0 = params['S0'].value

    modelValues = S0+A*np.exp(-((xi-x0)/sigma)**2/2)
    return (modelValues - data)

def fit(params, xi, data):
    """fits signal vs x data to gaussian model"""
    result = lmfit.minimize(model, params, args=(xi, data))
    final = data + result.residual
    return result

