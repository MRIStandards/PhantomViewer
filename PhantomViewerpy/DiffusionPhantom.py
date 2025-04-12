
# -*- coding: utf-8 -*-
"""
Created on Jan 19, 2015
module containing parameters describing NIST/ISMRM MRI system phantom
@author: stephen russek
"""
import numpy as np
import VPhantom
try:
  import Image
except:
  pass
 
nROIs = 13    # number of ROIs
angle = 60*np.pi/180
innerArrayRadius = 35.       #contrast arrays have 10 elements equally spaced on a circle or radius array radius in mm
outerArrayRadius = 60. #4 inner contrast spheres are on a square of size ArrayInnerSpacing in mm
DifArrayY = 0        #y position of arrays in mm
PhantomRadius = 98. #radius of phantom in mm
tubeRadius = 15. #inner radius of fiducial spheres in mm
WaterR1 = 0.3      #water R1 in s^-1
WaterR2 = 0.5      #water R2 in s^-1
Tref = 215.05   #water diffusion reference temperature (K) from Holz 2000
D0=1.635E-8     #water reference diffusion coeff m2/s
gammawater= 2.062
Tphantom = 0    #phantom temperature in C
DataTypes={"ADC":0,"T1":1,"T2":2}
DifArrayConcentration=[0,0,0,10,10,20,20,30,30,40,40,50,50]   #%PVP by weight
DifArrayT1=[2047.7,1448.1,1024.,724.1,512.,362.1,256.,181.,128.,90.5,64.,45.3,32.]    #T1 in ms
DifArrayADC=[1.099,1.099,1.099,0.860,0.860,0.58,0.58,0.42,0.42,0.275,0.275,0.2, 0.2]  #ADC in 10^-3 mm2/s
DiffusionPhantomImage = r"..\icons\DiffusionPhantom.jpg"

class DiffusionPhantom(VPhantom.VPhantom):
  """A virtual phantom that describes NIST-RSNA Isotropic Diffusion Phantom"""
  def __init__(self):
    VPhantom.VPhantom.__init__(self)
    self.phantomName = "NIST-RSNA Isotropic Diffusion Phantom"
    self.ROIsetdict = {'DifArray':0}
    self.DifROIs =self.SetDefaultContrastROIs("Dif")    #sets roi positions
    self.DifROIs.Comment = "PVP array"
    for roi in self.DifROIs.ROIs:    #sets T1, concentration
      roi.T1=DifArrayT1[roi.Index-1]
      roi.ADC=DifArrayADC[roi.Index-1]/1000   #note ADCs are in units of mm2/s
      roi.Concentration=DifArrayConcentration[roi.Index-1]
    self.ROIsets.append(self.DifROIs)
    self.phantomImage=DiffusionPhantomImage
    
    
  def SetDefaultContrastROIs(self,ptype):
  #set initial system phantom contrast ROIs
      ROIRadius = tubeRadius / 2 #set initial ROI radius in mm
      r=VPhantom.ROISet(ptype)
      r.ROIName =  ptype + "Array"
      r.Field = 1.5
      r.Temperature = 0.
      r.nROIs=nROIs

      for I in range (1, nROIs +1): #make 13 ROIs
          r.ROIs.append(VPhantom.ROI())
          r.ROIs[-1].Name = ptype + "-" +str(I)
          r.ROIs[-1].Index = I    #note the index starts at one while the list index starts at 0
          r.ROIs[-1].ROIRadius = ROIRadius
          r.ROIs[-1].Xcenter = 0.
          r.ROIs[-1].Zcenter = 0.
          r.ROIs[-1].Ycenter = DifArrayY
  
  #manually set ROI positions    
      r.ROIs[0].Xcenter = 0.
      r.ROIs[0].Zcenter = 0.
      
      r.ROIs[1].Xcenter = innerArrayRadius*np.cos(angle/2)
      r.ROIs[1].Zcenter = innerArrayRadius*np.sin(angle/2)
      
      r.ROIs[2].Xcenter = outerArrayRadius*np.cos(0)
      r.ROIs[2].Zcenter = outerArrayRadius*np.sin(0)
      
      r.ROIs[3].Xcenter = innerArrayRadius*np.cos(-angle/2)
      r.ROIs[3].Zcenter = innerArrayRadius*np.sin(-angle/2)
      
      r.ROIs[4].Xcenter = outerArrayRadius*np.cos(-angle)
      r.ROIs[4].Zcenter = outerArrayRadius*np.sin(-angle)
      
      r.ROIs[5].Xcenter = innerArrayRadius*np.cos(-1.5*angle)
      r.ROIs[5].Zcenter = innerArrayRadius*np.sin(-1.5*angle)
      
      r.ROIs[6].Xcenter = outerArrayRadius*np.cos(-2*angle)
      r.ROIs[6].Zcenter = outerArrayRadius*np.sin(-2*angle)
      
      r.ROIs[7].Xcenter = innerArrayRadius*np.cos(-2.5*angle)
      r.ROIs[7].Zcenter = innerArrayRadius*np.sin(-2.5*angle)
      
      r.ROIs[8].Xcenter = outerArrayRadius*np.cos(3*angle)
      r.ROIs[8].Zcenter = outerArrayRadius*np.sin(3*angle)
      
      r.ROIs[9].Xcenter = innerArrayRadius*np.cos(2.5*angle)
      r.ROIs[9].Zcenter = innerArrayRadius*np.sin(2.5*angle) 
      
      r.ROIs[10].Xcenter = outerArrayRadius*np.cos(2*angle)
      r.ROIs[10].Zcenter = outerArrayRadius*np.sin(2*angle)  
      
      r.ROIs[11].Xcenter = innerArrayRadius*np.cos(1.5*angle)
      r.ROIs[11].Zcenter = innerArrayRadius*np.sin(1.5*angle) 
      
      r.ROIs[12].Xcenter = outerArrayRadius*np.cos(angle)
      r.ROIs[12].Zcenter = outerArrayRadius*np.sin(angle) 
      
      return r
    
    
