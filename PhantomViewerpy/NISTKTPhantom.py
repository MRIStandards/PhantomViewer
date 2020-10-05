'''
Created on Jan 23, 2015

@author: stephen russek
'''
# -*- coding: utf-8 -*-

import numpy as np
import VPhantom 
nROIs = 8    # number of ROIs
angle = 60*np.pi/180
a = 21.2       #lattice constant
PhantomRadius = 60. #radius of phantom in mm
arrayRadius = 32.
tubeDiameter = 30. #tube diameter
innerROIDiameter = 6.
WaterR1 = 0.3      #water R1 in s^-1
WaterR2 = 0.5      #water R2 in s^-1
DataTypes={"ADC":0,"T1":1,"T2":2}
KTArrayConcentration=[1.0 for i in range(nROIs)]
KTArrayT1=[100.0 for i in range(nROIs)]
nSNRROIs = 1
snrROIdiameter =15
SNRColor = "y"
snrX = -50.
snrY=-50.
snrZ=0.

class KTPhantom(VPhantom.VPhantom):
  """A virtual phantom that describes NIST hcp phantom"""
  def __init__(self):
    VPhantom.VPhantom.__init__(self)
    self.phantomName = "NIST KT Phantom"
    self.ROIsetdict = {'KTArray':0}
    self.KTROIs =self.SetDefaultROIs("KT")    #sets roi positions
    self.KTROIs.Comment = "6 element 30 mm tube array"
    for roi in self.KTROIs.ROIs:    #sets T1, concentration
      roi.T1=KTArrayT1[roi.Index-1]
      roi.ADC=KTArrayT1[roi.Index-1]
      roi.Concentration=KTArrayConcentration[roi.Index-1]
    self.ROIsets.append(self.KTROIs)
    self.SNRROIs=self.SetSNRROIs('SNR')
    self.ROIsets.append(self.SNRROIs)

    
    
  def SetDefaultROIs(self,ptype):
  #set initial system phantom contrast ROIs
      ROIdiameter = tubeDiameter / 2 #set initial ROI radius in mm
      r=VPhantom.ROISet(ptype)
      r.ROIName =  ptype + "Array"
      r.Field = 1.5
      r.Temperature = 20.
      r.nROIs=nROIs

      for I in range (1, nROIs-1): #make ROIs
          r.ROIs.append(VPhantom.ROI())
          r.ROIs[-1].Name = ptype + "-" +str(I)
          r.ROIs[-1].Index = I    #note the index starts at one while the list index starts at 0
          r.ROIs[-1].d1 = ROIdiameter
          r.ROIs[-1].Xcenter = 0.
          r.ROIs[-1].Zcenter = 0.
          r.ROIs[-1].Ycenter = 0
          r.ROIs[-1].Xcenter = arrayRadius*np.sin((I-1)*angle)
          r.ROIs[-1].Ycenter = arrayRadius*np.cos((I-1)*angle)
  #manually set ROI positions    
     
      r.ROIs.append(VPhantom.ROI())
      r.ROIs[-1].Xcenter = -9
      r.ROIs[-1].Ycenter = 0.
      r.ROIs[-1].d1=innerROIDiameter
      r.ROIs.append(VPhantom.ROI())
      r.ROIs[-1].Xcenter =11
      r.ROIs[-1].Ycenter = 0.
      r.ROIs[-1].d1=innerROIDiameter
      
      return r
    
  def SetSNRROIs(self,ptype):
  #set initial system phantom contrast ROIs
      r=VPhantom.ROISet(ptype)
      r.ROIName =  ptype + "Array"
      r.Field = 1.5
      r.Temperature = 20.
      r.nROIs=nSNRROIs
      r.ROIColor = "y"
      r.ROIs.append(VPhantom.ROI())
      r.ROIs[-1].Name = ptype + "-" + "1"
      r.ROIs[-1].Index = 1    #note the index starts at one while the list index starts at 0
      r.ROIs[-1].d1 = snrROIdiameter
      r.ROIs[-1].Xcenter = snrX
      r.ROIs[-1].Zcenter = snrZ
      r.ROIs[-1].Ycenter = snrY
      return r
