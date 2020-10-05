# -*- coding: utf-8 -*-
"""
Created on Wed Jan 01 10:59:53 2014
module containing parameters describing NIST/ISMRM MRI system phantom
@author: stephen russek
"""
import numpy as np
import VPhantom
import Image 

ArrayRadius = 50.       #contrast arrays have 10 elements equally spaced on a circle or radius array radius in mm
ArrayInnerSpacing = 20. #4 inner contrast spheres are on a square of size ArrayInnerSpacing in mm
nSPArray = 2       #number of ROIs in slice profile arrays
nCArrayMax = 14   #maximum number of ROIs in T1,T2, and PD contrast arrays
nFAArrayMax = 57
nGradDirectionsMax = 12    #maximum number of gradient directions for diffusion imaging
nCArray = nCArrayMax
nT1Array = nCArrayMax      #actual number in T1 array
nT2Array = nCArrayMax
nPDArray = nCArrayMax
nFAArray = nFAArrayMax
nSNRROI = 1
SNRROIxi = 80.
SNRROIyi = 80.
SNRROIzi = 80.
SNRROIri = 10
T1ArrayY = -56.5        #y position of arrays in mm
T2ArrayY = -16.5
PDArrayY = 23.5
GaussiantoRectGaussian = 0.66
FASpacing = 40.     # spacing between fiducial array elements in mm
PhantomRadius = 100. #radius of phantom in mm
FASphereRadius = 5. #inner radius of fiducial spheres in mm
ContrastSphereDiameter = 20. #outside diameter of contrast spheres in mm
WaterR1 = 0.3      #water R1 in s^-1
WaterR2 = 0.5      #water R2 in s^-1
P0NiCl2 = 0.63738
P1NiCl2 = -0.01835
P2NiCl2 = 0.00894
P0MnCl2 = 32.1
P1MnCl2 = 27.3
DataTypes={"FA":0,"T1":1,"T2":2,"PD":3,"SP":4,"RI":5, "SNR":6}
spXCenter1 = 14.
spYCenter1 = -32.
spZCenter1 = 0.
spXCenter2 = 22.
spYCenter2 = -32.
spZCenter2 = 0.
spWidth = 4.        #width of slice profile ROIs in mm
spLength = 60.      #Length of slice profile ROI in mm ***has to be an integer for averaging and binning**
SPangle = 10. * np.pi / 180.
T1ArrayConcentration=[0.299,0.620,1.074,1.716,2.624,3.908,5.724,8.293,11.93,17.06,24.33,34.60,49.13,69.67]
T1ArrayT1=[2047.7,1448.1,1024.,724.1,512.,362.1,256.,181.,128.,90.5,64.,45.3,32.,22.6]
T2ArrayConcentration=[0.01206,0.01989,0.03076,0.04663,0.06879,0.1001,0.14442,0.20705,0.29565,0.42094,0.59814,0.84873,1.20313,1.70431]
T2ArrayT2=[724.1,512.,364.,256.,181.,128.,90.5,64.,45.3,32.,22.6,16.,11.3,8.0]
PDArrayPD=[5.,10.,15.,20.,25.,30.,35.,40.,50.,60.,70.,80.,90.,100.]
systemPhantomImage = Image.open("MRI_System_Phantom.jpg")

#  def SetDefaultFiducialROIs(self):
#    for I in range (-2,3):
#      for j in range(-2,3):
#        for K in range(-2,3):
#          If (I ^ 2 + j ^ 2 + K ^ 2) ^ 0.5 * ROIProperties.FASpacing < (PhantomRadius - FASphereRadius - 2):
#            FAROIsi(nf).ROIRadius = FASphereRadius
#            FAROIsi(nf).Xcenter = I * FASpacing
#            FAROIsi(nf).Ycenter = j * FASpacing
#            FAROIsi(nf).Zcenter = K * FASpacing
#            nf = nf + 1
#    for I in range ( 1, nFAArray)
#        FAROIs(I) = FAROIsi(I)

class SystemPhantom(VPhantom.VPhantom):
  """A virtual phantom that describes NIST ISMRM MRI system phantom"""
  def __init__(self):
    VPhantom.VPhantom.__init__(self)
    self.phantomName = "NIST-ISMRM MRI System Phantom"
    self.phantomImage=systemPhantomImage
    self.ROIsetdict = {'T1Array':0, 'T2Array':1, 'PDArray':2,'SPArray':3, 'ResolutionInset':4}
    
    self.T1ROIs =self.SetDefaultContrastROIs("T1")    #sets roi positions
    self.T1ROIs.Comment = "NiCl2 array with root 2 decrease in T1"
    for roi in self.T1ROIs.ROIs:    #sets T1, concentration
      roi.T1=T1ArrayT1[roi.Index-1]
      roi.Concentration=T1ArrayConcentration[roi.Index-1]
    self.ROIsets.append(self.T1ROIs)
    
    self.T2ROIs =self.SetDefaultContrastROIs("T2")
    self.T2ROIs.Comment = "MnCl2 array with root 2 decrease in T2" 
    for roi in self.T2ROIs.ROIs:    #sets T2, concentration
      roi.T2=T2ArrayT2[roi.Index-1]
      roi.Concentration=T2ArrayConcentration[roi.Index-1]
    self.ROIsets.append (self.T2ROIs)
    
    self.PDROIs =self.SetDefaultContrastROIs("PD")    #sets roi positions
    self.PDROIs.Comment = "D20 array 0% to 95%"
    for roi in self.PDROIs.ROIs:    #sets PD, concentration
      roi.PD=PDArrayPD[roi.Index-1]
    self.ROIsets.append(self.PDROIs)
    
  def SetDefaultContrastROIs(self,ptype):
  #set initial system phantom contrast ROIs
      ROIdiameter = ContrastSphereDiameter/2 #set initial ROI radius in mm
      Phase = 0
      ArrayY=T1ArrayY
      nArray=nCArrayMax
      r=VPhantom.ROISet(ptype)
      r.ROIName =  ptype + "Array"
      r.Field = 1.5
      r.Temperature = 20.
      r.nROIs=14

      for I in range (1, nArray - 3): #set properties of 10 element circular array
          r.ROIs.append(VPhantom.ROI())
          r.ROIs[-1].Name = ptype + "-" +str(I)
          r.ROIs[-1].Index = I    #note the index starts at one while the list index starts at 0
          r.ROIs[-1].d1 = ROIdiameter
          r.ROIs[-1].Xcenter = ArrayRadius * np.sin(2 * np.pi * (I - 1) / (nArray - 4) + Phase)
          r.ROIs[-1].Zcenter = ArrayRadius * np.cos(2 * np.pi * (I - 1) / (nArray - 4) + Phase)
          r.ROIs[-1].Ycenter = ArrayY
  
      r.ROIs.append(VPhantom.ROI())
      r.ROIs[-1].Name = ptype + "-" +str(nArray - 3)
      r.ROIs[-1].Index = nArray - 3
      r.ROIs[-1].d1 = ROIdiameter
      r.ROIs[-1].Xcenter = -ArrayInnerSpacing
      r.ROIs[-1].Zcenter = ArrayInnerSpacing
      r.ROIs[-1].Ycenter = ArrayY
  
      r.ROIs.append(VPhantom.ROI())
      r.ROIs[-1].Name = ptype + "-" +str(nArray - 2)
      r.ROIs[-1].Index = nArray - 2
      r.ROIs[-1].d1 = ROIdiameter
      r.ROIs[-1].Xcenter = ArrayInnerSpacing
      r.ROIs[-1].Zcenter = ArrayInnerSpacing
      r.ROIs[-1].Ycenter = ArrayY
  
      r.ROIs.append(VPhantom.ROI())
      r.ROIs[-1].Name = ptype + "-" +str(nArray - 1)
      r.ROIs[-1].Index = nArray - 1
      r.ROIs[-1].d1 = ROIdiameter
      r.ROIs[-1].Xcenter = ArrayInnerSpacing
      r.ROIs[-1].Zcenter = -ArrayInnerSpacing
      r.ROIs[-1].Ycenter = ArrayY
      
      r.ROIs.append(VPhantom.ROI())
      r.ROIs[-1].Name = ptype + "-" +str(nArray)
      r.ROIs[-1].Index = nArray 
      r.ROIs[-1].d1 = ROIdiameter
      r.ROIs[-1].Xcenter = -ArrayInnerSpacing
      r.ROIs[-1].Zcenter = -ArrayInnerSpacing
      r.ROIs[-1].Ycenter = ArrayY
    
      return r
      
  
  #  def SetDefaultSliceProfileROIs(self):    
  #    'set slice profile ROIs
  #    SPROIsi(1).Xcenter = spXCenter1
  #    SPROIsi(1).Ycenter = spYCenter1
  #    SPROIsi(1).Zcenter = spZCenter1
  #    SPROIsi(2).Xcenter = spXCenter2
  #    SPROIsi(2).Ycenter = spYCenter2
  #    SPROIsi(2).Zcenter = spZCenter2
  #    
  #    SPROIs(1) = SPROIsi(1)
  #    SPROIs(2) = SPROIsi(2)
  #    
  #    SNRROI(1).ROIRadius = SNRROIri
  #    SNRROI(1).Xcenter = SNRROIxi
  #    SNRROI(1).Zcenter = SNRROIyi
  #    SNRROI(1).Ycenter = SNRROIzi