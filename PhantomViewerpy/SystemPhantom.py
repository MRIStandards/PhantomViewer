# -*- coding: utf-8 -*-
"""
Created on Wed Jan 01 10:59:53 2014
module containing parameters describing NIST/ISMRM MRI system phantom
@author: stephen russek
"""
import numpy as np
import VPhantom


ArrayRadius = 50.       #contrast arrays have 10 elements equally spaced on a circle or radius array radius in mm
ArrayInnerSpacing = 20. #4 inner contrast spheres are on a square of size ArrayInnerSpacing in mm
nSPArray = 2       #number of ROIs in slice profile arrays
nCArrayMax = 14   #maximum number of ROIs in T1,T2, and PD contrast arrays
nFAArrayMax = 57
nGradDirectionsMax = 12    #maximum number of gradient directions for diffusion imaging
nCArray = nCArrayMax
nNiArray = nCArrayMax      #actual number in T1 array
nMnArray = nCArrayMax
nPDArray = nCArrayMax
nFAArray = nFAArrayMax
nSNRROIs = 2    #ROIs to determine background counts and noise
bgX = 80.
bgY = 80.
bgZ = 80.
snrX = 0.
snrY = 0.
snrZ = 0.
snrROIdiameter = 30
bgROIdiameter = 20
NiArrayY = -56.5        #y position of arrays in mm
MnArrayY = -16.5
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
DataTypes={"FA":0,"T1":1,"T2":2,"PD":3,"SP":4,"RI":5, "SNR":6,"LCTherm":6}
spXCenter1 = 14.    #nominal slice profile wedge locations
spYCenter1 = -32.
spZCenter1 = 0.
spXCenter2 = 22.
spYCenter2 = -32.
spZCenter2 = 0.
spWidth = 4.        #width of slice profile ROIs in mm
spLength = 54.      #Length of slice profile ROI in mm 
SPangle = 10. * np.pi / 180.
commercialSystemPhantom=True  #commercial phantom flag 
phantomReferenceData='Keenan ISMRM 3290 2016 @20C' 
NiArrayConcentration=[0.299, 0.623,1.072,1.72,2.617,3.912,5.731,8.297,11.936,17.07,24.326,34.59,49.122,69.68]
NiArrayT1_1p5T=[2033,1489,1012,730.8,514.1,367.9,260.1,184.6,132.7,92.7,65.4,46.32,32.45,22.859]
NiArrayT1_3T=[1989,1454,984.1,706,496.7,351.5,247.13,175.3,125.9,89,62.7,44.53,30.84,21.719]
NiArrayT2_1p5T=[1669,1244,859.3,628.5,446.3,321.2,227.7,161.9,117.1,81.9,57.7,41,28.7,20.2]
NiArrayT2_3T=[1465,1076,717.9,510.1,359.6,255.5,180.8,127.3,90.3,64.3,45.7,31.86,22.38, 15.83]
MnArrayConcentration=[0.013,0.021,0.031,0.047,0.069,0.101,0.145,0.207,0.296,0.421,0.599,0.849,1.104,1.704]
MnArrayT1_1p5T=[2376,2183,1870,1539,1237,1030,752.2,550.2,413.4,292.9,194.9,160.2,106.4,83.33]
MnArrayT1_3T=[2480,2173,1907,1604,1332,1044,801.7,608.6,458.4,336.5,244.2,176.6,126.9,90.9]
MnArrayT2_1p5T=[939.4,594.3,416.5,267,184.9,140.6,91.76,64.84,45.28,30.62,19.76,15.99,10.47,8.15]
MnArrayT2_3T=[581.3,403.5,278.1,190.94,133.27,96.89,64.07,46.42,31.97,22.56,15.813,11.237,7.911,5.592]
# NiArrayConcentration=[0.299,0.620,1.074,1.716,2.624,3.908,5.724,8.293,11.93,17.06,24.33,34.60,49.13,69.67]
# NiArrayT1=[2047.7,1448.1,1024.,724.1,512.,362.1,256.,181.,128.,90.5,64.,45.3,32.,22.6]
# MnArrayConcentration=[0.01206,0.01989,0.03076,0.04663,0.06879,0.1001,0.14442,0.20705,0.29565,0.42094,0.59814,0.84873,1.20313,1.70431]
# MnArrayT2=[724.1,512.,364.,256.,181.,128.,90.5,64.,45.3,32.,22.6,16.,11.3,8.0]
PDArrayPD=[5.,10.,15.,20.,25.,30.,35.,40.,50.,60.,70.,80.,90.,100.]

comment = "http://collaborate.nist.gov/mriphantoms/bin/view/MriPhantoms/MRISystemPhantom"
fiducialT1=100.0
fiducialT2=50.0
fiducialConcentration=1.0
fidColor1={'center':(102, 102, 102), 'left':(0, 114, 189), 'right':(217, 83, 25),
           'anterior':(237, 177,32 ), 'posterior':(126, 47, 32),
           'superior':(119, 47, 224), 'inferior':(155, 232, 234)}  # KKKolors
fidColor2={'center':(255, 165, 0), 'left':'r', 'right':'g',
           'anterior':'b', 'posterior':'c',
           'superior':'m', 'inferior':'y'} # orgbcmy
fidColor=fidColor1 
fidSymbol={'center':'o', 'left':'t', 'right':'s',
           'anterior':'p', 'posterior':'h',
           'superior':'star', 'inferior':'d'} 
#resolution inset parameters
RI1UR = np.array([[0.0,0.0],[10.58,.0],[19.76,0],[27.56,0], [33.96,0]] )    # upper right coordinate in mm
RI1d=np.array([0.8,0.7,0.6,0.5,0.4])  #hole diameters in mm
RI1angle=10.0   #angle of angled array in degrees
RInholes=155    #number of holes in ACR resolution inset
RIfrw=53.5      #resolution inset frame width
RIfrh=15.5        #resolution inset frame height
RIframeLL=np.array([-RIfrw/2,-RIfrh/2])
RIsize=np.array([RIfrw,RIfrh])
RIxoffset=15.0
RIzoffset=-0.5
#thermometer
TempArray=[15.0, 16.0,17.0,18.0,19.0,20.0,21.0,22.0,23.0,24.0]  #transition temperatures of 10-element liquid crystal thermometer
ThermCellDiameter=5   #5mm cells
ThermCellX=[24.3,15.3,10.2,5.1,0,-5.1,-10.2,-15.3,-20.4,-25.5]  #Position of the liquid crystal thermometer cells
ThermCellZ=[-4,4,-4,4,-4,4,-4,4,-4,4]



class SystemPhantom(VPhantom.VPhantom):
  """A virtual phantom that describes NIST ISMRM MRI system phantom"""
  def __init__(self):
    VPhantom.VPhantom.__init__(self)
    self.phantomName = "NIST-ISMRM MRI System Phantom"
    self.phantomReferenceData=phantomReferenceData
    self.Comment = comment
    self.phantomImage = "..\icons\MRISystemPhantom.jpg"
    self.ROIsetdict = {'NiArray':0, 'MnArray':1, 'PDArray':2,'FiducialArray':3, 'ResolutionInset':4}
    
    self.NiROIs =self.SetDefaultContrastROIs("Ni")    #sets roi positions
    self.NiROIs.Comment = "NiCl2 array with root 2 decrease in T1"
    self.NiROIs.ROIColor = 'g'
#     NiArrayT1=NiArrayT1_1p5T
#     NiArrayT2=NiArrayT2_1p5T
#     MnArrayT1=MnArrayT1_1p5T
#     MnArrayT2=MnArrayT2_1p5T    
    NiArrayT1=NiArrayT1_3T
    NiArrayT2=NiArrayT2_3T
    MnArrayT1=MnArrayT1_3T
    MnArrayT2=MnArrayT2_3T        
    for roi in self.NiROIs.ROIs:    #sets T1, concentration
      roi.T1=NiArrayT1[roi.Index-1]
      roi.T2=NiArrayT2[roi.Index-1]
      roi.Concentration=NiArrayConcentration[roi.Index-1]
    self.ROIsets.append(self.NiROIs)
    
    self.MnROIs =self.SetDefaultContrastROIs("Mn")
    self.MnROIs.Comment = "MnCl2 array with root 2 decrease in T2" 
    self.MnROIs.ROIColor = 'r'
    self.MnROIs.showBackgroundROI=False    #flag to show background ROI to determine if an ROI should be discarded
    for roi in self.MnROIs.ROIs:    #sets T2, concentration
      roi.T1=MnArrayT1[roi.Index-1]
      roi.T2=MnArrayT2[roi.Index-1]
      roi.Concentration=MnArrayConcentration[roi.Index-1]
    self.ROIsets.append (self.MnROIs)
    
    self.PDROIs =self.SetDefaultContrastROIs("PD")    #sets roi positions
    self.PDROIs.Comment = "D20 array 0% to 95%"
    self.PDROIs.ROIColor = 'y'
    self.PDROIs.showBackgroundROI=True    #flag to show background ROI to determine if an ROI should be discarded
    self.PDROIs.showSNRROI=True           #flag to show a noise ROI to determine noise from an image subtraction
    for roi in self.PDROIs.ROIs:    #sets PD, concentration
      roi.PD=PDArrayPD[roi.Index-1]
    self.ROIsets.append(self.PDROIs)
    
    self.SNRROIs=self.SetSNRROIs('SNR')
    
    self.fiducialROIs=self.SetDefaultFiducialROIs('Fiducial')
    self.fiducialROIs.Comment = "57 x 1cm spheres on 4cm grid"
    for roi in self.fiducialROIs.ROIs:    #sets T2, concentration
      roi.T1=fiducialT1
      roi.T2=fiducialT2
      roi.Concentration=fiducialConcentration
    self.ROIsets.append (self.fiducialROIs)
    
    self.RIROIs=self.SetRIROIs('RI')
    self.ROIsets.append (self.RIROIs)

    self.SPROIs=self.SetSliceProfileROIs('SP')
    self.ROIsets.append (self.SPROIs)
    
    
    self.TempROIs=self.SetThermometerROIs('LCTherm')
    self.ROIsets.append (self.TempROIs)
        
  def SetDefaultContrastROIs(self,ptype):
  #set initial system phantom contrast ROIs
      ROIdiameter = ContrastSphereDiameter/2 #set initial ROI radius in mm
      Phase = 0
      ArrayY=NiArrayY
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
      
  def SetSNRROIs(self,ptype):
  #set Baseline ROIs
      r=VPhantom.ROISet(ptype)
      r.ROIName =  ptype + "Array"
      r.Field = 1.5
      r.Temperature = 20.
      r.nROIs=nSNRROIs
      r.ROIColor = "b"
#       r.ROIs.append(VPhantom.ROI())   #ROI in non signal region to determine background signal used to exclude data near background
#       r.ROIs[-1].Name = ptype + "-" + "1"
#       r.ROIs[-1].Index = 1    #note the index starts at one while the list index starts at 0
#       r.ROIs[-1].d1 = bgROIdiameter
#       r.ROIs[-1].Xcenter = bgX
#       r.ROIs[-1].Zcenter = bgZ
#       r.ROIs[-1].Ycenter = bgY
      r.ROIs.append(VPhantom.ROI()) # large ROI in featureless region to determine noise via subtracting identical images
      r.ROIs[-1].Name = ptype + "-" + "1"
      r.ROIs[-1].Index = 1    #note the index starts at one while the list index starts at 0
      r.ROIs[-1].d1 = snrROIdiameter
      r.ROIs[-1].Xcenter = snrX
      r.ROIs[-1].Zcenter = snrZ
      r.ROIs[-1].Ycenter = snrY
      return r
    
  def SetSliceProfileROIs(self, ptype):    
      '''set 2 rectangular slice profile ROIs'''
      r=VPhantom.ROISet(ptype)
      r.ROIName =  ptype + "Array"
      r.Field = 1.5
      r.Temperature = 20.
      r.nROIs=2
      r.ROIColor = "y"
      roi=VPhantom.ROI()    #add a rectangular slice profile ROI
      roi.Name = ptype + "SP1"
      roi.Type='Rectangle'
      roi.dx=spWidth
      roi.dy=spLength
      roi.Xcenter = spXCenter1
      roi.Zcenter = spZCenter1
      roi.Ycenter = spYCenter1
      roi.Label = 'nolabel'
      roi.color='r'
      r.ROIs.append(roi)
      roi=VPhantom.ROI()    #add a rectangular slice profile ROI
      roi.Name = ptype + "SP21"
      roi.Type='Rectangle'
      roi.dx=spWidth
      roi.dy=spLength
      roi.Xcenter = spXCenter2
      roi.Zcenter = spZCenter2
      roi.Ycenter = spYCenter2
      roi.Label = 'nolabel'
      roi.color='b'
      r.ROIs.append(roi)
      r.nROIs=len(r.ROIs)
      return r
  
  def SetDefaultFiducialROIs(self, ptype):
      '''Defines 57 fiducial spheres, step first in x (LR), then z(SI), then, y(AP)
       from Plate1 to Plate5; commercial system phantoms do not have FId. 57 which is next to fill port'''
      r=VPhantom.ROISet(ptype)
      #ri=VPhantom.ROISet(ptype)
      r.ROIName =  ptype + "Array"
      r.Field = 1.5
      r.Temperature = 20.

      r.ROIColor = "m"
      nindex=1
      for J in range (-2,3):
          for K in range(-2,3):
              for I in range(-2,3):
                  if (I**2 + J**2 + K**2)**0.5 * FASpacing < (PhantomRadius - FASphereRadius - 2):
                       r.ROIs.append(VPhantom.ROI())
                       r.ROIs[-1].Name = ptype + "-" +str(I) + ',' + str(J) + ',' + str(K)
                       r.ROIs[-1].Index = nindex
                       r.ROIs[-1].Label = str(nindex)
                       r.ROIs[-1].d1 = FASphereRadius*2.0
                       r.ROIs[-1].Xcenter = I * FASpacing
                       r.ROIs[-1].Zcenter = K * FASpacing
                       r.ROIs[-1].Ycenter = -J * FASpacing
                       #set identical initial ROIs which do not change
                       r.initialROIs.append(VPhantom.ROI())
                       r.initialROIs[-1].Name = ptype + "-" +str(I) + ',' + str(J) + ',' + str(K)
                       r.initialROIs[-1].Index = nindex
                       r.initialROIs[-1].Label = str(nindex)
                       r.initialROIs[-1].d1 = FASphereRadius*2.0
                       r.initialROIs[-1].Xcenter = I * FASpacing
                       r.initialROIs[-1].Zcenter = K * FASpacing
                       r.initialROIs[-1].Ycenter = -J * FASpacing
                       r.ROIs[-1].color=fidColor['center']
                       r.ROIs[-1].symbol=fidSymbol['center']
                       if r.ROIs[-1].Xcenter==80: #Left
                         r.ROIs[-1].color=fidColor['left']
                         r.ROIs[-1].symbol=fidSymbol['left']
                       if r.ROIs[-1].Xcenter==-80:  #right
                         r.ROIs[-1].color=fidColor['right']
                         r.ROIs[-1].symbol=fidSymbol['right']
                       if r.ROIs[-1].Ycenter==80:   #Posterior
                         r.ROIs[-1].color=fidColor['posterior']
                         r.ROIs[-1].symbol=fidSymbol['posterior']
                       if r.ROIs[-1].Ycenter==-80:  #Anterior
                         r.ROIs[-1].color=fidColor['anterior']
                         r.ROIs[-1].symbol=fidSymbol['anterior']
                       if r.ROIs[-1].Zcenter==80:   #Superior
                         r.ROIs[-1].color=fidColor['superior']
                         r.ROIs[-1].symbol=fidSymbol['superior']
                       if r.ROIs[-1].Zcenter==-80:  #Inferior
                         r.ROIs[-1].color=fidColor['inferior']
                         r.ROIs[-1].symbol=fidSymbol['inferior']
                         
                       nindex +=1
      if commercialSystemPhantom:
        del r.ROIs[54]
        del r.initialROIs[54]
      r.nROIs=len(r.ROIs)
      return r
    
  def SetRIROIs(self, ptype):
      '''reads in resolution set hole positions and defines an ROI set'''
      riholes=self.createResArray()   #make an array of hole positions
      r=VPhantom.ROISet(ptype)
      r.ROIName =  ptype + "Array"
      r.Field = 1.5
      r.Temperature = 20.
      r.ROIColor = 'y'
      for i, hole in enumerate(riholes):    #make array of hole ROIs in the XZ plane
        roi=VPhantom.ROI()
        roi.Name = ptype + "-" +str(i) 
        roi.Index = i+1
        roi.d1 = hole[2]
        roi.Xcenter = hole[0]
        roi.Zcenter = hole[1]
        roi.Ycenter = 0
        roi.Label = 'nolabel'
        r.ROIs.append(roi)
      roi=VPhantom.ROI()    #add a rectangular array that enscribes RI coffin
      roi.Name = ptype + "BoundingBox"
      roi.Type='Rectangle'
      roi.dx=RIfrw
      roi.dy=RIfrh
      roi.Xcenter = 0
      roi.Zcenter = 0
      roi.Ycenter = 0
      roi.Label = 'nolabel'
      r.ROIs.append(roi)
      r.nROIs=len(r.ROIs)
      return r

  def SetThermometerROIs(self, ptype):
      '''Defines 10 fliquid crytal thermometer cells'''
      r=VPhantom.ROISet(ptype)
      #ri=VPhantom.ROISet(ptype)
      r.ROIName =  ptype + "Array"
      r.Field = 1.5
      r.Temperature = 20.
      for i in range (10):
        r.ROIs.append(VPhantom.ROI())
        r.ROIs[-1].Name = ptype + "-" +str(TempArray[i])
        r.ROIs[-1].Index = i+1
        r.ROIs[-1].Label = str(TempArray[i])
        r.ROIs[-1].Value = TempArray[i]
        r.ROIs[-1].d1 = ThermCellDiameter
        r.ROIs[-1].Xcenter = ThermCellX[i]
        r.ROIs[-1].Zcenter = ThermCellZ[i]
        r.ROIs[-1].Ycenter = 0
      r.nROIs=len(r.ROIs)
      return r
        
  def createResArray(self):
    '''creates a hole array consisting of five 4x4+rotated(4x4) holes of different sizes'''
    holearray=[]  #array of holes (x,y,diameter) all in mm
    for n,d in enumerate(RI1d):
      for i, j in np.ndindex((4,4)):
        holearray.append(np.array([-i*2*d+RI1UR[n,0]-RIxoffset,-(j*2*d+RI1UR[n,1]-RIzoffset),d]))
      for i, j in np.ndindex((4,4)):
        x,y=self.rotateArray(i*2*d,j*2*d,RI1angle+90.0)
        if not(j==0 and i==0):
          holearray.append(np.array([x+RI1UR[n,0]-RIxoffset,-(y+RI1UR[n,1]-RIzoffset),d]))
    return holearray
  
  def rotateArray(self, x,y,angle): 
    theta=angle*np.pi/180.0
    xr=np.cos(theta)*x+np.sin(theta)*y
    yr=-np.sin(theta)*x+np.cos(theta)*y
    return xr,yr 
    
  def changeB0(self,b0=1.5):
    self.B0=b0
    self.NiROIs.field=b0
    self.MnROIs.field=b0       
    if b0 > 1.4 and b0 < 1.6 :
        NiArrayT1=NiArrayT1_1p5T
        NiArrayT2=NiArrayT2_1p5T
        MnArrayT1=MnArrayT1_1p5T
        MnArrayT2=MnArrayT2_1p5T
   
    if b0 > 2.9 and b0 < 3.1:
        NiArrayT1=NiArrayT1_3T
        NiArrayT2=NiArrayT2_3T
        MnArrayT1=MnArrayT1_3T
        MnArrayT2=MnArrayT2_3T        
    for roi in self.NiROIs.ROIs:    #sets T1, concentration
      roi.T1=NiArrayT1[roi.Index-1]
      roi.T2=NiArrayT2[roi.Index-1]
    
    for roi in self.MnROIs.ROIs:    #sets T2, concentration
      roi.T1=MnArrayT1[roi.Index-1]
      roi.T2=MnArrayT2[roi.Index-1]
                      