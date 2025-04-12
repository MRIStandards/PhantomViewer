'''
Created on Feb 25, 2015
MSK t1rho virtual phantom
@author: stephen russek
'''

# -*- coding: utf-8 -*-

import numpy as np
import VPhantom 
nROIs = 12    # number of ROIs
angle = 60*np.pi/180
a = 24       #lattice constant
PhantomRadius = 60. #radius of phantom in mm
sphereDiameter = 12. #tube diameter
WaterR1 = 0.3      #water R1 in s^-1
WaterR2 = 0.5      #water R2 in s^-1
DataTypes={"ADC":0,"T1":1,"T2":2}
T1ArrayConcentration=[0.299,0.620,1.074,1.716,2.624,3.908,5.724,8.293,11.93,17.06,24.33,34.60,49.13,69.67]
T1ArrayT1=[2047.7,1448.1,1024.,724.1,512.,362.1,256.,181.,128.,90.5,64.,45.3,32.,22.6]
T2ArrayConcentration=[0.01206,0.01989,0.03076,0.04663,0.06879,0.1001,0.14442,0.20705,0.29565,0.42094,0.59814,0.84873,1.20313,1.70431]
T2ArrayT2=[724.1,512.,364.,256.,181.,128.,90.5,64.,45.3,32.,22.6,16.,11.3,8.0]
PDArrayPD=[5.,10.,15.,20.,25.,30.,35.,40.,50.,60.,70.,80.,90.,100.]
nSNRROIs = 1
snrROIdiameter =10
SNRColor = "y"
snrX = -50.
snrY=0.
snrZ=-50.

class MSKPhantom(VPhantom.VPhantom):
  """A virtual phantom that describes MSK phantom"""
  def __init__(self):
    VPhantom.VPhantom.__init__(self)
    self.phantomName = "MSK Phantom"
    self.ROIsetdict = {'hcpArray':0}
    self.hcpROIs =self.SetDefaultROIs("hcp")    #sets roi positions
    self.hcpROIs.Comment = "contrast array"
    for roi in self.hcpROIs.ROIs:    #sets T1, concentration
      roi.T1=T1ArrayT1[roi.Index-1]
      roi.Concentration=T1ArrayConcentration[roi.Index-1]
    self.ROIsets.append(self.hcpROIs)
    self.SNRROIs=self.SetSNRROIs('SNR')
    self.ROIsets.append(self.SNRROIs)
    
    
  def SetDefaultROIs(self,ptype):
  #set initial system phantom contrast ROIs
      ROIdiameter = sphereDiameter/2  #set initial ROI radius in mm
      r=VPhantom.ROISet(ptype)
      r.ROIName =  ptype + "Array"
      r.Field = 1.5
      r.Temperature = 0.
      r.nROIs=nROIs

      for I in range (1, nROIs +1): #make ROIs
          r.ROIs.append(VPhantom.ROI())
          r.ROIs[-1].Name = ptype + "-" +str(I)
          r.ROIs[-1].Index = I    #note the index starts at one while the list index starts at 0
          r.ROIs[-1].d1 = ROIdiameter
          r.ROIs[-1].Xcenter = 0.
          r.ROIs[-1].Zcenter = 0.
          r.ROIs[-1].Ycenter = 0
          r0=-1.5*a*np.array([0,1]) #position of first ROI
          rd=a*np.array([0,1])    #vector from first to second ROI
          rx=a*np.array([-np.sin(angle),np.cos(angle)])
          rnx=a*np.array([np.sin(angle),np.cos(angle)])
  #manually set ROI positions    
      r.ROIs[0].Xcenter = r0[0]
      r.ROIs[0].Zcenter = r0[1]
      r1=r0+rx
      r.ROIs[1].Xcenter = r1[0]
      r.ROIs[1].Zcenter = r1[1]
      r2=r1+rx
      r.ROIs[2].Xcenter = r2[0]
      r.ROIs[2].Zcenter = r2[1]
      r3=r2+rd
      r.ROIs[3].Xcenter = r3[0]
      r.ROIs[3].Zcenter = r3[1]
      r4=r3+rnx
      r.ROIs[4].Xcenter = r4[0]
      r.ROIs[4].Zcenter = r4[1]
      r5=r4+rnx
      r.ROIs[5].Xcenter = r5[0]
      r.ROIs[5].Zcenter = r5[1]      
      r6=r5-rx
      r.ROIs[6].Xcenter = r6[0]
      r.ROIs[6].Zcenter = r6[1]
      r7=r6-rx
      r.ROIs[7].Xcenter = r7[0]
      r.ROIs[7].Zcenter = r7[1] 
      r8=r7-rd
      r.ROIs[8].Xcenter = r8[0]
      r.ROIs[8].Zcenter = r8[1]
      r9=r8-rnx
      r.ROIs[9].Xcenter = r9[0]
      r.ROIs[9].Zcenter = r9[1]
      r10=r9+rx
      r.ROIs[10].Xcenter = r10[0]
      r.ROIs[10].Zcenter = r10[1] 
      r11=r10+rx
      r.ROIs[11].Xcenter = r11[0]
      r.ROIs[11].Zcenter = r11[1]     
      r12=r11+rnx
      r.ROIs[12].Xcenter = r12[0]
      r.ROIs[12].Zcenter = r12[1]     
      r13=r12-rx
      r.ROIs[13].Xcenter = r13[0]
      r.ROIs[13].Zcenter = r13[1]                                
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
