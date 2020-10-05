'''
Created on Jan 23, 2015

@author: stephen russek
'''
# -*- coding: utf-8 -*-

import numpy as np
import VPhantom 
nROIs = 17    # number of ROIs
angle = 30*np.pi/180
a = 21.2       #lattice constant
PhantomRadius = 60. #radius of phantom in mm
tubeDiameter = 15. #tube diameter
WaterR1 = 0.3      #water R1 in s^-1
WaterR2 = 0.5      #water R2 in s^-1
DataTypes={"ADC":0,"T1":1,"T2":2}
hcpArrayConcentration=[0,0,0,10,10,20,20,30,30,40,40,50,50,0,0,0,0]
hcpArrayT1=[100.0 for i in range(17)]
hcpArrayADC=[1.0 for i in range(17)]    #ADC in 10^-3 mm2/s
nSNRROIs = 1
snrROIdiameter =20
SNRColor = "y"
snrX = -50.
snrY=-50.
snrZ=0.

class hcpPhantom(VPhantom.VPhantom):
  """A virtual phantom that describes NIST hcp phantom"""
  def __init__(self):
    VPhantom.VPhantom.__init__(self)
    self.phantomName = "NIST hcp Phantom"
    self.ROIsetdict = {'hcpArray':0}
    self.hcpROIs =self.SetDefaultROIs("hcp")    #sets roi positions
    self.hcpROIs.Comment = "PVP array"
    for roi in self.hcpROIs.ROIs:    #sets T1, concentration
      roi.T1=hcpArrayT1[roi.Index-1]
      roi.ADC=hcpArrayADC[roi.Index-1]
      roi.Concentration=hcpArrayConcentration[roi.Index-1]
    self.ROIsets.append(self.hcpROIs)
    self.SNRROIs=self.SetSNRROIs('SNR')
    self.ROIsets.append(self.SNRROIs)
    
    
  def SetDefaultROIs(self,ptype):
  #set initial system phantom contrast ROIs
      ROIdiameter = tubeDiameter / 2 #set initial ROI radius in mm
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
          r0=2*a*np.array([np.sin(angle),np.cos(angle)])
          rx=a*np.array([1,0])
          rd=a*np.array([np.sin(angle),np.cos(angle)])
  #manually set ROI positions    
      r.ROIs[0].Xcenter = r0[0]
      r.ROIs[0].Ycenter = r0[1]
      r1=r0-rx
      r.ROIs[1].Xcenter = r1[0]
      r.ROIs[1].Ycenter = r1[1]
      r2=r1-rx
      r.ROIs[2].Xcenter = r2[0]
      r.ROIs[2].Ycenter = r2[1]
      r3=r2-rd
      r.ROIs[3].Xcenter = r3[0]
      r.ROIs[3].Ycenter = r3[1]
      r4=r3+rx
      r.ROIs[4].Xcenter = r4[0]
      r.ROIs[4].Ycenter = r4[1]
      r5=r4+rx
      r.ROIs[5].Xcenter = r5[0]
      r.ROIs[5].Ycenter = r5[1]      
      r6=r5+rx
      r.ROIs[6].Xcenter = r6[0]
      r.ROIs[6].Ycenter = r6[1]
      r7=r6-rd
      r.ROIs[7].Xcenter = r7[0]
      r.ROIs[7].Ycenter = r7[1] 
      r8=r7-rx
      r.ROIs[8].Xcenter = r8[0]
      r.ROIs[8].Ycenter = r8[1]
      r9=r8-rx
      r.ROIs[9].Xcenter = r9[0]
      r.ROIs[9].Ycenter = r9[1]
      r10=r9-rd
      r.ROIs[10].Xcenter = r10[0]
      r.ROIs[10].Ycenter = r10[1] 
      r11=r10+rx
      r.ROIs[11].Xcenter = r11[0]
      r.ROIs[11].Ycenter = r11[1]     
      r12=r11+rx
      r.ROIs[12].Xcenter = r12[0]
      r.ROIs[12].Ycenter = r12[1]     
      r13=r12+rx
      r.ROIs[13].Xcenter = r13[0]
      r.ROIs[13].Ycenter = r13[1]
      r14=r13-rd
      r.ROIs[14].Xcenter = r14[0]
      r.ROIs[14].Ycenter = r14[1]     
      r15=r14-rx
      r.ROIs[15].Xcenter = r15[0]
      r.ROIs[15].Ycenter = r15[1]     
      r16=r15-rx
      r.ROIs[16].Xcenter = r16[0]
      r.ROIs[16].Ycenter = r16[1]                                   
      return r
    
  def SetSNRROIs(self,ptype):
  #set Baseline ROIs
      r=VPhantom.ROISet(ptype)
      r.ROIName =  ptype
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