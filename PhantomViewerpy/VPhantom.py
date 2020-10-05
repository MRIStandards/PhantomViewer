# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 13:11:14 2013
Contains classes that constitute a virtual phantom
A virtual phantom contains elements, ROIs, 
@author: stephen russek
"""
import numpy as np



class VPhantom():
  """Class that defines a virtual phantom with all properties required to simulate image of a real object"""
  def __init__(self):
    self.phantomName = "New phantom"
    self.B0 = 1.5        #B0 field in Tesla
    self.Temperature = 20   #Temperature in C"
    self.Comment = "New phantom"
    self.units = "mm, ms, C, mM"
    self.nROISets = 0
    self.ROIsets= []   #list of ROI sets included in the phantom
    self.ROIsetdict = {}     #dictionary that lists all ROI sets in a particular phantom
    self.Elements = []
    self.phantomImage =  ""    #image of the phantom
    
  def printROIinfo(self):
      """Returns a string describing a virtual phantom"""
      ps = "Phantom Viewer ROI file" + "\n"
      ps = ps + "PhantomName= " + self.phantomName + "\n"
      ps = ps + "Units= distances in mm, times in ms, fields in T, Temperatures in C, concentrations in mM" + "\n"
      ps = ps + "PhantomComment= " + self.Comment + "\n"
      ps = ps + "B0= " + str(self.B0) + "\n"
      ps = ps + "Temperature= " + str(self.Temperature) + "\n" + "\n"
      ps = ps + "NumberofROIsets= " + str(len(self.ROIsets)) + "\n" + "\n"
      for rs in self.ROIsets:
        ps = ps + "ROISetName="  + rs.ROIName + "\n"
        ps = ps + "ROISetComment="  + rs.Comment + "\n"
        ps = ps + "ROIColor="  + rs.ROIColor + "\n"
        ps = ps + "nROIs=" + str(rs.nROIs) + "\n"
        ps = ps + "Index="  + str([o.Index for o in rs.ROIs]).replace('[','').replace(']','') + "\n"
        ps = ps + "Concentration="  + str(["{:.3f}".format(o.Concentration) for o in rs.ROIs]).replace('[','').replace(']','').replace("\'",'') + "\n"
        ps = ps + "T1="  + str(["{:.2f}".format(o.T1) for o in rs.ROIs]).replace('[','').replace(']','').replace("\'",'') + "\n"
        ps = ps + "T2="  + str(["{:.2f}".format(o.T2) for o in rs.ROIs]).replace('[','').replace(']','').replace("\'",'') + "\n"
        ps = ps + "PD="  + str(["{:.2f}".format(o.PD) for o in rs.ROIs]).replace('[','').replace(']','').replace("\'",'') + "\n"
        ps = ps + "Xcenter="  + str(["{:.2f}".format(o.Xcenter) for o in rs.ROIs]).replace('[','').replace(']','').replace("\'",'') + "\n"
        ps = ps + "Ycenter="  + str(["{:.2f}".format(o.Ycenter) for o in rs.ROIs]).replace('[','').replace(']','').replace("\'",'') + "\n"
        ps = ps + "Zcenter="  + str(["{:.2f}".format(o.Zcenter) for o in rs.ROIs]).replace('[','').replace(']','').replace("\'",'') + "\n"
        ps = ps + "\n"
      return ps
    
  def readPhantomFile(self, fileName=''):
    '''Reads in and parses phantom file'''
    f = open(str(fileName), 'r')
    for line in f:
      parameter = line[:line.find('=')]
      values=line[line.find('=')+1:-1]
      if parameter == "PhantomName":
        self.phantomName = values
      if parameter == "B0":
        self.B0 = float(values)
      if parameter == "Temperature":
        self.Temperature = float(values)
      if parameter == "PhantomComment":
        self.Comment = values
      if parameter == "NumberofROIsets":
        self.nROISets = int(values)
      if parameter == "ROISetName":
        newROISet=ROISet(values)    #create new ROI set and append to list of ROISets
        newROISet.ROIName = values
        self.ROIsets.append(newROISet)
        if newROISet.ROIName=='SNR':   #sets flag if phantom has SNR ROIs to determine background
          self.SNRROIs=newROISet
      if parameter == "nROIs":
        newROISet.nROIs = int(values)
        newROISet.ROIs= [ROI() for i in range(newROISet.nROIs)]
      if parameter == "ROISetComment":
        newROISet.Comment = values
      if parameter == "ROIColor":
        newROISet.Color = values        
      if parameter == "Index":
        for i in range (newROISet.nROIs):
          newROISet.ROIs[i].Index = int(values.split(",")[i])                    
      if parameter == "T1":
        for i in range (newROISet.nROIs):
          newROISet.ROIs[i].T1 = float(values.split(",")[i]) 
      if parameter == "T2":
        for i in range (newROISet.nROIs):
          newROISet.ROIs[i].T2= float(values.split(",")[i])
      if parameter == "PD":
        for i in range (newROISet.nROIs):
          newROISet.ROIs[i].PD= float(values.split(",")[i])
      if parameter == "ADC":
        for i in range (newROISet.nROIs):
          newROISet.ROIs[i].ADC= float(values.split(",")[i])                  
      if parameter == "Concentration":
        for i in range (newROISet.nROIs):
          newROISet.ROIs[i].Concentration = float(values.split(",")[i])  
      if parameter == "Xcenter":
        for i in range (newROISet.nROIs):
          newROISet.ROIs[i].Xcenter = float(values.split(",")[i])  
      if parameter == "Ycenter":
        for i in range (newROISet.nROIs):
          newROISet.ROIs[i].Ycenter = float(values.split(",")[i])
      if parameter == "Zcenter":
        for i in range (newROISet.nROIs):
          newROISet.ROIs[i].Zcenter = float(values.split(",")[i])          

  def changeB0(self,b0=1.5):
          self.B0=b0
        
class Element():
    pass

class Material():
    pass

class ROISet():
  """Set of ROIs"""
  def __init__(self, name):
    self.ROIName  = name
    self.Comment = ""
    self.ROIColor = "g"   #color that the ROIs ill be plotted, value set = "rgbcmykwrg"
    self.ROIs = []    #List of ROIs, can change as ROIs are moved to fit image"
    self.initialROIs=[] #fixed list of initial ROIs that cannot change
    self.showBackgroundROI=False    #flag to show background ROI to determine if an ROI should be discarded
    self.showSNRROI=False           #flag to show a noise ROI to determine noise from an image subtraction
    dROI=ROI()
    if name == "":  #If no name is given create one default ROI
        self.ROIs.append(dROI)
        self.nROIs = 1
    else:
        self.nROIs = 0
        
  def SetROIsParameter(self, n, param, value):
    nROI = self.ROIs[n]
    nROI.SetROIParameter(param,value)

  def translate(self,vector):
      for roi in self.ROIs:
          roi.translate(vector)
          
  def rotate(self,axis,theta):
      for roi in self.ROIs:
          roi.rotate(axis,theta)
          
  def printROIs (self):   #for debugging 
      print ("ROI name=" + self.ROIName +  " ,ROI number=" + str(self.nROIs)) 
      for roi in self.ROIs:
        print ("ROI= " + str(roi.Index) + ", ROI position= " + str(roi.Xcenter) + ", " + str(roi.Ycenter) + ", " + str(roi.Zcenter)) 
  

      
class ROI():
  ''' Class to define regions of interest (ROI)
      Each ROI has a position, geometrical shape, size, material properties 
  '''
  def __init__(self):
    self.Name = "default"
    self.Index = 0
    self.Type = "Sphere"  #Sphere or Rectangle, 
    self.Label = 'index'   #label to be printed next to ROI
    self.Xcenter = 0.0  #position in mm
    self.Ycenter = 0.0
    self.Zcenter = 0.0
    self.dx=1.0   #length and width of rectangle
    self.dy=1.0
    self.theta=0    #angle of rectangle in degrees
    self.d1 = 10.0  #diameter in mm
    self.T1 = 200.  #ROI prescribed T1 in ms, will be a function of temperature and field
    self.T2 = 100.  #ROI prescribed T2 in ms
    self.PD = 100   #prescribed proton density in percent
    self.Concentration = 100
    self.ADC = 1.1    #Apparent diffusion constant in 10^-3 mm^2/s
    self.SignalAve = 0.  #average value within ROI
    self.SignalRMS = 0.
    self.array=np.zeros((10,10))        #array of voxel values
    self.color='k'    #ROI color to differentiate ROIs when data sets are plotted, color set "rgbcmykwrg" or (r,g,b)
    self.symbol='o'   #symbol to represent ROI in data plots
    
  def SetROIParameter(self, param, value):
    if param == "Index":
            self.Index = value
    if param == "Type":
            self.Type = value
    if param == "Xcenter":
            self.Xcenter = value
    if param == "Ycenter":
            self.Ycenter = value
    if param == "Zcenter":
            self.Zcenter = value
    if param == "d1":
            self.d1 = value
    if param == "T1":
            self.T1 = value
    if param == "T2":
            self.T2 = value 
    if param == "PD":
            self.PD = value   
    if param == "ADC":
            self.ADC = value                     
    if param == "Concentration":
            self.Concentration = value 
 
  def Rcenter(self):
      return np.array([self.Xcenter , self.Ycenter, self.Zcenter])
    
  def translate(self, vector):
      self.Xcenter += vector[0]
      self.Ycenter += vector[1]
      self.Zcenter += vector[2]
      
  def rotate(self, axis, dtheta):
      '''rotate ROI position and orientation by theta in radians'''
      r= rotateVector(self.Rcenter(), axis, dtheta)
      self.Xcenter = r[0]
      self.Ycenter = r[1]
      self.Zcenter = r[2]
      self.theta+=dtheta*180.0/np.pi     

  def ROIUpperLeft(self):
      pass   
    
def rotateVector( v, axis, theta):
    """
    rotates vector v by  rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians. Euler-Rodriques formula   to rotate v np.dot(rotation_matrix(axis,theta), v)
    """
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis/np.sqrt(np.dot(axis, axis))
    a = np.cos(theta/2)
    b, c, d = -axis*np.sin(theta/2)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    R = np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]]) 
    return     np.dot(R, v)