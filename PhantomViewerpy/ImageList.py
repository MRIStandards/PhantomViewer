'''
Created on Jun 2, 2014
Class to import a variety of magnetic resonance image files and generate an image list
The image list contains a set of lists that contains header information, image parameters, and image data
ImageList has a default null image at list position 0, the first image is then at index 1.  
  This is so if you delete all image there still is the default image to be displayed. 
This class contains a method to write an ImageList to DICOM, not a good way to modify DICOM files since not all of the header information is preserved
To write out animated gifs requires visvis "pip install visvis" :Visvis is a pure Python library for visualization of 1D to 4D data in an object oriented way.
Last modification 6-3-14
@author: stephen russek
'''
import numpy as np
import datetime, time
import pydicom    #pydicom is used to import DICOM images  pydicom.UID
from pydicom.dataset import Dataset, FileDataset
from fdftopy import VarianData
try:
  import Image  #imports tif gif
except:
  pass     
try:
  import imageio
except:
   print ('Can not import imageio needed for animated GIFs')
import struct

class ImageList():
  '''Class to contain a generic MR Image list; can accommodate DICOM, tif, fdf (Varian) files; first element in ImageList is a default image '''
  def __init__(self ,parent = None):    #all attributes are list of the same length, **all attributes have default value at 0**
    #Image stack parameters
    self.ImageStackName='default' 
    #image stack arrays
    self.bValue = []
    self.bValue.append(0.0)
    self.ColumnDirection=[]
    self.ColumnDirection.append(np.array([0.,1.,0.]))  #list of unit vectors along columns, default columns along y axis 
    self.Columns=[]   #list of number of columns in each image (int)
    self.Columns.append(128)
    self.Comment = []
    self.Comment.append("") 
    self.DataType = []    #Data type string can be "real", "imag", "absval", "phase" or "complex"
    self.DataType.append("")       
    self.FA =[] #Flip angle in degrees
    self.FA.append(0.)  
    self.FileName = []
    self.FileName.append( "" )   #first element of list is null data to clear all images and labels
    self.FileType = []    #file type is  dcm  for DICOM; fdf for Varian flexible format; fid for Varian free induction decay; .tiff for tiff
    self.FileType.append("")    
    self.FoVX =[]       #Field of View in the horizontal direction in mm
    self.FoVX.append(128.)
    self.FoVY =[]       #Field of View in the vertical direction in mm
    self.FoVY.append(128.)
    self.header = []    #String containing full file header
    self.header.append("")
    self.InPlanePhaseEncodingDirection=[]
    self.InPlanePhaseEncodingDirection.append("")
    self.InstitutionName= []
    self.InstitutionName.append("") 
    self.ImagePosition= []      #upper left corner coordinates
    self.ImagePosition.append(np.array([0,0,0]))
    self.ImageCenter= []      #Image center coordinates
    self.ImageCenter.append(np.zeros(3))
    self.ImagingFrequency= []      #Image center coordinates
    self.ImagingFrequency.append(0.0)

    self.MagneticFieldStrength= []
    self.MagneticFieldStrength.append(0.0)    
    self.Manufacturer= []
    self.Manufacturer.append("")
    self.PatientName= []
    self.PatientName.append("")    
    self.ProtocolName= []
    self.ProtocolName.append("")        
    self.PixelBandwidth= []
    self.PixelBandwidth.append("")        
    self.PA =[] #pixel array list, list of numpy arrays
    self.PA.append(np.zeros([128,128]))
    self.PixelSpacingX =[]
    self.PixelSpacingX.append(1.)
    self.PixelSpacingY =[]
    self.PixelSpacingY.append(1.)    
    self.ReceiveCoilName= []
    self.ReceiveCoilName.append("")    
    self.Rows=[]    #list of number of rows in each image (int)
    self.Rows.append(128)
    self.RowDirection=[]
    self.RowDirection.append(np.array([1.,0.,0.]))    #list of unit vectors along rows, default rows along x axis 
    self.ScaleSlope = []    #Phillip's scaling parameters
    self.ScaleSlope.append (1.0)  #default value is 1
    self.ScaleIntercept = []
    self.ScaleIntercept.append (0.0)    #default value is 0
    self.SeriesDescription= []
    self.SeriesDescription.append("") 
    self.SliceThickness =[]
    self.SliceThickness.append(0.)
    self.SliceLocation =[]
    self.SliceLocation.append(np.nan)
    self.StudyDate = []
    self.StudyDate.append("")
    self.TR =[] #Repetition time in ms
    self.TR.append(0.)
    self.TE =[] #Echo time in ms
    self.TE.append(0.)
    self.TI =[] #Inversion time in ms
    self.TI.append(0.)

  def deleteImage(self,i):  #deletes ith image from list
    for at in self.__dict__.keys():
      if len(getattr(self,at)) >= i+1 and isinstance(getattr(self,at), list):   #make sure the item is a list and it has an ith element to delete
        del getattr(self,at)[i]
    return 'Deleted image ' + str(i)

  def addImage(self,pa):  #adds image to end of list
    for at in self.__dict__.keys():
      slist=getattr(self,at)
      if  isinstance(slist, list):   #make sure the item is a list and it has an ith element to delete
        slist.append(slist[-1])
    self.PA[-1]=pa
    return 'Added image to stack'

  def sort_list(self, list1, list2, reverse=False):
    '''sorts list1 using lit2, returns sorted list1''' 
    zipped_pairs = zip(list2, list1) 
    z = [x for _, x in sorted(zipped_pairs, reverse=reverse)]
    return z 
  
  def sortImageList(self,rlist, reverse=False):  #reorders all lists according to rlist
    for at in self.__dict__.keys():
      slist=getattr(self,at)
      if  isinstance(slist, list):   #make sure the item is a list
        l=self.sort_list(slist,rlist, reverse=reverse)
        setattr(self, at,l)
    return 
          
  def np3dArray(self):  #needs work
    nim=len(self.PA)
    if self.PA[1].ndim==3:
      return self.PA[1]
    else:
      try: 
        np3d=np.stack(self.PA[1:nim], axis=0)
      except:
        np3d=0
      return np3d
            
  def addFile(self,fileName):  
    extension = str(fileName.split(".")[-1])
    try:
      if extension.lower() == "dcm" or extension==fileName or extension.lower() == "ima": # if .dcm , .ima or  no extension try dicom
          self.unpackImageFile (pydicom.read_file(str(fileName)), fileName, "dcm")
      elif extension.lower() == "tif":
          #self.ds.PA.append(Image.open(str(fileName)))
          self.unpackImageFile (Image.open(str(fileName)), fileName, "tif")
      elif extension.lower() == "fdf":
          VData=VarianData()
          self.unpackImageFile (VData.read(str(fileName)), fileName, "fdf")
      else:
          self.unpackImageFile (pydicom.read_file(str(fileName)), fileName, "dcm")    #if the extension cannot be recognized try DICOM
    except:
       return [False,"Error: Image file cannot be opened:" + fileName]
    return [True, fileName]
        
  def    writeDicomFiles(self, filename):         
        for i in range(1,len(self.PA)):     #Image 0 is a null image for default display and will not be written out
            fileName = filename + str(i) + ".dcm"
#         print ("printing dicom file " + fileName)         
            # Populate required values for file meta information
            file_meta = Dataset()
            file_meta.MediaStorageSOPClassUID = 'Secondary Capture Image Storage'
            file_meta.MediaStorageSOPInstanceUID = '1.3.6.1.4.1.9590.100.1.1.111165684411017669021768385720736873780'
            file_meta.ImplementationClassUID = '1.3.6.1.4.1.9590.100.1.0.100.4.0'
            # Create the FileDataset instance (initially no data elements, but file_meta supplied)
            ds = FileDataset(fileName, {}, file_meta=file_meta, preamble="\0"*128)
            ds.Modality = 'MR'
            ds.ContentDate = str(datetime.date.today()).replace('-','')
            ds.ContentTime = str(time.time()) #milliseconds since the epoch
            ds.StudyInstanceUID =    '1.3.6.1.4.1.9590.100.1.1.124313977412360175234271287472804872093'
            ds.SeriesInstanceUID = '1.3.6.1.4.1.9590.100.1.1.369231118011061003403421859172643143649'
            ds.SOPInstanceUID =        '1.3.6.1.4.1.9590.100.1.1.111165684411017669021768385720736873780'
            ds.SOPClassUID = 'MR Image Storage'
            ds.SecondaryCaptureDeviceManufctur = 'Python 2.7.3'
        ## These are the necessary imaging components of the FileDataset object.
            ds.SamplesPerPixel = 1
            ds.PhotometricInterpretation = "MONOCHROME2"
            ds.PixelRepresentation = 0
            ds.HighBit = 15
            ds.BitsStored = 16
            ds.BitsAllocated = 16
            ds.SmallestImagePixelValue = '\\x00\\x00'
            ds.LargestImagePixelValue = '\\xff\\xff'
            #Add MRI data elements 
            ds.bValue= self.bValue[i]
            ds.Columns=self.Columns[i]
            ds.FoVX=self.FoVX[i]
            ds.FoVY=self.FoVY[i]
            ds.ImageOrientationPatient=self.ImagePosition[i].tolist()
            ds.InstitutionName=self.InstitutionName[i]
            ds.PixelBandwidth= self.PixelBandwidth[i] 
            ds.PixelSpacing =[self.PixelSpacingY[i],self.PixelSpacingX[i]]         
            ds.Rows= self.Rows[i]
            ds.PatientName = self.PatientName[i]
            ds.PatientID = "123456"
            ds.ProtocolName=self.ProtocolName[i]
            ds.RepetitionTime = self.TR[i]
            ds.SeriesDescription=self.SeriesDescription[i]
            ds.EchoTime = self.TE[i]
            ds.FlipAngle= self.FA[i]
            ds.InversionTime=self.TI[i] 
            ds.ImageOrientationPatient[3:6]=self.ColumnDirection[i].tolist() 
            ds.ImageOrientationPatient[:3]=self.RowDirection[i].tolist()
            ds.MagneticFieldStrength = self.MagneticFieldStrength[i]
            ds.Manufacturer = self.Manufacturer[i]            
            pixel_array=np.transpose(self.PA[i])        #numpy    pixel array
            if self.DataType[i].lower() == "phase" :
                pixel_array=(pixel_array+np.pi) *10000    # phase data -pi to pi is converted to 0 to 10000pi
                pixel_array = pixel_array.astype(np.uint16)
                #print "Adjusting phase to 16 bit integer"
            if pixel_array.dtype != np.uint16:
                scale=65536/np.amax(pixel_array)
                pixel_array*=scale       #rescale so that max value is 12 bits                             
                pixel_array = pixel_array.astype(np.uint16)
            ds.PixelData = pixel_array.tostring()    # image byte data
            # Set the transfer syntax
            ds.is_little_endian = True
            ds.is_implicit_VR = True
            ds.save_as(fileName)

  def writeAnimatedGIF(self, filename, duration=0.2):         
        filename = filename  + ".gif"
        images=[]
        for im in self.PA[1:]:
           im=np.transpose(im)
           images.append(im)
        imagesa=np.asarray(images)   
        imagesa-=imagesa.min()    #set lowestvalue to 0
        im=imagesa.astype(float)
        im/=(im.max()/255.0)  #set highest value to 255
        imageio.mimsave(filename, im.astype(np.uint8), format='GIF', duration=duration)        

        
  def unpackImageFile(self, ImageFile, FileName, fType): 
    """Unpacks several types of image files and appends all attributes into image stack  lists"""
    self.FileName.append(FileName)
    if fType == "dcm": 
      self.header.append(str(ImageFile))
    elif fType == "fdf":
      self.header.append(ImageFile.header)
    elif fType == "tif":
      self.header.append("tif")
    else:
      self.header.append("file not recognized")
    b = 0.0
    if hasattr(ImageFile,"bValue"):   #b value can be in several tags
      b=ImageFile.bValue 
    try:    #Siemens Trio   
      b=float(ImageFile[0x0019, 0x100C].value)
    except:
      pass   
    if hasattr(ImageFile,"Diffusionbvalue"): #Phillips tag (0018, 9087
      b=ImageFile.Diffusionbvalue
    try:    #GE b values often at (0043, 1039) in first of 5 element string array   
      b=float(ImageFile[0x0043, 0x1039].value[0])
    except:
      pass 
    self.bValue.append(b)
    self.Columns.append(ImageFile.Columns) if hasattr(ImageFile,"Columns") else self.Columns.append(0) 
#ImageOrientationPatient is a list of strings
    self.ColumnDirection.append(np.asfarray(ImageFile.ImageOrientationPatient[3:6])) if hasattr(ImageFile,"ImageOrientationPatient") else self.ColumnDirection.append(np.array([0.,1.,0.]))
    self.Comment.append(ImageFile.Comment) if hasattr(ImageFile,"Comment") else self.Comment.append("") 
    self.DataType.append(ImageFile.DataType) if hasattr(ImageFile,"DataType") else self.DataType.append("")
    self.Rows.append(ImageFile.Rows) if hasattr(ImageFile,"Rows") else self.Rows.append(0)
    self.StudyDate.append(ImageFile.StudyDate) if hasattr(ImageFile,"StudyDate") else self.StudyDate.append("")
    self.Manufacturer.append(ImageFile.Manufacturer) if hasattr(ImageFile,"Manufacturer") else self.Manufacturer.append("")
    self.SeriesDescription.append(ImageFile.SeriesDescription) if hasattr(ImageFile,"SeriesDescription") else self.SeriesDescription.append("")
    self.InstitutionName.append(ImageFile.InstitutionName) if hasattr(ImageFile,"InstitutionName") else self.InstitutionName.append("")     
    self.MagneticFieldStrength.append(ImageFile.MagneticFieldStrength) if hasattr(ImageFile,"MagneticFieldStrength") else self.MagneticFieldStrength.append(0.0) 
    self.ImagingFrequency.append(ImageFile.ImagingFrequency) if hasattr(ImageFile,"ImagingFrequency") else self.ImagingFrequency.append(0.0)         
    self.PatientName.append(ImageFile.PatientName) if hasattr(ImageFile,"PatientName") else self.PatientName.append("")    
    self.PixelBandwidth.append(ImageFile.PixelBandwidth) if hasattr(ImageFile,"PixelBandwidth") else self.PixelBandwidth.append(0.0)        
    self.PixelSpacingX.append(ImageFile.PixelSpacing[1]) if hasattr(ImageFile,"PixelSpacing") else self.PixelSpacingX.append(1.0)   #column spacing, distance between neighboring columns
    self.PixelSpacingY.append(ImageFile.PixelSpacing[0]) if hasattr(ImageFile,"PixelSpacing") else self.PixelSpacingY.append(1.0)  #row spacing, distance between neighboring rows  
    self.ProtocolName.append(ImageFile.ProtocolName) if hasattr(ImageFile,"ProtocolName") else self.ProtocolName.append("")    
    self.ImagePosition.append(np.asfarray(ImageFile.ImagePositionPatient)) if hasattr(ImageFile,"ImagePositionPatient") else self.ImagePosition.append(np.array([1.,0.,0.]))
    self.ImagePosition.append(ImageFile.ImagePositionPatient) if hasattr(ImageFile,"ImagePositionPatient") else self.ImagePosition.append(np.array([1.,0.,0.]))
    self.ReceiveCoilName.append(ImageFile.ReceiveCoilName) if hasattr(ImageFile,"ReceiveCoilName") else self.ReceiveCoilName.append("")
    self.RowDirection.append(np.asfarray(ImageFile.ImageOrientationPatient[:3])) if hasattr(ImageFile,"ImageOrientationPatient") else self.RowDirection.append(np.array([1.,0.,0.]))
    self.TR.append(ImageFile.RepetitionTime) if hasattr(ImageFile,"RepetitionTime") else self.TR.append(0.0) 
    self.TE.append(ImageFile.EchoTime) if hasattr(ImageFile,"EchoTime") else self.TE.append(0.0)
    self.FA.append(ImageFile.FlipAngle) if hasattr(ImageFile,"FlipAngle") else self.FA.append(0.0)       
    self.InPlanePhaseEncodingDirection.append(ImageFile.InPlanePhaseEncodingDirection) if hasattr(ImageFile,"InPlanePhaseEncodingDirection") else self.InPlanePhaseEncodingDirection.append('none')    
    self.TI.append(ImageFile.InversionTime) if hasattr(ImageFile,"InversionTime") else self.TI.append(0.0) 
    self.SliceThickness.append(ImageFile.SliceThickness) if hasattr(ImageFile,"SliceThickness") else self.SliceThickness.append(0.0)    
    self.SliceLocation.append(ImageFile.SliceLocation) if hasattr(ImageFile,"SliceLocation") else self.SliceLocation.append(0.0)

# Phillips scaling corrections FP = (PV-SI)/SS = PV/SS where FP is the floating point value, PV is the pixel value, SI and SS are the scaled slope and intercept
# SS is in [0x2005, 0x100E] and si is in [0x2005, 0x100D], both are single precision floats
    try:
      ss=struct.unpack('f',ImageFile[0x2005, 0x100E].value)[0]
    except:
      ss=1.0
    try:
      si=struct.unpack('f',ImageFile[0x2005, 0x100D].value)[0]
    except:
      si=0.0
    self.ScaleIntercept.append(si)  
    self.ScaleSlope.append(ss) 
# Derived quantities
    self.FoVX.append(ImageFile.FoVX) if hasattr(ImageFile,"FoVX")  else self.FoVX.append(self.PixelSpacingX[-1]*self.Columns[-1])
    self.FoVY.append(ImageFile.FoVY) if hasattr(ImageFile,"FoVY")  else self.FoVY.append(self.PixelSpacingY[-1]*self.Rows[-1]) 
#    self.ImageCenter.append(self.ImagePosition[-1]+self.FoVX[-1]/2*self.RowDirection[-1]+self.FoVY[-1]/2*self.ColumnDirection[-1])
#add pixel arrays
    if fType == "dcm": #if data set is a dicom file display its pixel array
        self.PA.append(np.transpose(ImageFile.pixel_array.astype(float)))  #Convert to float arrays, ???why transpose????
        
    elif fType == "tif":
        self.PA.append(np.array(ImageFile))
    elif fType == "fdf":
      if len(ImageFile.matrix) == 2:  #Varian file has 2d array    
        self.PA.append(ImageFile.PA)
      if len(ImageFile.matrix) == 3:  #Varian file has 3d array
        nim = int(ImageFile.matrix[2])
        self.PA.append(ImageFile.PA[:,:,0])
        for i in range(1,nim):
          self.addImage(ImageFile.PA[:,:,i])  #adds image to end of list
    else:
        self.PA.append(np.zeros([128,128]))    
    
