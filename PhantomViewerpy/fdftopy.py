# -*- coding: utf-8 -*-
"""
Created on Wed Dec 04 12:19:26 2013

@author: stephen russek
"""

import os
import re
import numpy as np
import struct

class VarianData:
  ''' Unpacks Varian fdf files'''
  def __init__(self):
    self.bValue = 0.0
    self.Columns = 0
    self.ColumnDirection= [0,0,0]
    self.DataType = ""
    self.EchoTime = 1.0
    self.FileType= "fdf"
    self.FlipAngle = 0.0   
    self.FoVX = 50.
    self.FoVY = 50.   
    self.ImageOrientationPatient = []
    self.InversionTime= 0.0
    self.fmt = ""
    self.header= ""
    self.Manufacturer = "Agilent"
    self.matrix = []
    self.PA=np.zeros([128,128])  #pixel array
    self.PixelSpacing = [1.0,1.0]
    self.ProtocolName = ""   
    self.RepetitionTime = 1.0   
    self.Rows = 0
    self.RowDirection= [1,0,0]
    self.StudyDate = ""
    self.SeriesDescription = ""
    self.SliceLocation = ""   #given bugs in fdf output slicelocation is set to slice number
    self.ro = 0      #number of readout points
    self.pe = 0      #number of phase encode points


  def read( self, filename ):
    if filename.endswith(".fdf"):
      data = self.readFDF( filename )
    elif filename.endswith(".img"):
      data = self.readIMG( filename )
    else:
      print ("Unknown filename %s " % (filename))
    return data

  def readFDF(self, filename ):
    fdfImage = VarianData()
    fp = open( filename, "rb" )

    xsize = -1
    ysize = -1
    zsize = 1
    bigendian = -1
    done = False

    while not done :
      line = fp.readline()  #read in as bytes
      line=line.decode('utf-8')   #modified to work in Python3
      fdfImage.header += line 

      if( len( line ) >= 1 and line[0] == chr(12) ):
        break
    
      if( len( line ) >= 1 and line[0] != chr(12) ):    #unpack header    
          if( line.find("bigendian") > 0 ):
             endian = line.split("=")[-1].rstrip("\n; ").strip(" ")       
          if( line.find("bvalue") > 0 ):
             fdfImage.bValue = float(line.split("=")[-1].rstrip("\n; ").strip(" "))    
          if( line.find("*type") > 0 ):
             fdfImage.DataType = line.split("=")[-1].rstrip("\n; ").strip(" ").replace('"', '')            
          if( line.find("echos") > 0 ):
             nechoes = line.split("=")[-1].rstrip("\n; ").strip(" ")
          if( line.find("orientation") > 0 ):
             orient= line.split("=")[-1].rstrip("\n; ").strip(" ")
             fdfImage.ImageOrientationPatient = orient.replace("{"," ").replace("}"," ").split(",")
             #print fdfImage.ImageOrientationPatient
          if( line.find("studyid") > 0 ):
             fdfImage.SeriesDescription = line.split("=")[-1].rstrip("\n; ").strip(" ").replace('"', '')
          if( line.find("sequence") > 0 ):
             fdfImage.ProtocolName = line.split("=")[-1].rstrip("\n; ").strip(" ").replace('"', '')
          if( line.find("span") > 0 ):
             span = line.split("=")[-1].rstrip("\n; ").strip(" ")
             span=span.replace("{"," ").replace("}"," ")
             fdfImage.FoVX = float(span.split(",")[0]) * 10.        #asuming cm and converts to mm, needs work
             fdfImage.FoVY = float(span.split(",")[1]) * 10.
          if( line.find("TR =") > 0 ):
             fdfImage.RepetitionTime = float(line.split("=")[-1].rstrip("\n; ").strip(" "))             
          if( line.find("TE =") > 0 ):
             fdfImage.EchoTime = float(line.split("=")[-1].rstrip("\n; ").strip(" "))             
          if( line.find("TI =") > 0 ):
             fdfImage.InversionTime = float(line.split("=")[-1].rstrip("\n; ").strip(" ")) 
          if( line.find("ro_size") > 0 ):
             fdfImage.ro = line.split("=")[-1].rstrip("\n; ").strip(" ")
             #print "ro found=" + str(fdfImage.ro)             
          if( line.find("pe_size") > 0 ):
             fdfImage.pe = line.split("=")[-1].rstrip("\n; ").strip(" ")
             #print "pe found=" + str(fdfImage.pe)
          if( line.find("echo_no") > 0 ):
             echo_no = line.split("=")[-1].rstrip("\n; ").strip(" ")
          if( line.find("nslices") > 0 ):
             nslices = line.split("=")[-1].rstrip("\n; ").strip(" ")        
          if( line.find("slice_no") > 0 ):
             sl = line.split("=")[-1].rstrip("\n; ").strip(" ")
          if( line.find("location") > 0 ):
             location = line.split("=")[-1].rstrip("\n; ").strip(" ")  
             fdfImage.SliceLocation = float(location[1:-2].split(',')[2]) *10   #last element in location string is slice location in cm      
          if( line.find("matrix") > 0 ):
             fdfImage.matrix = re.findall("(\d+)", line.rstrip())        
             if len(fdfImage.matrix) == 2:
               xsize, ysize = int(fdfImage.matrix[0]), int(fdfImage.matrix[1])
               #print "xsize,ysize=" +str(xsize) + str(ysize)
             elif len(fdfImage.matrix) == 3:
               xsize, ysize, zsize = int(fdfImage.matrix[0]), int(fdfImage.matrix[1]), int(fdfImage.matrix[2])
               #print "xsize,ysize,zsize=" +str(xsize) + str(ysize)+ str(zsize)
    fp.seek(-xsize*ysize*zsize*4,2) #set files current position xsize*ysize*zsize*4bytes from end of file
    
    if bigendian == 1:
        fdfImage.fmt = ">%df" % (xsize*ysize*zsize)
    else:
        fdfImage.fmt = "<%df" % (xsize*ysize*zsize)    #our images ar bigendian=0
    
    #print "fmt=" + str(fdfImage.fmt)
    fdfImage.Rows=float(fdfImage.pe)       #phase encode is normally along row
    fdfImage.Columns=float(fdfImage.ro)    #read out is normally along column
    fdfImage.PixelSpacing = [fdfImage.FoVY/fdfImage.Columns,fdfImage.FoVX/fdfImage.Rows]   #conform to DICOM standard
    data = struct.unpack(fdfImage.fmt, fp.read(xsize*ysize*zsize*4))
#    data = np.array( data ).reshape( [xsize, ysize, zsize ] ).squeeze()
    if len(fdfImage.matrix) == 2:
      fdfImage.PA = np.transpose(np.resize(data, [ysize,xsize]))
    if len(fdfImage.matrix) == 3:
      fdfImage.PA = np.transpose(np.resize(data, [ysize,xsize,zsize]))      
    #print "PA shape=" + str(fdfImage.PA.shape)
    fp.close()
    return fdfImage

  def readIMG(self, directory):
    
    # Get a list of all the FDF files in the directory
     try:
       files = os.listdir(directory)
     except:
       print ("Could not find the directory %s" % directory)
     return
    
     files = [ file for file in files if file.endswith('.fdf') ]
    
     data = []
     for file in files:
       data.append( self.readFDF( directory+"/"+file ) )
    
     data = transpose( array( data ), (1,2,0) ) 
    
     return data
