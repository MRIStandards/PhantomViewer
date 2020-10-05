

import os
import string
from __main__ import vtk, qt, ctk, slicer
from DICOMLib import DICOMPlugin
from DICOMLib import DICOMLoadable


#
# This is the plugin to handle translation of DICOM objects
# that can be represented as multivolume objects
# from DICOM files into MRML nodes.  It follows the DICOM module's
# plugin architecture.
#


class DICOMPhilipsRescalePluginClass(DICOMPlugin):
  """ MV specific interpretation code
  """


  def __init__(self,epsilon=0.01):
    super(DICOMPhilipsRescalePluginClass,self).__init__()
    self.loadType = "PhilipsRescaledVolume"


    self.tags['seriesInstanceUID'] = "0020,000E"
    self.tags['seriesDescription'] = "0008,103E"
    self.tags['position'] = "0020,0032"
    self.tags['studyDescription'] = "0008,1030"


    # tags used to rescale the values
    self.philipsVolumeTags = {}
    self.philipsVolumeTags['PrivateScaleIntercept'] = "2005,100d"
    self.philipsVolumeTags['PrivateScaleSlope'] = "2005,100e"
    self.philipsVolumeTags['PrivateImageType'] = "2005,1011"
    self.philipsVolumeTags['ScaleIntercept'] = "0028,1052"
    self.philipsVolumeTags['ScaleSlope'] = "0028,1053"
    self.philipsVolumeTags['Manufacturer'] = "0008,0070"


    for tagName,tagVal in self.philipsVolumeTags.iteritems():
      self.tags[tagName] = tagVal


    self.epsilon = epsilon


  def examine(self,fileLists):
    """ Returns a list of DICOMLoadable instances
    corresponding to ways of interpreting the 
    fileLists parameter.
    """
    loadables = []
    allfiles = []
    scalarVolumePlugin = slicer.modules.dicomPlugins['DICOMScalarVolumePlugin']()
    for files in fileLists:
      fileSuitable = True
      for f in files:
        manufacturer = slicer.dicomDatabase.fileValue(f, self.tags['Manufacturer'])
        imageType = slicer.dicomDatabase.fileValue(f,self.tags['PrivateImageType'])
        if string.find(manufacturer, 'Philips') == -1 or imageType != "M":
          fileSuitable = False
          break
      if not fileSuitable:
        continue
      loadables += scalarVolumePlugin.examine([files])
      loadables[-1].name = loadables[-1].name+' - rescaled based on Philips private tags'
      loadables[-1].warning = 'Warning: Do not use for reading derived maps such as ADC, FA etc.'
      loadables[-1].confidence = 0.15 # make sure standard scalar volume reader takes preference


    return loadables


  def load(self,loadable):
    """Load the selection as a scalar volume, but rescale the values
    """


    scalarVolumePlugin = slicer.modules.dicomPlugins['DICOMScalarVolumePlugin']()
    vNode = scalarVolumePlugin.loadFilesWithArchetype(loadable.files, loadable.name)


    if vNode:
      intercept = float(slicer.dicomDatabase.fileValue(loadable.files[0], self.tags['ScaleIntercept']))
      slope = float(slicer.dicomDatabase.fileValue(loadable.files[0], self.tags['ScaleSlope']))
      privateIntercept = float(slicer.dicomDatabase.fileValue(loadable.files[0], self.tags['PrivateScaleIntercept']))
      privateSlope = float(slicer.dicomDatabase.fileValue(loadable.files[0], self.tags['PrivateScaleSlope']))




      print('Slope: '+str(slope)+' private: '+str(privateSlope)+' Intercept: '+str(intercept)+' private: '+str(privateIntercept))
      
      # vtkImageShiftScale first shifts, then scales


      # First, revert the scaling applied on ITK-level read
      reverseShift = vtk.vtkImageShiftScale()
      reverseShift.SetShift(-1.*intercept)
      reverseShift.SetScale(1)
      reverseShift.SetInput(vNode.GetImageData())
      reverseShift.SetOutputScalarTypeToFloat()


      reverseScale = vtk.vtkImageShiftScale()
      reverseScale.SetShift(0)
      reverseScale.SetScale(1./slope)
      reverseScale.SetInput(reverseShift.GetOutput())
      reverseScale.SetOutputScalarTypeToFloat()


      # Second, apply scaling using the private tags information
      rescale = vtk.vtkImageShiftScale()
      rescale.SetShift(-1.*privateIntercept)
      rescale.SetScale(1./privateSlope)
      rescale.SetOutputScalarTypeToFloat()
      rescale.SetInput(reverseScale.GetOutput())
      rescale.Update()


      imageData = vtk.vtkImageData()
      imageData.DeepCopy(rescale.GetOutput())


      # Note: the assumption here is that intercept/slope are identical for all
      # slices in the series. According to Tom Chenevert, this is typically the
      # case: "The exception is when there are multiple image types in a series,
      # such as real, imaginary, magnitude and phase images all stored in the
      # series.  But this is not common."
      vNode.SetAndObserveImageData(imageData)


    return vNode


#
# DICOMPhilipsRescalePlugin
#


class DICOMPhilipsRescalePlugin:
  """
  This class is the 'hook' for slicer to detect and recognize the plugin
  as a loadable scripted module
  """
  def __init__(self, parent):
    parent.title = "DICOM Philips Volume Rescale+Import Plugin"
    parent.categories = ["Developer Tools.DICOM Plugins"]
    parent.contributors = ["Andrey Fedorov, BWH; Tom Chenevert, Univeristy of Michigan"]
    parent.helpText = """
    This plugin addresses an issue with some images produced by Philips
    scanners, where the values stored in PixelData need to be rescaled using
    the information saved in the private tags to obtain quantitative
    measuremenes. The rescale formula is the following:
    QuantitativeValue = [SV-ScInt] / ScSlp,
    where SV = stored DICOM pixel value, ScInt = Scale Intercept =
    (2005,100d), ScSlp = Scale Slope = (2005,100e)
    This information was provided by Tom Chenevert, U.Michigan, as part of NCI
    Quantitative Imaging Network Bioinformatics Working Group activities.
    """
    parent.acknowledgementText = """
    This DICOM Plugin was developed by 
    Andrey Fedorov, BWH,
    and was partially funded by NIH grant U01CA151261.
    """


    # don't show this module - it only appears in the DICOM module
    parent.hidden = True


    # Add this extension to the DICOM module's list for discovery when the module
    # is created.  Since this module may be discovered before DICOM itself,
    # create the list if it doesn't already exist.
    try:
      slicer.modules.dicomPlugins
    except AttributeError:
      slicer.modules.dicomPlugins = {}
    slicer.modules.dicomPlugins['DICOMPhilipsRescalePlugin'] = DICOMPhilipsRescalePluginClass


#
#


class DICOMPhilipsRescalePluginWidget:
  def __init__(self, parent = None):
    self.parent = parent
    
  def setup(self):
    # don't display anything for this widget - it will be hidden anyway
    pass


  def enter(self):
    pass
    
  def exit(self):
    pass

 