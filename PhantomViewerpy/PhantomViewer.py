# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 16:30:54 2013
PhantomViewer: Main class to display and analyze phantom data, mostly MRI data
modules required (most can be installed using 'pip install package' or  'pip install package --upgrade':
verify module version using "pip show module"
   Anaconda 3 with Python 3.7    (64 bit, Python, PyQT5, Numpy, Scipy, Pillow )
   pydicom Version: 1.4.2:       (for DICOM file import/export)
   pyqtgraph Version: 0.10.0     (for plotting and image windows)
   lmfit Version: 1.0.0          (for nonlinear least squares fitting routines)
   unwrap Version: 0.1.1
   scikit-image  Version: 0.11.3: uses unwrap_phase from skimage.resoration
   PyOpenGL                      (Used for 3D plotting)

will work with Anaconda 2 with Python 2.7.11, PYQT4 (not recommended)    
Uses PhantomViewerGui created from PhantomViewer.ui by Qt Designer
  convert ui file to python by executing:   
      "designer\pyuic4 designer\PhantomViewerGui.ui -o PhantomViewerpy\PhantomViewerGui4.py" 
  or  "designer\pyuic5 designer\PhantomViewerGui.ui -o PhantomViewerpy\PhantomViewerGui5.py" 
  from system shell to regenerate PhantomViewerGui.py from PhantomViewerGui.ui
@author: Stephen Russek
Units: times in ms, distances in mm, ADC in mm2/s, Temperature in C,

VPhantom: class to describe phantoms, several virtual phantoms, including the NIST/ISMRM system phantom are available
"""

import sys
import os             #operating system file/directory names
import copy           #copies objects eg ROIsets
import time           # use for time, date and sleep
import numpy as np
import numpy.ma as ma   #masked arrays eg ROI arrays
import threading      # to run tasks such as generating maps in background     
from pyqt import *         #imports required PyQt modules, tries PyQT4 then PyQt5
import pydicom        #import pydicom to read DICOM data sets and DICOMDIR
from pydicom.filereader import read_dicomdir
import pyqtgraph as pg    #uses pyqtgraph PlotWidget for graphs and ImageView windows for images, http://www.pyqtgraph.org/documentation/introduction.html
import pyqtgraph.opengl as gl
import pyqtgraph.functions as fn
import lmfit  #used for nonlinear least square fits,builds on Levenberg-Marquardt algorithm of scipy.optimize.leastsq(),
from astropy.nddata.utils import block_reduce #used for down sampling images
#from boto.sdb.db.sequence import double
#from mpl_toolkits.axes_grid1.anchored_artists import arr
from PIL import Image           #used for resizing 2d images
from scipy.interpolate import splrep, splev # Bspline representation of a curve and Bspline evaluation
from scipy.ndimage import zoom   #used for interpolation of images
from scipy import signal
try:
    from unwrap import unwrap   #package for phase unwrapping https://pypi.python.org/pypi/unwrap 
except:
  pass
from skimage.restoration import unwrap_phase

#PhantomViewer modules
if pyqtVersion==4:
  from PhantomViewerGui4 import Ui_PhantomViewerGui    #main window Gui
if pyqtVersion==5:
  from PhantomViewerGui5 import Ui_PhantomViewerGui    #main window Gui
import ROIProperties, ROIInfo, PhantomProperties
import T1IRabs, T1VFA, T1VTR, T2SE, DifModel, LCPowerLawFit # models for data fitting
import T1IRmap,  T2SEmap, DifMap     #modules to make parameter maps
import Resolution, Fiducial   #modules to perform automated resolution and geometric distortion measurements
import Info
import FitPlots   #Plot window to compare data to fits
import ImageList  #class to make an image list from a stack of image files
import DICOMDIRlist
import VPhantom, SystemPhantom, DiffusionPhantom,NISThcpPhantom, NISTKTPhantom, NISThcpCoronalPhantom

class PhantomViewer(QMainWindow):
  """Main ROI viewing window, customized for different analysis eg geometric distortion, T1, T2, SNR etc"""
  def __init__(self, dt,parent = None):
    super(PhantomViewer, self).__init__()
    pg.setConfigOption('background', 0.2)   #Background on plots 0 = black, 1 = white
    pg.setConfigOption('foreground', 'w')
    self.setAttribute(Qt.WA_NativeWindow, True)
    self.ui = Ui_PhantomViewerGui()
    self.ui.setupUi(self)
    self.modDate = '4/16/2020'
    self.wTitle = "PhantomViewer "  + self.modDate + ': ' 
    self.setWindowTitle(self.wTitle)
    self.nImages=0
    self.nCurrentImage=0
    self.ims=imageStackWindow(self)   #define separate window object  for the image stack
    self.imswin=self.ims.win      #image stack window
    self.imv = self.ims.imv    #image stack pyqtgraph image view object

    self.imv.ui.roiBtn.setText("Line scan/ROI")
#    self.imv.ui.histogram.plot.setLogMode(None,True)    #Not working in pyqtgraph 0.11!!set the histogram y axis to a log scale  
#    self.imv.ui.histogram.plot.setLogMode(None,True)    #set the histogram y axis to a log scale      
    self.imv.vLine = pg.InfiniteLine(angle=90, movable=False)   #cross hair
    self.imv.hLine = pg.InfiniteLine(angle=0, movable=False)
    self.imv.addItem(self.imv.vLine, ignoreBounds=True)
    self.imv.addItem(self.imv.hLine, ignoreBounds=True)
    self.proxy = pg.SignalProxy(self.imv.view.scene().sigMouseMoved, rateLimit=60, slot=self.mouseMoved)
    self.proxy2 = pg.SignalProxy(self.imv.view.scene().sigMouseClicked, rateLimit=60, slot=self.mouseClicked)
    self.rdPlot=self.ui.gvRawData   #plot widget for raw data
    self.resultsPlot=self.ui.gvResults  #plot widget for results of data fitting
    self.dicomHeader = "DICOM Header"
    self.ui.lblnImages.setText((str(self.nImages)))
    self.ui.lblCurrentImage.setText("none")
    self.ds= ImageList.ImageList()                         #list of data sets, can be dicom, tiff, fdf
    self.stackList=[]       #list of image stacks, first element is the current image stack
    self.seriesFileNames = []
    self.image3D = np.array ([1,1]) 
#signals and slots
# Menu items
#Images
    self.ui.actionSelectImages.triggered.connect(self.openFile)
    self.ui.actionSelectImages.setToolTip('Opens DICOM, DicomDir, tif, fdf ')  
    self.ui.actionClear_All_Images.triggered.connect(self.clearImages)
    self.ui.actionDelete_Current_Image.setShortcut('Ctrl+D')
    self.ui.actionDelete_Current_Image.triggered.connect(self.deleteCurrentImage)
#ImageStacks
    self.ui.actionSort_on_Slice_Loc.triggered.connect(self.sortOnSliceLoc)
    self.ui.actionSort_on_Slice_Loc_Reverse.triggered.connect(self.sortOnSliceLocReverse)
    self.ui.actionSort_on_TE.triggered.connect(self.sortOnTE)
    self.ui.actionShow_Image_Stack.triggered.connect(self.showImswin)
    self.ui.actionDuplicate_Image_Stack.triggered.connect(self.duplicateImageStack)
    self.ui.actionChange_Image_Stack.triggered.connect(self.changeImageStack)
    self.ui.actionWrite_Animated_GIF.triggered.connect(self.writeAnimatedGIF)
    self.ui.actionSave_ROI_arrays.triggered.connect(self.saveROIarrays)
#Phantoms    
    self.ui.actionOpenPhantomFile.triggered.connect(self.openPhantomFile)
    self.ui.actionSystem_Phantom.triggered.connect(self.SystemPhantom)
    self.ui.actionDiffusion_Phantom.triggered.connect(self.DiffusionPhantom)
    self.ui.actionBreast_Phantom.triggered.connect(self.BreastPhantom)
    self.ui.actionNIST_hcp_Phantom.triggered.connect(self.NISThcpPhantom)
    self.ui.actionNIST_hcp_Coronal_Phantom.triggered.connect(self.NISThcpCoronalPhantom)
    self.ui.actionNIST_KT_Phantom.triggered.connect(self.NISTKTPhantom)
    self.ui.actionShowPhantomInfo.triggered.connect(self.showPhantomProperties)
    self.ui.actionSave_phantom_file.triggered.connect(self.savePhantomFile) 
#ROIs       
    self.ui.actionOpen_ROI_File.triggered.connect(self.openROIFile)    
    self.ui.actionShow_ROI_Properties.triggered.connect(self.showROIInfo)
    self.ui.actionShowHide_ROIs.triggered.connect(self.toggleShowROIs)
    self.ui.actionReset_ROIs.triggered.connect(self.resetROIs)
    self.ui.actionThin.triggered.connect(self.ROIThin)
    self.ui.actionMedium.triggered.connect(self.ROIMedium) 
    self.ui.actionThick.triggered.connect(self.ROIThick) 
    self.ui.actionROI_Color.triggered.connect(self.ROIColor)
    self.ui.actionLabel_Color.triggered.connect(self.labelColor)   
#Analysis
    self.ui.actionT1_Analysis.triggered.connect(self.T1Analysis)
    self.ui.actionT2_Analysis.triggered.connect(self.T2Analysis)
    self.ui.actionProton_density_SNR.triggered.connect(self.PDSNRAnalysis)
    self.ui.actionFiducial.triggered.connect(self.fiducialAnalysis)
    self.ui.actionFiducial_Array_Analysis.triggered.connect(self.openFiducialAnalysisWindow)
    self.ui.actionResolution_Inset.triggered.connect(self.riAnalysis)
    self.ui.actionSlice_Profile.triggered.connect(self.sliceProfileAnalysis)
    self.ui.actionSlice_Profile_Setup.triggered.connect(self.sliceProfileSetup)
    self.ui.actionLC_Thermometer.triggered.connect(self.LCThermometerAnalysis)
    self.ui.actionDiffusion.triggered.connect(self.diffusionAnalysis)
    self.ui.actionPhase_Unwrap.triggered.connect(self.unWrapCurrentImage)
    self.ui.actionPlane_Background_Subtract.triggered.connect(self.planeBackgroundSubtract)
    self.ui.actionParabola_Background_Subtract.triggered.connect(self.ParabolaBackgroundSubtract)    
    self.ui.actionImage_Subtraction.triggered.connect(self.imageSubtraction)
    self.ui.actionDown_Sample.triggered.connect(self.downSample)
    self.ui.actionInterpolate.triggered.connect(self.interpolate)
    self.ui.actionCrop_to_ROI.triggered.connect(self.cropToROI)
    self.ui.actionPlot_ROI_vs_Index.triggered.connect(self.plotROIvsIndex)
    self.ui.actionSetScaleSlopeto1_0.triggered.connect(self.setSlopeOffset)
#Fit Options
    self.ui.actionMake_Map.triggered.connect(self.makeMap)
#Tools
    self.ui.action3dViewer.triggered.connect(self.View3d)
    self.ui.action3dViewer.setToolTip('Creates 3d rendering of image stack')
    self.ui.actionView3D_color.triggered.connect(self.view3DColor)
    self.ui.actionView3D_transparency.triggered.connect(self.view3DTransparency)
    self.ui.actionView3D_invert_image.triggered.connect(self.view3Dinvert)
    self.ui.actionWrite_to_DICOM.triggered.connect(self.writeDicomFiles)
    self.ui.actionFitting_Results.triggered.connect(self.viewReport)
    self.ui.actionSave_Results.triggered.connect(self.saveReport)
    self.ui.actionClear_Report.triggered.connect(self.clearReport)
    self.ui.actionSaveDataFiles.triggered.connect(self.saveData)
    self.ui.txtResults.setHidden(False)

#  push buttons and sliders
    self.ui.hsROISet.valueChanged.connect(self.changeROISet)
    self.ui.vsImage.valueChanged.connect(self.imageSlider)
    self.ui.pbReflectX.clicked.connect(self.reflectX)
    self.ui.pbReflectY.clicked.connect(self.reflectY)
    self.ui.pbReflectZ.clicked.connect(self.reflectZ)
    self.ui.pbRotate90X.clicked.connect(self.rotate90X)
    self.ui.pbRotate90Y.clicked.connect(self.rotate90Y)
    self.ui.pbRotate90Z.clicked.connect(self.rotate90Z)
    self.ui.pbShowRawData.clicked.connect(self.showRawData)
    self.ui.pbFitData.clicked.connect(self.fitData) 
    self.ui.pbViewFits.clicked.connect(self.viewFits)

    self.ui.pbSliceProfileSetup.clicked.connect(self.sliceProfileSetup) 
    self.ui.pbSliceProfileAnalize.clicked.connect(self.sliceProfileAnalysis)    
    
    self.ui.hsAngle.valueChanged.connect(self.rotateROIs) 
    self.ui.hsSize.valueChanged.connect(self.reSizeROIs)
    self.ui.rbEditROIset.clicked.connect(self.editROISet)
    self.ui.rbEditSingleROI.clicked.connect(self.editSingleROI)
    self.ui.rbAllSlices.clicked.connect(self.allSlices)
    self.ui.rbSelectedSlice.clicked.connect(self.currentSlice)
    self.ui.rbUseROIValues.clicked.connect(self.useROIValues)   #flag self.useROIValues = True  use nominal values as initial guess
    self.ui.rbUseBestGuess.clicked.connect(self.useBestGuess)    #flag self.useROIValues = False
    self.ui.rbResultsY.clicked.connect(self.replotResults)   #f 
    self.ui.rbResultsError.clicked.connect(self.replotResults)
    self.ui.rbViewDicomHeader.toggled.connect(self.viewDicomHeader)
    self.headerWindow=messageWindow(wtitle='Image Header')
    #self.ui.chShowBackgroundROIs.clicked.connect(self.showROIs)
    #self.ui.rbViewMessages.toggled.connect(self.viewMessages)
    self.ui.lblTE.textChanged.connect(self.TEChanged)
    self.TEChangedFlag = True   #flag to input new TE values
    self.ui.lblTR.textChanged.connect(self.TRChanged)
    self.TRChangedFlag = True   #flag to input new TR values
       
#  setup regions of interest (ROIs)
    self.dataType = str(dt)    #string indicating data type "T1", "T2", "PD-SNR", "Dif" ;"SliceProfile" determines what fitting models are accessible
    self.setDataType(self.dataType)
    self.dataHeader=''   #string to document data type and origin, used when outputting data to files
    self.ADCmap = False
    self.Phantom = VPhantom.VPhantom()
    self.ui.lblPhantomType.setText(self.Phantom.phantomName)
    self.ROIset = "NiArray"   #determines which ROI set to use via ROIsetdict 
    self.InitialROIs = VPhantom.ROISet("default")     #Current ROIs are initial ROIs with global rotation and translation
    self.currentROIs=copy.deepcopy(self.InitialROIs)
    self.currentROI = 1   #currently selected ROI
    self.useROIValues = False    #flag to instruct fits to take initial guesses from the current ROI values
    self.pgROIs=[]  #list of pyqtgraph ROIs
    self.pgROIlabels = []   #List of labels that go with ROIs
    self.bShowROIs = False    #Flag to determine if ROIs are shown or hidden
    self.bRoiLabel=True       #Flag to determine if ROI labels are shown
    self.roiPen = pg.mkPen('g', width=2) #one of: r, g, b, c, m, y, k, w    or R, G, B, [A]    integers 0-255  or (R, G, B, [A])  tuple of integers 0-255
    self.roiPenf = pg.mkPen('g', width=1)
    self.labelPen=pg.mkPen('y', width=1)
    self.lblColor=(255,204,0)
    self.lblFont=QFont()
    self.lblFont.setPixelSize(18)
    self.bShowSNRROIs = False    #flag to show SNR ROIs to measure background signal; used to determine points that are in the background
    self.bSNRROIs  = False      #flag to indicate SNR ROIs are plotted
    self.snrPen = pg.mkPen('c', width=3)
    self.ROItrans=np.array([0,0,0], float)  #ROI translation
    self.bResetROIs = True  #flag to reset ROIs to initials positions
    self.relativeOffseth = 0.0  #horizontal and vertical offsets of the ROI positions
    self.relativeOffsetv = 0.0
    self.theta = 0.0                 #current rotation of the ROIs in radians
    self.bEditROISet = True   #flag to specify whether to edit ROI set or an individual ROI
    self.bAllSlices = True   #flag to specify if all slices will be analyzed or just currently selected one
    self.view3DColor = QColor(255, 255 ,255 , alpha=10)
    self.view3DBackground = QColor(155, 155 ,255 , alpha=10)
    self.view3DTransparency = 1.75   #set transparency scaling for 3dview, 1 = transparency set by voxel value
    self.view3Dinvert = False    #flag to invert contrast in 3d image

#global data structures
    self.rdy = np.array([0.,0.]) # 2d numpy array of raw data (usually signal averaged i= ROI j=image number in stack)
    self.resy = np.array([0.,0.]) # 2d numpy array of residuals data-fit (usually signal averaged i= ROI j=image number in stack)
    self.rdx = np.array([0.]) # 1d numpy array of imaging parameter that is varied ie TE, TI, FA etc; usually a 1d array (independent variable)
    self.bdArray = np.array([0.]) # 2d numpy array of background around each ROI in each image
    self.rdBackground = np.array([0.])   #average background counts in signal free region in reduced image array
    self.background = 0.0    #background averaged over images
    self.ui.lblBackGround.setText(str(self.background))
    self.noisefactor= 1.0   #data below noisefactor*background will not be fit
    self.fity=np.zeros((14,100)) # 2d numpy array of fitting model using optimized parameters
    self.fitx = np.arange(100) # 100 element numpy array of imaging parameter spanning range of image data (independent variable)
    self.clickArray = []   # list of numpy array of points I, J, value generated from mouse clicks on the images
    self.clickArraymm = []   # Same as clickArray except in mm x,y, value 
    self.T1Params = []  #list of lmfit dictionaries
    self.report=""      # String to record all output to allow printing or saving a report
    self.imageDirectory = '' #last opened image directory
    self.showReferenceValues = True   #flag to show reference parameter values listed in phantom type or ROI file
    self.sliceList = []    #ordered list of slices in the current image stack
    self.reverseSliceOrder = False    #flag to reverse slice order
    self.reverseTEOrder = False    #flag to reverse TE order
    self.resultsPlotY = True        #flag to plot fit parameter in Results plot
    self.resultsPlotError = False        #flag to plot error relative to reference in Results plot
    self.messageLogOn=False         #flag to log messages (msgPrint) into messageLog
    self.messageLog=''
    self.msgPrint('Python='+ sys.version + '; pydicom=' + pydicom._version.__version__ + '; PyQt=' +PYQT_VERSION_STR +'\n')    #print Python version being used
    self.bgXSpacing=14.0   #distance between center of ROI and location where background is measured"
    self.bgYSpacing=12.5

#********************Setup phantom parameters************************************    
  def setupPhantom(self, phantom):
    self.Phantom=phantom
    self.ui.lblPhantomType.setText(self.Phantom.phantomName)
    self.currentROIs = self.Phantom.ROIsets[0]     #Current ROIs are initial ROIs with global rotation and translation
    self.ui.hsROISet.setMaximum(len(self.Phantom.ROIsets)-1)
    self.ui.hsROISet.setValue(0)
    self.ui.lblROISet.setText(self.Phantom.ROIsets[0].ROIName)
    self.ROIset = self.Phantom.ROIsets[0].ROIName   #determines which ROI set to use via ROIsetdict   
    self.InitialROIs=copy.deepcopy(self.currentROIs)    #make copy to reset if necessary
    self.roiPen = pg.mkPen(self.currentROIs.ROIColor, width=3)
    self.useROIValues = False   #default use best guess for initial fitting parameters, if true use ROI values
    self.ui.rbUseBestGuess.setChecked(True)   #does not seem to work
    self.showReferenceValues = True
    try:  #if user or image file sets a B0 field, use this to obtain correct reference values
        self.Phantom.changeB0(b0=float(self.ui.lblField.toPlainText()))
    except:
      pass
    self.theta=0.0
    self.ui.lblPhantomAngle.setText("0")
    self.ui.hsAngle.setValue(0)
    self.relativeOffseth = 0.0  #horizontal and vertical offsets of the ROI positions
    self.relativeOffsetv = 0.0
    self.redrawROIs()
        
#**************** Standard Phantoms *******************
  def SystemPhantom (self):
    self.setupPhantom(SystemPhantom.SystemPhantom())
   
  def DiffusionPhantom (self):
    self.setupPhantom(DiffusionPhantom.DiffusionPhantom())
    self.dataType = "Dif"
    self.setDataType(self.dataType)

  def BreastPhantom (self):
    pass

#*************NIST specific phantoms******************  
  def NISThcpPhantom (self):
    self.setupPhantom(NISThcpPhantom.hcpPhantom())
      
  def NISThcpCoronalPhantom (self):
    self.setupPhantom(NISThcpCoronalPhantom.hcpCoronalPhantom())
        
  def NISTKTPhantom (self):
    self.setupPhantom(NISTKTPhantom.KTPhantom())
    
#   def showPhantomInfo(self):
#     '''Shows phantom information and image'''
#     if hasattr(self,"InfoWindow")==False:
#         self.InfoWindow = Info.InfoWindow()
#     self.InfoWindow.show()
#     self.InfoWindow.setWindowTitle(self.Phantom.phantomName )
#     self.InfoWindow.ui.lblInfo.setPixmap(QPixmap('..\icons\MRISystemPhantom.jpg'))
#     self.InfoWindow.ui.lblInfo.show()
        
  def showPhantomProperties(self):
    '''Shows phantom image'''
    if hasattr(self,"PhantomProperties")==False:
        self.PhantomProperties = PhantomProperties.PhantomProperties(self.Phantom)
    self.PhantomProperties.show()
    self.PhantomProperties.setWindowTitle(self.Phantom.phantomName )
    self.PhantomProperties.ui.lblPhantomImage.setPixmap(QPixmap(self.Phantom.phantomImage))


#*****************************ImageFile  Commands**************************               
  def openDICOMdir (self,filename):
    """Opens DICOMDIR and selected image series"""
    d1= []  #d1 is 3d data stack for 3d images
    dcmdir = read_dicomdir(str(filename))
    if hasattr(self,"DICOMDIRGui")==False:
      self.DICOMDIRlist=DICOMDIRlist.DICOMDIRlist(self)
      self.DICOMDIRlist.setWindowTitle("DICOMDIR list" )
    self.DICOMDIRlist.ui.lblDICOMDIR.setText(filename)
    dv=self.DICOMDIRlist.ui.listDICOMDIR
    for patrec in dcmdir.patient_records:
      s = '' #"Patient: {0.PatientID}: {0.PatientsName}".format(patrec)
      studies = patrec.children
      for study in studies:
            s= s + "    Study {0.StudyID}: {0.StudyDate}".format(study)
            try:
              s= s + "    Study description {0.StudyDescription}".format(study)
            except:
                pass
            dv.addItem(s)
            all_series = study.children
            for series in all_series:
                nImages = len(series.children)
                plural = ('', 's')[nImages > 1]
                if not 'SeriesDescription' in series:
                    series.SeriesDescription = "N/A"
                s= "Series={0.SeriesNumber},  {0.Modality}: {0.SeriesDescription}"  " ({1} image{2})".format(series, nImages, plural)
                dv.addItem(s)
                image_records = series.children
                image_filenames = [os.path.join(self.imageDirectory, *image_rec.ReferencedFileID) for image_rec in image_records]
      if self.DICOMDIRlist.exec_():   #dialog to return list of selected DICOM series to open
        nSeries=self.DICOMDIRlist.selectedSeries()
      else:
        return
      nImages = 0
      for study in studies:
          all_series = study.children
          for series in all_series:
            if str(series.SeriesNumber) in nSeries:
              nImages += len(series.children)
              image_records = series.children
              image_filenames = [os.path.join(self.imageDirectory, *image_rec.ReferencedFileID) for image_rec in image_records]
              self.seriesFileNames.extend(image_filenames)
              nFiles= len(image_filenames)
              dialogBox = QProgressDialog(labelText = 'Importing Files...',minimum = 0, maximum = nFiles)
              dialogBox.setCancelButton(None)
              dsets= [pydicom.read_file(image_filename)for image_filename in image_filenames]
              for i,dcds in enumerate(dsets):   #unpack DICOM data sets (dcds)
                self.ds.unpackImageFile (dcds, image_filenames[i], "dcm")
                #d1.append(self.ds.PA[i+1])
                dialogBox.setValue(i)
    self.nImages += nImages
    self.ui.lblnImages.setText(str(self.nImages))
    self.ui.vsImage.setMinimum(1)       #set slider to go from 1 to the number of images
    self.ui.vsImage.setMaximum(self.nImages)
    self.nCurrentImage=1
    self.ui.vsImage.setValue(self.nCurrentImage)
    self.displayCurrentImage()
    ##self.image3D= np.dstack(d1)
    self.imswin.show()
                                             
  def openFile (self):
    '''opens image file, set of highlighted files, or a DICOM directory'''
    f = QFileDialog.getOpenFileNames(self,"Open Image Files  or DICOMDIR",self.imageDirectory )
    if not f:  #if cancel is pressed return
      return None     
    if type(f)==tuple:    #passes  string with PyQt4 and a tuple with PyQt5
      self.fileNames=f[0]
    else:
      self.fileNames=f
    self.imageDirectory=os.path.dirname(str(self.fileNames[0]))   #Save current directory   
    if len(self.fileNames) == 1:    #check to see if file is a DICOMDIR
      filename=self.fileNames[0]
      if "DICOMDIR" in filename: 
          self.openDICOMdir(filename)
          return None                                            
    self.seriesFileNames.extend(self.fileNames)     #concatenate new file list with previous file list
    nFiles= len(self.fileNames)
    dialogBox = QProgressDialog( 'Importing Files...',None, 0, nFiles, self)
    dialogBox.setCancelButton(None)
    dialogBox.setWindowModality(Qt.WindowModal)
    for i in range(nFiles):
      fileName = self.fileNames[i]
      fstatus=self.ds.addFile(fileName)
      if fstatus[0]:
        dialogBox.setValue(i+1)
        QApplication.processEvents()
      else:
        self.msgPrint(fstatus[1])   #Could not open or read file
    self.sliceList = sorted(set(self.ds.SliceLocation))   #make an ordered list of slices
    self.nImages=len(self.ds.FileName)-1
    self.ui.lblnImages.setText(str(self.nImages))
    if self.nImages < 1 :
      limage = 0
    else:
      limage = 1
    self.ui.vsImage.setMinimum(limage)       #set slider to go from 1 to the number of images
    self.ui.vsImage.setMaximum(self.nImages)
    self.nCurrentImage=limage
    self.ui.vsImage.setValue(self.nCurrentImage)
    self.displayCurrentImage()
    #self.image3D= self.ds.np3dArray()
    self.ims.win.setWindowTitle(self.fileNames[0])
    self.imswin.show()

  def showImswin(self):
    self.imswin.show()
    
  def writeDicomFiles (self):
    f = QFileDialog.getSaveFileName(parent=None, caption="Dicom File Name")
    if not f:  #if cancel is pressed return
      return None     
    if type(f)==tuple:    #passes  string with PyQt4 and a tuple with PyQt5
      fileName=f[0]
    else:
      fileName=f
    self.ds.writeDicomFiles(str(fileName))   #write current image list in DICOM format to filename+ imagenumber + .dcm

  def writeAnimatedGIF (self):
    f = QFileDialog.getSaveFileName(self, caption="GIF File Name")
    if not f:  #if cancel is pressed return
      return None     
    if type(f)==tuple:    #passes  string with PyQt4 and a tuple with PyQt5
      fileName=f[0]
    else:
      fileName=f
    self.ds.writeAnimatedGIF(fileName)   #write current image list to GIF
             
  def changeROISet (self):
    self.InitialROIs =self.Phantom.ROIsets[self.ui.hsROISet.value()]
    self.ui.lblROISet.setText(self.InitialROIs.ROIName)
    self.resetROIs()
    
 # *******************Image Display Methods*****************************     
  def imageSlider (self):
    self.nCurrentImage=self.ui.vsImage.value()
    self.displayCurrentImage()
      
  def displayCurrentImage (self):
    '''Displays current image as set by self.nCurrentImage and associated header parameters'''
    i=self.nCurrentImage
    self.ui.lblCurrentImage.setText(str(self.nCurrentImage)) 
    self.ui.lblDate.setText(format(self.ds.StudyDate[i])) 
    self.ui.lblDataType.setText(format(self.ds.DataType[i])) 
    self.ui.lblFileName.setText(self.ds.FileName[i]) 
    self.ui.lblManufacturer.setText(self.ds.Manufacturer[i]) 
    self.ui.lblSeries.setText(self.ds.SeriesDescription[i]) 
    self.ui.lblInstitution.setText(self.ds.InstitutionName[i]) 
    self.ui.lblField.setText(str(self.ds.MagneticFieldStrength[i]))
    self.ui.lblReceiveCoil.setText(str(self.ds.ReceiveCoilName[i]))    
    self.ui.lblPatient.setText(str(self.ds.PatientName[i])) 
    self.ui.lblProtocol.setText(str(self.ds.ProtocolName[i])) 
    self.ui.lblBW.setText(str(self.ds.PixelBandwidth[i])) 
    self.TEvaries=not(self.checkEqual(self.ds.TE))
    if self.TEvaries:
      self.T2Analysis()
    self.ui.lblTE.setStyleSheet("background-color: yellow") if self.TEvaries else self.ui.lblTE.setStyleSheet("background-color: white")  
    self.ui.lblTE.setText(str(self.ds.TE[i]))
    self.TRvaries = not(self.checkEqual(self.ds.TR))
    if self.TRvaries:
      self.ui.tabT1.setCurrentIndex(2)
    self.ui.lblTR.setStyleSheet("background-color: yellow") if self.TRvaries else self.ui.lblTR.setStyleSheet("background-color: white")
    self.ui.lblTR.setText(str(self.ds.TR[i])) 
    self.ui.lblColumns.setText(str(self.ds.Columns[i]))  
    self.ui.lblRows.setText(str(self.ds.Rows[i]))
    self.TIvaries = not(self.checkEqual(self.ds.TI))
    if self.TIvaries:
      self.ui.tabT1.setCurrentIndex(0)
    self.ui.lblTI.setStyleSheet("background-color: yellow") if self.TIvaries else self.ui.lblTI.setStyleSheet("background-color: white")    
    self.ui.lblTI.setText(str(self.ds.TI[i])) 
    self.ui.lblSliceThickness.setText("{:.2f}".format(self.ds.SliceThickness[i]))
    self.sliceLocationVaries = not(self.checkEqual(self.ds.SliceLocation))
    self.ui.lblSliceLocation.setStyleSheet("background-color: yellow") if self.sliceLocationVaries else self.ui.lblSliceLocation.setStyleSheet("background-color: white")
    self.ui.lblSliceLocation.setText("{:.2f}".format(self.ds.SliceLocation[i])) 
    self.ui.lblPixelSpacingRow.setText("{:.3f}".format(self.ds.PixelSpacingX[i])) 
    self.ui.lblPixelSpacingCol.setText("{:.3f}".format(self.ds.PixelSpacingY[i]))
    self.FAvaries =  not(self.checkEqual(self.ds.FA))
    if self.FAvaries:
      self.ui.tabT1.setCurrentIndex(1)
    self.ui.lblFA.setStyleSheet("background-color: yellow") if self.FAvaries else self.ui.lblFA.setStyleSheet("background-color: white")    
    self.ui.lblFA.setText(str(self.ds.FA[i]))
    self.ui.lblPhaseEncodeDirection.setText(str(self.ds.InPlanePhaseEncodingDirection[i])) 
    self.ui.lblFoVX.setText("{:.3f}".format(self.ds.FoVX[i]))
    self.ui.lblFoVY.setText("{:.3f}".format(self.ds.FoVY[i]))    
    self.ui.lblbValue.setText("{:.2f}".format(self.ds.bValue[i]))
    self.bvaries=not(self.checkEqual(self.ds.bValue))
    if self.bvaries:
      self.diffusionAnalysis()
      self.ui.lblbValue.setStyleSheet("background-color: yellow")
    else:
      self.ui.lblFA.setStyleSheet("background-color: white")    
    data = self.ds.PA[i]  
    xscale =self.ds.PixelSpacingX[i] if (self.ds.PixelSpacingX[i] > 0.) else 1
    yscale = self.ds.PixelSpacingY[i] if (self.ds.PixelSpacingY[i] > 0.) else 1
    xmin = -self.ds.FoVX[i]/2   #set origin to center of image, need to upgrade to set by DICOM tag
    ymin = -self.ds.FoVY[i]/2    
    self.ui.lblUpperLeft.setText("UL=" + "{:.1f}".format(self.ds.ImagePosition[i][0]) + "," + "{:.1f}".format(self.ds.ImagePosition[i][1]) + "," + "{:.1f}".format(self.ds.ImagePosition[i][2]))
    self.imv.setImage(data,pos = (xmin,ymin), scale = (xscale,yscale),)
    self.headerWindow.text.setText(self.ds.header[i])
#    self.ui.lbldX.setText(str(self.ROItrans[0]))
#    self.ui.lbldY.setText(str(self.ROItrans[1]))
#    self.ui.lbldZ.setText(str(self.ROItrans[2]))
    self.imv.getView().setLabel('bottom',self.DirectionLabel(self.ds.RowDirection[i]),"mm")
    self.imv.getView().setLabel('left',self.DirectionLabel(self.ds.ColumnDirection[i]),"mm")
    self.ui.lblScaleSlope.setText("{:.3e}".format(self.ds.ScaleSlope[i]))
    self.ui.lblScaleIntercept.setText("{:.3e}".format(self.ds.ScaleIntercept[i]))
    self.imv.activateWindow()
  #    self.showROIInfo()  #update ROI display                             

#******************End Image Methods************************************************

#******************ROI Methods******************************************************         
  def toggleShowROIs(self):
      self.bShowROIs =not self.bShowROIs
      self.showROIs()
      
  def showROIs(self, roipen=None, clearrois=None, showrois=None):  
    """Displays ROIs if bShowROI is True or if showrois is true, erase ROIs if false"""
    if clearrois:
      self.bShowROIs=False
    if showrois:
      self.bShowROIs=True
    if self.bShowROIs :
        if roipen==None:
          self.roiPen = pg.mkPen(self.currentROIs.ROIColor, width=1)
        else:
          self.roiPen=roipen
        self.ui.lbldh.setText(str(self.relativeOffseth))
        self.ui.lbldv.setText(str(self.relativeOffsetv))
        self.ui.lblCurrentROIs.setText(str(self.currentROIs.ROIName))
        self.ui.lblnROIs.setText(str(self.currentROIs.nROIs))        
        self.pgROIs = []
        #self.bShowSNRROIs = self.ui.chShowBackgroundROIs.isChecked()
        for roi in self.currentROIs.ROIs:
            r=np.array([roi.Xcenter,roi.Ycenter,roi.Zcenter])
            imCoord=self.GlobaltoRel(r,self.nCurrentImage)
            lab=roi.Label            
            if roi.Label =='nolabel':
              lab=''
            if roi.Label =='index':
              lab=str(roi.Index)
            if roi.Type=="Sphere":
              pgroi=fCircleROI(self,[imCoord[0]-roi.d1/2, imCoord[1]-roi.d1/2], [roi.d1, roi.d1],lab, pen=self.roiPen)  #needs work
            if roi.Type=="Rectangle":
              pgroi=fRectROI(self,[imCoord[0]-roi.dx/2, imCoord[1]-roi.dy/2], [roi.dx, roi.dy],lab,angle=roi.theta, pen=self.roiPen)  #needs work

            pgroi.Index=roi.Index
            self.pgROIs.append(pgroi)
        for roi in self.pgROIs:
            self.imv.getView().addItem(roi)
            if self.bRoiLabel:
              self.imv.getView().addItem(roi.label)
        if self.currentROIs.showBackgroundROI:
            roi=self.Phantom.SNRROIs.ROIs[0]
            r=np.array([roi.Xcenter,roi.Ycenter,roi.Zcenter])
            imCoord=self.GlobaltoRel(r,self.nCurrentImage)
            bgroi=fCircleROI(self,[imCoord[0]-roi.d1/2, imCoord[1]-roi.d1/2], [roi.d1, roi.d1],"SNR", pen=self.snrPen)
            self.imv.getView().addItem(bgroi)
            self.bgROI=bgroi
            self.bSNRROIs=True
        if self.currentROIs.showSNRROI:
            roi=self.Phantom.SNRROIs.ROIs[0]
            r=np.array([roi.Xcenter,roi.Ycenter,roi.Zcenter])
            imCoord=self.GlobaltoRel(r,self.nCurrentImage)
            snrroi=fCircleROI(self,[imCoord[0]-roi.d1/2, imCoord[1]-roi.d1/2], [roi.d1, roi.d1],"SNR", pen=self.snrPen)
            self.imv.getView().addItem(snrroi)
            self.snrROI=snrroi
            self.bSNRROIs=True
    else:   #remove all ROIs from the images
        if hasattr(self,"pgROIs"):
          for roi in self.pgROIs:
            self.imv.getView().removeItem(roi)
            try:
              self.imv.getView().removeItem(roi.label)
            except:
              pass
        if self.bSNRROIs:
          try:
            self.imv.getView().removeItem(self.snrROI)
            self.imv.getView().removeItem(self.bgROI)
          except:
            pass
          self.bSNRROIs=False    
        self.pgROIs = []
        self.ui.lblnROIs.setText("")
        self.ui.lblCurrentROIs.setText("")        

  def showROIInfo(self):
      if self.ui.cbViewROIInfo.isChecked():   
        if hasattr(self,"ROIInfo")==False:
            self.ROIInfo = ROIInfo.ROIInfoWindow()
            self.ROIInfo.setWindowTitle(" ROI Info")
            self.ROIInfo.setWindowFlags(Qt.WindowStaysOnTopHint)
            point        = self.rect().topLeft()
            self.ROIInfo.move(point) #QPoint(self.width(), 0))       
        self.ROIInfo.update(self.currentROIs,self.currentROI, self)
        self.ROIInfo.show()
        
  def editROISet(self): 
    self.bEditROISet = True
    
  def editSingleROI(self): 
    self.bEditROISet = False
    
  def allSlices(self): 
    self.bAllSlices = True
    
  def currentSlice(self): 
    self.bAllSlices = False
    
  def useROIValues(self): 
    self.useROIValues = True 

  def toggleRoiLabels(self):
    self.bRoiLabel = not self.bRoiLabel
    self.redrawROIs()
      
  def useBestGuess(self): 
    self.useROIValues = False
       
  def resetROIs(self):
        self.currentROIs=copy.deepcopy(self.InitialROIs)
        self.theta=0.0
        self.ui.lblPhantomAngle.setText("0")
        self.ui.hsAngle.setValue(0)
        self.relativeOffseth = 0.0  #horizontal and vertical offsets of the ROI positions
        self.relativeOffsetv = 0.0
        self.redrawROIs()
        
  def redrawROIs(self, roipen=None):
        self.bShowROIs = False    #erase and redraw initial ROIs
        self.showROIs()  
        self.bShowROIs = True
        self.showROIs(roipen=roipen)  
           
  def ROIColor(self):  
    col = QColorDialog.getColor()
    self.roiPen.setColor(col)
    self.redrawROIs(roipen=self.roiPen)
    
  def labelColor(self):  
    col = QColorDialog.getColor()
    self.lblColor=col
    self.redrawROIs()
           
  def ROIThin(self):  
    self.roiPen.setWidth(1)
    self.redrawROIs(roipen=self.roiPen)
  def ROIMedium(self):  
    self.roiPen.setWidth(2)
    self.redrawROIs(roipen=self.roiPen)
  def ROIThick(self):  
    self.roiPen.setWidth(3)
    self.redrawROIs(roipen=self.roiPen)
                    
  def openROIFile(self):
        if hasattr(self,"ROIprop")==False:
            self.ROIprop = ROIProperties.ROIProperties(self.Phantom)
            self.ROIprop.setWindowTitle("ROI Properties")      
        self.ROIprop.openROIFile(self.imageDirectory)
        self.bShowROIs=True       
        self.redrawROIs() 
        
  def reflectX(self):
    for roi in self.currentROIs.ROIs:
      roi.Xcenter=-roi.Xcenter
    self.redrawROIs()
    
  def reflectY(self):
    for roi in self.currentROIs.ROIs:
      roi.Ycenter=-roi.Ycenter
    self.redrawROIs()

  def reflectZ(self):
    for roi in self.currentROIs.ROIs:
      roi.Zcenter=-roi.Zcenter
    self.redrawROIs() 

  def rotate90X(self):
    for roi in self.currentROIs.ROIs:
      roi.Xcenter,roi.Ycenter,roi.Zcenter=roi.Xcenter,-roi.Zcenter,roi.Ycenter
    self.redrawROIs()
    
  def rotate90Y(self):
    for roi in self.currentROIs.ROIs:
      roi.Ycenter=-roi.Ycenter
    self.redrawROIs()

  def rotate90Z(self):
    for roi in self.currentROIs.ROIs:
      roi.Xcenter,roi.Ycenter,roi.Zcenter=roi.Ycenter,-roi.Xcenter,roi.Zcenter
    self.redrawROIs() 
    
  def translateROIs(self,tvector,snap, roiindex):
        self.relativeOffseth += tvector[0]
        self.relativeOffsetv += tvector[1]
        self.ui.lbldh.setText(str(self.relativeOffseth))
        self.ui.lbldv.setText(str(self.relativeOffsetv))
        r=self.reltoGlobal(tvector[0],tvector[1],self.nCurrentImage)
        if self.bEditROISet == True:  #translate ROI set
          for roi in self.pgROIs:
            roi.translate(tvector, snap=False, finish=False)
            roi.label.setPos(roi.pos())
          self.currentROIs.translate(r)  
        else:   #translate single ROI with roiindex
          roi = self.pgROIs[roiindex-1]
          roi.translate(tvector, snap=False, finish=False)
          roi.label.setPos(roi.pos())
          self.currentROIs.ROIs[roiindex-1].translate(r)
                
  def rotateROIs(self):
    if hasattr(self,"currentROIs"):
      t=float(self.ui.hsAngle.value())/10 #rotation angle in degrees
      self.ui.lblPhantomAngle.setText("{:.1f}".format(t))
      thetanew=float(t * np.pi / 180.)
      dtheta=thetanew-self.theta
      perpAxis=np.cross(self.ds.RowDirection[self.nCurrentImage], self.ds.ColumnDirection[self.nCurrentImage])
      self.theta=thetanew
      self.currentROIs.rotate(perpAxis, dtheta)
      for i, roi in enumerate(self.pgROIs):   #change position of ROIs in pyqtgraph imageview object
            r=np.array([self.currentROIs.ROIs[i].Xcenter,self.currentROIs.ROIs[i].Ycenter,self.currentROIs.ROIs[i].Zcenter])
            (h,v) = self.GlobaltoRel(r,self.nCurrentImage)
            if self.currentROIs.ROIs[i].Type=="Sphere":
              roi.setPos((h-self.currentROIs.ROIs[i].d1/2,v-self.currentROIs.ROIs[i].d1/2))
              roi.label.setPos(h-self.currentROIs.ROIs[i].d1/2,v-self.currentROIs.ROIs[i].d1/2)
            if self.currentROIs.ROIs[i].Type=="Rectangle": 
              roi.setPos((h-self.currentROIs.ROIs[i].dx/2,v-self.currentROIs.ROIs[i].dy/2))
              roi.setAngle(self.currentROIs.ROIs[i].theta)
            
  def deleteROI(self,roi):
    if roi.Index >> 0:
      del self.currentROIs.ROIs[roi.Index-1]
      for i,roi in enumerate(self.currentROIs.ROIs):    #rename ROI indexes
        roi.Index=i+1
    self.currentROIs.nROIs = len(self.currentROIs.ROIs)
    self.bShowROIs = False
    self.showROIs()        
    self.bShowROIs = True
    self.showROIs() 
    
  def addROI(self,pos): 
    newROI = VPhantom.ROI()
    r=self.reltoGlobal(pos.x(),pos.y(),self.nCurrentImage)
    newROI.Xcenter=r[0]
    newROI.Ycenter=r[1]
    newROI.Zcenter=r[2]
    newROI.d1=10
    newROI.Index=len(self.currentROIs.ROIs) +1
    newROI.Name = "def-" + str(newROI.Index)
    self.currentROIs.ROIs.append(newROI)
    self.bShowROIs = False
    self.showROIs()        
    self.bShowROIs = True
    self.showROIs() 
    
  def reSizeROIs(self):
    if hasattr(self,"currentROIs"):
      size=float(self.ui.hsSize.value())/2
      self.ui.lblROISize.setText("{:.1f}".format(size))
      if self.bEditROISet == True:  #change size in ROI set
        for roi in self.currentROIs.ROIs: #initial ROIs remain unchanged
          roi.d1=size
      else:
        self.currentROI.d1=size
      for roi in self.pgROIs:   #erase ROIs
          self.imv.getView().removeItem(roi)
          self.imv.getView().removeItem(roi.label)  
      self.pgROIs = []
      self.showROIs()   #redraw ROIs 

#**********************End ROI Methods**************************************************
                               
  def savePhantomFile(self):
      f = QFileDialog.getSaveFileName(self, "Save ROI File", self.imageDirectory, ".dat")
      if not f:  #if cancel is pressed return
        return None     
      if type(f)==tuple:    #passes  string with PyQt4 and a tuple with PyQt5
        fileName=f[0]
      else:
        fileName=f
      file= open(fileName, 'w')
      self.Phantom.ROIsets.append(self.currentROIs)
      s = self.Phantom.printROIinfo()
      file.write(s)
      file.close()
      
  def openPhantomFile(self):
      self.Phantom=VPhantom.VPhantom()
      f = QFileDialog.getOpenFileName(self,"Open Phantom File",self.imageDirectory )
      if not f:  #if cancel is pressed return
        return None     
      if type(f)==tuple:    #passes  string with PyQt4 and a tuple with PyQt5
        self.phantomFileName=f[0]
      else:
        self.phantomFileName=f
      self.Phantom.readPhantomFile(self.phantomFileName)
      self.InitialROIs =self.Phantom.ROIsets[0]
      self.ui.lblROISet.setText(self.InitialROIs.ROIName)
      self.ui.hsROISet.setMaximum(self.Phantom.nROISets)
      self.resetROIs()  
               
  def checkEqual(self, lst):    #returns True if all elements (except the 0th element) of the list are equal
    return lst[2:] == lst[1:-1]  

  def clearImages (self):  #Deletes all images except default image at index 1
    self.ds= ImageList.ImageList()                         #list of data sets, can be dicom, tiff, fdf
    del self.seriesFileNames[:]
    self.nCurrentImage=0
    self.nImages=0
    self.ui.txtResults.clear()
    #self.image3D.zeros [1,1,1]
    self.displayCurrentImage()
    self.ui.lblnImages.setText(str(self.nImages))
    self.ui.vsImage.setMaximum(0)   
    self.rdPlot.setLabel('left', "Counts", units='A')
    self.rdPlot.setLabel('bottom', "X Axis", units='s')
    self.rdPlot.setTitle("Raw Data")
    self.rdPlot.clear()
    self.resultsPlot.clear()
    self.resultsPlot.setLabel('bottom', "X Axis", units='s')
    self.resultsPlot.setTitle("Results")
    self.ims.win.setWindowTitle("Image Stack")
#    self.resetROIs()

  def deleteCurrentImage(self):
    if self.nCurrentImage > 0:
      self.ds.deleteImage(self.nCurrentImage)
      self.nImages -= 1
      self.ui.lblnImages.setText(str(self.nImages))
      self.ui.vsImage.setMinimum(1)       #set slider to go from 1 to the number of images
      self.ui.vsImage.setMaximum(self.nImages)
      if self.nImages == 0:
          self.nCurrentImage=0
      else:
          self.nCurrentImage = 1
      self.ui.vsImage.setValue(self.nCurrentImage)
      self.displayCurrentImage()


  def sortOnSliceLoc(self):
    self.ds.sortImageList(self.ds.SliceLocation, reverse=False)
    self.displayCurrentImage()
    self.msgPrint('Images reorderd by slice location\n')
    
  def sortOnSliceLocReverse(self):
    self.ds.sortImageList(self.ds.SliceLocation, reverse=True)
    self.displayCurrentImage()
    self.msgPrint('Images reorderd by slice location\n')
         
  def sortOnTE(self):
    self.ds.sortImageList(self.ds.TE)
    self.displayCurrentImage()
    self.msgPrint('Images reorderd by TE\n')
    self.reverseTEOrder = not self.reverseTEOrder  #toggle order 
         
  def unWrapCurrentImage(self):
    '''unwraps a phase image using skimage.restoration.unwrap'''
    self.operateOnAll=self.ui.cbOperateOnAllImages.isChecked()
    if self.nCurrentImage > 0:
      if self.operateOnAll == False:
        self.ds.PA[self.nCurrentImage] = unwrap(self.ds.PA[self.nCurrentImage],wrap_around_axis_0=False, wrap_around_axis_1=False,wrap_around_axis_2=False)
      if self.operateOnAll == True:
        for i in range(1,len(self.ds.PA)):
          #self.ds.PA[i] = unwrap(self.ds.PA[i],wrap_around_axis_0=False, wrap_around_axis_1=False,wrap_around_axis_2=False)
          self.ds.PA[i] =unwrap_phase(self.ds.PA[i])    #From SKIMAGE should be similar to unwrap
      self.displayCurrentImage()

  def setSlopeOffset(self):
      '''set scaling slope and offset to 1,0 respectively to ignore Phillips scaling protocol'''
      for i in range(1,len(self.ds.PA)):
        self.ds.ScaleSlope[i] = 1.0
        self.ds.ScaleIntercept[i]= 0.0
      self.displayCurrentImage() 
            
  def rotateIm90(self):
      '''rotate all images 90 degrees'''
      for i in range(1,len(self.ds.PA)):
        self.ds.PA[i]=np.rot90(self.ds.PA[i])
      self.displayCurrentImage()
      
  def flipIm0(self):
      '''reflect images across axis 0'''
      for i in range(1,len(self.ds.PA)):
        self.ds.PA[i]=np.flip(self.ds.PA[i], axis=0)
      self.displayCurrentImage() 
  def flipIm1(self):
      '''reflect images across axis 1'''
      for i in range(1,len(self.ds.PA)):
        self.ds.PA[i]=np.flip(self.ds.PA[i], axis=1)
      self.displayCurrentImage()
                      
  def setbgROI(self,i,j):
      '''Finds average background around ROI i for image j'''
      try:
          dx=int(self.bgXSpacing/self.ds.PixelSpacingX[j])
          dy=int(self.bgYSpacing/self.ds.PixelSpacingY[j])
          roi =self.pgROIs[i]
          x=roi.pos()[0]+roi.size()[0]/2
          y=roi.pos()[1]+roi.size()[1]/2
          Xindex = int((x+self.ds.FoVX[j]/2)/self.ds.PixelSpacingX[j]) #if self.ds.PixelSpacingX[self.nCurrentImage] > 0. else Xindex = int(mousePoint.x())
          Yindex = int((y+self.ds.FoVY[j]/2)/self.ds.PixelSpacingY[j]) #if self.ds.PixelSpacingY[self.nCurrentImage] > 0. else Yindex = int(mousePoint.y())
          value=  self.ds.PA[j][Xindex,Yindex] 
          bg=(self.ds.PA[j][Xindex+dx,Yindex] +self.ds.PA[j][Xindex-dx,Yindex]+self.ds.PA[j][Xindex,Yindex+dy]+self.ds.PA[j][Xindex,Yindex-dy])/4
          return bg
      except:
          return 0.0

  
  def findBackgroundAroundROIs(self):
      '''Finds average background around an ROI'''
      dx=int(self.bgXSpacing/self.ds.PixelSpacingX[self.nCurrentImage])
      dy=int(self.bgYSpacing/self.ds.PixelSpacingY[self.nCurrentImage])
      for i, roi in enumerate(self.pgROIs):
        x=roi.pos()[0]+roi.size()[0]/2
        y=roi.pos()[1]+roi.size()[1]/2
        Xindex = int((x+self.ds.FoVX[self.nCurrentImage]/2)/self.ds.PixelSpacingX[self.nCurrentImage]) #if self.ds.PixelSpacingX[self.nCurrentImage] > 0. else Xindex = int(mousePoint.x())
        Yindex = int((y+self.ds.FoVY[self.nCurrentImage]/2)/self.ds.PixelSpacingY[self.nCurrentImage]) #if self.ds.PixelSpacingY[self.nCurrentImage] > 0. else Yindex = int(mousePoint.y())
        value=  self.ds.PA[self.nCurrentImage][Xindex,Yindex] 
        bg=(self.ds.PA[self.nCurrentImage][Xindex+dx,Yindex] +self.ds.PA[self.nCurrentImage][Xindex-dx,Yindex]+self.ds.PA[self.nCurrentImage][Xindex,Yindex+dy]+self.ds.PA[self.nCurrentImage][Xindex,Yindex-dy])/4
        self.msgPrint('ROI' +str(i)+ ' BG= ' + "{:6.2f}".format(bg))
        
  def planeBackgroundSubtract(self):
    '''subtracts a plane defined by 3 points set by the previous mouse clicks'''   
    self.operateOnAll=self.ui.cbOperateOnAllImages.isChecked()
    r1=self.clickArray[-3]
    r2=self.clickArray[-2]
    r3=self.clickArray[-1]
    if self.operateOnAll == False:
      a=np.cross(r2-r1,r3-r1)
      an=np.linalg.norm(a)
      a=a/an
      plane=np.zeros([self.ds.Rows[self.nCurrentImage],self.ds.Columns[self.nCurrentImage]])
      for i in range(int(self.ds.Rows[self.nCurrentImage])):
        for j in range(int(self.ds.Columns[self.nCurrentImage])):
          plane[i,j]=(-a[0]*(i-r1[0])-a[1]*(j-r1[1]))/a[2]+r1[2]
      self.ds.PA[self.nCurrentImage]=self.ds.PA[self.nCurrentImage]-plane
    if self.operateOnAll == True:
        for ii in range(1,len(self.ds.PA)):    #Subtract planar background from all images
          r1[2]=self.ds.PA[ii][int(r1[0]),int(r1[1])] #set value of three in plane vectors
          r2[2]=self.ds.PA[ii][int(r2[0]),int(r2[1])]
          r3[2]=self.ds.PA[ii][int(r3[0]),int(r3[1])]
          a=np.cross(r2-r1,r3-r1)
          an=np.linalg.norm(a)
          a=a/an
          plane=np.zeros([self.ds.Rows[ii],self.ds.Columns[ii]])
          for i in range(int(self.ds.Rows[ii])):
            for j in range(int(self.ds.Columns[ii])):
              plane[i,j]=(-a[0]*(i-r1[0])-a[1]*(j-r1[1]))/a[2]+r1[2]
          self.ds.PA[ii]=self.ds.PA[ii]-plane
    self.displayCurrentImage()  

  def ParabolaBackgroundSubtract(self):   
    r1=self.clickArray[-3]
    r2=self.clickArray[-2]
    r3=self.clickArray[-1]
    z1=self.clickArray[-3][2]
    z2=self.clickArray[-2][2]
    z3=self.clickArray[-1][2]
    r1[2]=0
    r2[2]=0
    r3[2]=0    
    Sav=(z1+z2)/2
    d=np.linalg.norm(np.cross((r2-r1),(r3-r1))/np.linalg.norm(r2-r1))
    a=(Sav-z3)/d**2
    #print ("a= " + str(a))
    #print ("Sav= " + str(Sav))
    #print ("d= " + str(d))
    #print ( str(z1) + str(z2)+str(z3))
#    self.PBGSubtract.setText("r1,r2,r3 = " +str(r1) + str(r2)+ str(r3))
#    self.PBGSubtract.setModal(True)
#    self.PBGSubtract.show()
    parab=np.zeros([self.ds.Rows[self.nCurrentImage],self.ds.Columns[self.nCurrentImage]])
    for i in range(int(self.ds.Rows[self.nCurrentImage])):
        for j in range(int(self.ds.Columns[self.nCurrentImage])):
          r0=np.array([i,j,0])
          d=np.linalg.norm(np.cross((r2-r1),(r0-r1))/np.linalg.norm(r2-r1))
          parab[i,j]=a*d**2
    self.ds.PA[self.nCurrentImage]=self.ds.PA[self.nCurrentImage]+parab
    self.displayCurrentImage()  
 
  def imageSubtraction(self):
    '''Subtracts designaed image from all of the images in a stack'''
    text, ok = QInputDialog.getText(self, 'Image Subtraction', 'Image number to be subtracted:')
    iscale=1.
    i=int(text)
    if ok and i > 0 and i < len(self.ds.PA)-1:
      self.operateOnAll=self.ui.cbOperateOnAllImages.isChecked()
      refImage=iscale*self.ds.PA[i]
      if self.operateOnAll:
        for im, pa in enumerate(self.ds.PA[1:]):  #need to skip first dummy image
          self.ds.PA[im+1]=self.ds.PA[im+1]- refImage 
          if im <4:
            self.msgPrint('subtract ' + str(i) + ' from ' + str(im+1) +'; ')
          else:
            self.msgPrint( '.')
      else:
        self.ds.PA[self.nCurrentImage]=self.ds.PA[self.nCurrentImage]- iscale*self.ds.PA[i]
      self.displayCurrentImage()  
 
  def downSample(self):
    self.operateOnAll=self.ui.cbOperateOnAllImages.isChecked()
    text, ok = QInputDialog.getText(self, 'Down Sample', 'Down Sampling Factor:')
    iscale=1.
    i=int(text)
    if ok and i > 0 and not self.operateOnAll: #and i < len(self.ds.PA)-1:
      self.ds.PA[self.nCurrentImage]=block_reduce(self.ds.PA[self.nCurrentImage], i, func=np.mean)
      self.ds.Columns[self.nCurrentImage]=self.ds.Columns[self.nCurrentImage]/i
      self.ds.Rows[self.nCurrentImage]=self.ds.Rows[self.nCurrentImage]/i
      self.ds.PixelSpacingX[self.nCurrentImage]=self.ds.PixelSpacingX[self.nCurrentImage]*i
      self.ds.PixelSpacingY[self.nCurrentImage]=self.ds.PixelSpacingY[self.nCurrentImage]*i
      self.displayCurrentImage()
    if ok and i > 0 and self.operateOnAll:
      for im, pa in enumerate(self.ds.PA): 
        self.ds.PA[im]=block_reduce(self.ds.PA[im], i, func=np.mean)
        self.ds.Columns[im]=self.ds.Columns[im]/i
        self.ds.Rows[im]=self.ds.Rows[im]/i
        self.ds.PixelSpacingX[im]=self.ds.PixelSpacingX[im]*i
        self.ds.PixelSpacingY[im]=self.ds.PixelSpacingY[im]*i
        self.displayCurrentImage()
        
  def interpolate(self):
    '''scales and interpolates 2d images using PIL.image.resize, changes image array size
    3d image are scaled using scipy.ndimage.zoom'''
    interp= ['bicubic', 'nearest', 'bilinear', 'lanczos']
    self.operateOnAll=self.ui.cbOperateOnAllImages.isChecked()
    scalex,ok = QInputDialog.getDouble(self, "Scale X",'the number of rows will be increased (or decreased) by', 1.0, 0.1, 10.0)
    if not ok:
      return
    scaley,ok = QInputDialog.getDouble(self, "Scale Y",'the number of columns will be increased (or decreased) by', 1.0, 0.1, 10.0)
    if not ok:
      return
    dlg2 =  QInputDialog(self)                 
    dlg2.setInputMode(QInputDialog.TextInput) 
    dlg2.setComboBoxItems(interp)
    dlg2.setLabelText('Interpolation type:')                        
    dlg2.setGeometry(50,50,300,100)            
                 
    if self.ds.PA[self.nCurrentImage].ndim==2:
      ok = dlg2.exec_() 
      itype=dlg2.textValue()
      if itype== 'nearest':
        resamp=Image.NEAREST
      if itype== 'bicubic':
        resamp=Image.BICUBIC
      if itype== 'bilinear':
        resamp=Image.BILINEAR
      if itype== 'lanczos':
        resamp=Image.LANCZOS      
      if ok and not self.operateOnAll: #and i < len(self.ds.PA)-1:
        self.ds.Columns[self.nCurrentImage]=int(self.ds.Columns[self.nCurrentImage]*scalex)
        self.ds.Rows[self.nCurrentImage]=int(self.ds.Rows[self.nCurrentImage]*scaley)
        self.ds.PixelSpacingX[self.nCurrentImage]=self.ds.PixelSpacingX[self.nCurrentImage]/scalex
        self.ds.PixelSpacingY[self.nCurrentImage]=self.ds.PixelSpacingY[self.nCurrentImage]/scaley
        n=int(self.ds.Columns[self.nCurrentImage])
        m=int(self.ds.Rows[self.nCurrentImage])
        resizedArray=np.array(Image.fromarray(self.ds.PA[self.nCurrentImage]).resize(size=(n,m), resample=resamp))
        self.ds.PA[self.nCurrentImage]=resizedArray 
        self.displayCurrentImage()
        self.msgPrint('Scaled current image, scalex={:.2f}, scaley={:.2f} , interpolation type={}, image shape= {},{} \n'.format(scalex,scaley,itype,n,m))
      if ok  and self.operateOnAll:
        for im, pa in enumerate(self.ds.PA):
          if im>0:
            self.ds.Columns[im]=int(self.ds.Columns[im]*scalex)
            self.ds.Rows[im]=int(self.ds.Rows[im]*scaley)
            self.ds.PixelSpacingX[im]=self.ds.PixelSpacingX[im]/scalex
            self.ds.PixelSpacingY[im]=self.ds.PixelSpacingY[im]/scaley
            n=self.ds.Columns[im]
            m=self.ds.Rows[im]
            resizedArray=np.array(Image.fromarray(self.ds.PA[im]).resize(size=(n,m), resample=resamp))
            self.ds.PA[im]=resizedArray 
        self.displayCurrentImage()
        self.msgPrint('Scaled image stack using PIL.Image.resize, scalex={:.2f}, scaley={:.2f} , interpolation type={}, image shape= {},{} \n'.format(scalex,scaley,itype,n,m))
    if self.ds.PA[self.nCurrentImage].ndim==3:
      scalez,ok = QInputDialog.getDouble(self, "Scale Z",'the number of images in satck will be increased (or decreased) by', 1.0, 0.1, 10.0)
      if not ok:
        return
      resizedArray=zoom(self.ds.PA[self.nCurrentImage], (scalez,scalex,scaley),order=3)
      self.ds.PA[self.nCurrentImage]=resizedArray
      s=resizedArray.shape
      self.msgPrint('Scaled image stack using scipy.ndimage.zoom, scalex={:.2f}, scaley={:.2f}, scalez={:.2f} ,interpolation=cubic spline image shape= {},{},{} \n'.format(scalex,scaley,scalez,s[0],s[1],s[2]))
      self.displayCurrentImage()
                   
  def cropToROI(self):
    '''Crops image to the rectangular ROI in ImageView widget'''
    self.operateOnAll=self.ui.cbOperateOnAllImages.isChecked()  
    if not self.operateOnAll:    #crop single image or len(self.ds.PA) < 2
        r=self.imv.roi
        arr = r.getArrayRegion(self.imv.image, self.imv.getImageItem())
        self.ds.PA[self.nCurrentImage] = arr
        self.ds.FoVX[self.nCurrentImage]=arr.shape[0]*self.ds.PixelSpacingX[self.nCurrentImage]
        self.ds.FoVY[self.nCurrentImage]=arr.shape[1]*self.ds.PixelSpacingY[self.nCurrentImage]
        self.ds.Columns[self.nCurrentImage]=arr.shape[1]
        self.ds.Rows[self.nCurrentImage]=arr.shape[0]
        self.msgPrint('Cropped current image to ' + "{:.2f}".format(self.ds.FoVX[self.nCurrentImage]) + ' mm;' + "{:.2f}".format(self.ds.FoVY[self.nCurrentImage]) + 'mm \n')
    if self.operateOnAll: #crop multiple images and len(self.ds.PA)>=2
        r=self.imv.roi
        for i,pa in enumerate(self.ds.PA):
          if i>0:   #skip dummy image at image 0
            arr = r.getArrayRegion(pa, self.imv.getImageItem())
            self.ds.PA[i] = arr
            self.ds.FoVX[i]=arr.shape[0]*self.ds.PixelSpacingX[i]
            self.ds.FoVY[i]=arr.shape[1]*self.ds.PixelSpacingY[i]
            self.ds.Columns[i]=arr.shape[1]
            self.ds.Rows[i]=arr.shape[0]
        self.msgPrint('Cropped image stack to ' + "{:.2f}".format(self.ds.FoVX[1]) + ' mm;' + "{:.2f}".format(self.ds.FoVY[1]) + 'mm \n')
    self.displayCurrentImage()
 
  def duplicateImageStack(self):
    dlg =  QInputDialog(self)                 
    dlg.setInputMode(QInputDialog.TextInput) 
    dlg.setLabelText('Enter stack name')                        
    dlg.setGeometry(50,50,300,100)            
    ok = dlg.exec_() 
    stackname=dlg.textValue()
    ims=copy.deepcopy(self.ds)
    ims.ImageStackName=stackname
    self.stackList.append(ims)
    self.msgPrint('Made copy of current image stack, name= ' + stackname + '\n')
    self.displayCurrentImage()
 
  def changeImageStack(self):
    isNames=[str(x.ImageStackName) for x in self.stackList]
    dlg =  QInputDialog(self)                 
    dlg.setInputMode(QInputDialog.TextInput) 
    dlg.setComboBoxItems(isNames)
    dlg.setLabelText(", ".join(isNames))                        
    dlg.setGeometry(50,50,300,100)            
    ok = dlg.exec_() 
    newdataset=dlg.textValue()
    self.ds=self.stackList[isNames.index(newdataset)]
    self.msgPrint('Changed current image stack to ' + newdataset + '\n')
    self.ims.win.setWindowTitle(newdataset)
    self.displayCurrentImage()
       
  def mouseClicked(self, evt):
    '''Adds points to point array used for analysis and background subtraction, adds or subtracts ROI'''
    try:
      pos = evt[0]  ## using signal proxy turns original arguments into a tuple
      mousePoint = self.imv.view.vb.mapSceneToView(pos.scenePos())
      if abs(mousePoint.x()) < self.ds.FoVX[self.nCurrentImage]/2 and abs(mousePoint.y()) < self.ds.FoVY[self.nCurrentImage]/2:
        Xindex = int((mousePoint.x()+self.ds.FoVX[self.nCurrentImage]/2)/self.ds.PixelSpacingX[self.nCurrentImage]) #if self.ds.PixelSpacingX[self.nCurrentImage] > 0. else Xindex = int(mousePoint.x())
        Yindex = int((mousePoint.y()+self.ds.FoVY[self.nCurrentImage]/2)/self.ds.PixelSpacingY[self.nCurrentImage]) #if self.ds.PixelSpacingY[self.nCurrentImage] > 0. else Yindex = int(mousePoint.y())
        value=  self.ds.PA[self.nCurrentImage][Xindex,Yindex]      
        self.msgPrint ( "Add point: " +  "I= " + str(Xindex) + ", J=" + str(Yindex) + ", value=" + str(value) + ", x(mm)= " + "{:.2f}".format(mousePoint.x()) + ", y(mm)=" + "{:.2f}".format(mousePoint.y()) + "\n")
        pnt= np.array([float(Xindex),float(Yindex),float(value)])
        self.clickArray.append(pnt)
        pnt= np.array([mousePoint.x(),mousePoint.y(),float(value)])
        self.clickArraymm.append(pnt)
        dx=self.clickArraymm[-1][0]-self.clickArraymm[-2][0]
        dy=self.clickArraymm[-1][1]-self.clickArraymm[-2][1]
        self.msgPrint ( "dx(mm)= " +  "{:.2f}".format(dx) + ", dy(mm)= " +  "{:.2f}".format(dy)+ ", dr(mm)= " +  "{:.2f}".format((dx**2+dy**2)**0.5) + '\n')
        if self.ui.rbDeleteROI.isChecked():
          self.deleteROI(self.currentROI)
        if self.ui.rbAddROI.isChecked():
          self.addROI(mousePoint)
    except:
      pass

  def viewDicomHeader (self):        
    if self.ui.rbViewDicomHeader.isChecked():
      self.headerWindow.win.show()
      dh = str(self.ds.header[self.nCurrentImage])
      if dh == '':
          dh="DICOM Header"
      self.headerWindow.text.setText(dh)
    else:
      self.headerWindow.win.hide()
      

  def view3DColor(self):  
    self.view3DColor = QColorDialog.getColor()
    self.View3d()
    
  def view3DTransparency(self):  
    t , ok =  QInputDialog.getInt( None,'Transparency' , 'Enter value 0 (solid) to 10 transparent', value=self.view3DTransparency, min=0, max=10, step=1)
    if ok:
      self.view3DTransparency = t
      self.View3d() 
  
  def view3Dinvert(self):  
    self.view3Dinvert = not self.view3Dinvert
    self.View3d()
    
  def View3d(self):
      '''creates 3d rendering of current image stack'''
      if not hasattr(self,"view3Dwin"):   
        self.view3Dwin = gl.GLViewWidget()
        self.view3Dwin.opts['distance'] = 300
        self.view3Dwin.resize(800,800)
        self.view3Dwin.setWindowTitle('3D View ' )
      self.view3Dwin.show()
      try:
         self.view3Dwin.removeItem(self.image3DVol)
      except:
        pass
      ax = gl.GLAxisItem()
      self.view3Dwin.addItem(ax)
      g = gl.GLGridItem()
      g.scale(10, 10, 10)
#      self.view3Dwin.addItem(g) 
      self.image3D=self.ds.np3dArray()
      data=self.image3D.astype(float) /float(self.image3D.max())  #normalize data to 1
      if self.view3Dinvert :
        data=1-data
      d2 = np.empty(data.shape + (4,), dtype=np.ubyte)
      d2[..., 0] = data * self.view3DColor.red()
      d2[..., 1] = data * self.view3DColor.green()
      d2[..., 2] = data * self.view3DColor.blue()
      d2[..., 3] = (data)**self.view3DTransparency * 255.   #sets transparency  
      d2[:, 0:3, 0:3] = [255,0,0,20]   #draw axes at corner of box 
      d2[0:3, :, 0:3] = [0,255,0,20]
      d2[0:3, 0:3, :] = [0,0,255,20]    
      self.image3DVol=gl.GLVolumeItem(d2)
      self.image3DVol.translate(-128,-128,-128)
      self.view3Dwin.addItem(self.image3DVol)
      #self.view3Dwin.update(self.geometry())      
      #self.view3Dwin.repaint(self.geometry())
        

  def mouseMoved(self,evt): 
    '''mouse move event to move crosshairs and display location and values'''
    pos = evt[0]  ## using signal proxy turns original arguments into a tuple
    if self.imv.view.sceneBoundingRect().contains(pos):
        mousePoint = self.imv.view.vb.mapSceneToView(pos)
        self.ui.lblH.setText("{:.2f}".format(mousePoint.x()))
        self.ui.lblV.setText("{:.2f}".format(mousePoint.y()))
        if abs(mousePoint.x()) < self.ds.FoVX[self.nCurrentImage]/2 and abs(mousePoint.y()) < self.ds.FoVY[self.nCurrentImage]/2:
          Xindex = int((mousePoint.x()+self.ds.FoVX[self.nCurrentImage]/2)/self.ds.PixelSpacingX[self.nCurrentImage]) #if self.ds.PixelSpacingX[self.nCurrentImage] > 0. else Xindex = int(mousePoint.x())
          Yindex = int((mousePoint.y()+self.ds.FoVY[self.nCurrentImage]/2)/self.ds.PixelSpacingY[self.nCurrentImage]) #if self.ds.PixelSpacingY[self.nCurrentImage] > 0. else Yindex = int(mousePoint.y())
          value=  self.ds.PA[self.nCurrentImage][Xindex,Yindex]
          try:      
            self.ui.lblValue.setText("{:.2f}".format(value))
          except:
            pass
          rc= self.reltoGlobal(mousePoint.x(), mousePoint.y(), self.nCurrentImage)
          self.ui.lblX.setText("{:.2f}".format(rc[0]))
          self.ui.lblY.setText("{:.2f}".format(rc[1]))
          self.ui.lblZ.setText("{:.2f}".format(rc[2]))
        self.imv.vLine.setPos(mousePoint.x())
        self.imv.hLine.setPos(mousePoint.y()) 

  def reltoGlobal (self, h,v,n):   #given relative coordinate h,v of image n returns np vector of global coordinates 
    #rc= ((h+self.ds.FoVX[n]/2) * self.ds.RowDirection[n]+(v+self.ds.FoVX[n]/2)*self.ds.ColumnDirection[n])+self.ds.ImagePosition[n]
    rc= (h* self.ds.RowDirection[n]+v*self.ds.ColumnDirection[n])
    return rc

  def GlobaltoRel(self,r,n):    #Given r vector in global coordinates returns h,v in image plane of image n
    h=np.dot(r,self.ds.RowDirection[n])  
    v=np.dot(r,self.ds.ColumnDirection[n])
#    h=np.dot(r-self.ds.ImageCenter[n],self.ds.RowDirection[n])  
#    v=np.dot(r-self.ds.ImageCenter[n],self.ds.ColumnDirection[n])
    return [h,v]
   
  def DirectionLabel (self,Vector): #returns a direction label corresponding to input vector
      Label = ""
      if abs(Vector[0])> 0.01:
          Label = "{:.2f}".format(Vector[0])  + "X "
      if abs(Vector[1])> 0.01:
          Label += "+  " + "{:.2f}".format(Vector[1]) +"Y  "          
      if abs(Vector[2])> 0.01:
          Label += "+  " + "{:.2f}".format(Vector[2]) +"Z"          
      return Label
           
  def plotROIvsIndex(self):
    '''plots ROI average value vs image index in raw data plot window'''
    self.dataType = ""
    self.showRawData()

#******************************Data Extraction**************************************    
  def showRawData(self):
    '''Plots ROI signal vs relevant parameter; outputs data in self.rdx and self.rdy'''
    self.ui.txtResults.clear()
    self.msgPrint (self.imageDirectory + "\n")
    self.msgPrint ('Study date =' + self.ds.StudyDate[self.nCurrentImage])  
    self.msgPrint (', Manufacturer =' + self.ds.Manufacturer[self.nCurrentImage]) 
    self.msgPrint (', Series =' +self.ds.SeriesDescription[self.nCurrentImage]) 
    self.msgPrint (', Institution =' +self.ds.InstitutionName[self.nCurrentImage]) 
    self.msgPrint (', B0(T) =' +str(self.ds.MagneticFieldStrength[self.nCurrentImage]))
    self.msgPrint (', Receive Coil =' +str(self.ds.ReceiveCoilName[self.nCurrentImage]))    
    self.msgPrint (', Protocol =' +str(self.ds.ProtocolName[self.nCurrentImage])) 
    self.msgPrint (', Pixel Bandwidth =' +str(self.ds.PixelBandwidth[self.nCurrentImage])+'\n') 
    self.msgPrint ("Data Type = " + self.dataType) 
    self.msgPrint (", Raw data: " + time.strftime("%c") + os.linesep)
    self.rdPlot.clear()
    self.bAllSlices=self.ui.rbAllSlices.isChecked()
    if self.bAllSlices == True:   #analyze all slices together
      self.reducedImageSet= range(1,len(self.ds.PA))
      self.msgPrint ("Slice locations(mm)=" + str(self.ds.SliceLocation[1:]) + "\n")
    else:   # only analyze images which are at the current slice location
      currentSL = self.ds.SliceLocation[self.nCurrentImage]   
      self.reducedImageSet= [i for i, val in enumerate(self.ds.SliceLocation) if np.isclose(val, currentSL , rtol=1e-04, atol=1e-03)] #[0]] something changed, not returning an array now
      self.msgPrint ("Slice location(mm)=" + "{:6.1f}".format(self.ds.SliceLocation[self.nCurrentImage]) + "\n")   
    rd = np.zeros((len(self.pgROIs),len(self.reducedImageSet)))     #array containing raw data average signal in each ROI for each image
    std = np.zeros((len(self.pgROIs),len(self.reducedImageSet)))    #array containing raw data standard deviation of signal in each ROI for each image
    bgROI = np.zeros((len(self.pgROIs),len(self.reducedImageSet)))    #array containing background signal around each ROI for each signal
# Set independent variable (parameter that is being varied ie TI, TE, TR, b, T etc
# T1 data
    if self.dataType == "T1":
      if self.ui.tabT1.currentIndex() == 0: #T1 Inversion Recovery
        self.rdPlot.setLogMode(x=False,y=True)
        self.rdPlot.setLabel('bottom', "TI(ms)")
        self.rdx = np.array([self.ds.TI[i] for i in self.reducedImageSet])
        self.msgPrint ( "TI(ms)=")
        for ti in self.rdx: 
          self.msgPrint ( "{:12.1f}".format(ti))
      if self.ui.tabT1.currentIndex() == 1: #T1 VFA
        self.rdPlot.setLogMode(x=False,y=False)
        self.rdPlot.setLabel('bottom', "FA(deg)")
        self.rdx = np.array([self.ds.FA[j] for j in self.reducedImageSet])
        self.msgPrint ( "FA(deg)=")
        for fa in self.rdx: 
          self.msgPrint ( "{:12.1f}".format(fa))
      if self.ui.tabT1.currentIndex() == 2: #T1 VTR
        self.rdPlot.setLogMode(x=False,y=False)
        self.rdPlot.setLabel('bottom', "TR(ms)")
        self.rdx = np.array([self.ds.TR[j] for j in self.reducedImageSet])
        self.msgPrint ( "TR(ms)=")
        for tr in self.rdx: 
          self.msgPrint ( "{:12.1f}".format(tr))     
      self.msgPrint (os.linesep)
#T2 Data
    if self.dataType == "T2":
      self.rdPlot.setLogMode(x=False,y=True)
      self.rdPlot.setLabel('bottom', "TE(ms)")
      self.rdx = np.array([self.ds.TE[i] for i in self.reducedImageSet])
      self.msgPrint ( "TE(ms)=")
      for te in self.rdx: 
        self.msgPrint ( "{:12.1f}".format(te))
      self.msgPrint (os.linesep)
#General Data
    if self.dataType == "":
      self.rdPlot.setLogMode(x=False,y=False)
      self.rdPlot.setLabel('bottom', "index")
      self.rdx = np.arange(1,len(self.reducedImageSet)+1)
      self.msgPrint ( "index= ")
      for index in self.rdx: 
        self.msgPrint ( "{:12.1f}".format(index))
      self.msgPrint (os.linesep)
#Diffusion Data
    if self.dataType == "Dif":
      if self.ui.tabDif.currentIndex() == 0: #fit signal vs b-value
        self.ADCmap = False
        self.rdPlot.setLogMode(x=False,y=True)
        self.rdPlot.setLabel('bottom', "b(s/mm^2)")
        self.rdx = np.array([self.ds.bValue[i] for i in self.reducedImageSet])
        self.msgPrint ( "b(s/mm^2)=")
        for b in self.rdx: 
          self.msgPrint ( "{:12.1f}".format(b))
        self.msgPrint (os.linesep)
      if self.ui.tabDif.currentIndex() == 1: #ADC map
        self.ADCmap = True
        self.rdPlot.setLogMode(x=False,y=False)
        self.rdPlot.setLabel('bottom', "ROI")
        self.rdx = np.array([roi.Index for roi in self.currentROIs.ROIs])
#PD Data
    if self.dataType == "PD-SNR":
      self.rdPlot.setLogMode(x=False,y=False)
      self.rdx = np.array([roi.PD for roi in self.currentROIs.ROIs])
      self.msgPrint ( "PD(%)=")
      for pd in self.rdx: #note that the first image is a blank dummy
        self.msgPrint ( "{:12.1f}".format(pd))
      self.msgPrint (os.linesep)
    self.dataHeader=self.ui.txtResults.toPlainText()
#Thermometer Data
    if self.dataType == "LCTherm":
      self.rdPlot.setLogMode(x=False,y=False)
      self.rdx = np.array([roi.Value for roi in self.currentROIs.ROIs])
      self.msgPrint ( "LC Temperature Indicator (C)=")
      for T in self.rdx: #note that the first image is a blank dummy
        self.msgPrint ( "{:12.1f}".format(T))
      self.msgPrint (os.linesep)
    self.dataHeader=self.ui.txtResults.toPlainText()
#Set and Plot raw signal data  
    for i, roi in enumerate(self.pgROIs):
      self.msgPrint ("ROI-" +"{:02d}".format(i+1) + '    ') 
      for j, pa in enumerate([self.ds.PA[k] for k in self.reducedImageSet]):
        array = roi.getArrayRegion(pa,self.imv.getImageItem())  #gets masked voxel array for ROI region   
        rd[i ,j]= (array.mean()-self.ds.ScaleIntercept[self.reducedImageSet[j]])/self.ds.ScaleSlope[self.reducedImageSet[j]] #corrects for scaling in Phillips data
        self.msgPrint ( "{:12.3f}".format(rd[i,j]) )
        std[i ,j]=array.std()       #does not have Phillips scaling
        bgROI[i,j]=self.setbgROI(i,j+1) #i=roi, j=image starts at 1 
      c = self.rgb_to_hex(self.setPlotColor(i))
      if self.dataType in ["T1" , "T2", "Dif", ''] and not self.ADCmap:
        self.rdPlot.plot(self.rdx, rd[i,:],pen=self.setPlotColor(i),symbolBrush=self.setPlotColor(i), symbolPen='w', name=self.currentROIs.ROIs[i].Name)    
        self.ui.txtRdLegend.insertHtml('<font size="5" color=' + c + '>' + u"\u25CF" + '</font>' + self.currentROIs.ROIs[i].Name + '<BR>'  )  #u"\u25CF"  + '<BR>' 
      self.msgPrint (os.linesep)
    self.msgPrint('Standard deviations \n')
    for i, roi in enumerate(self.pgROIs):   #printing out standard deviations
      self.msgPrint ("ROI-" +"{:02d}".format(i+1) + '    ') 
      for j, pa in enumerate([self.ds.PA[k] for k in self.reducedImageSet]):
        self.msgPrint ( "{:12.3f}".format(std[i,j]) )
      self.msgPrint (os.linesep)
    self.msgPrint('ROI background signal \n')
    for i, roi in enumerate(self.pgROIs):   #printing out local background
      self.msgPrint ("ROI-" +"{:02d}".format(i+1) + '    ') 
      for j, pa in enumerate([self.ds.PA[k] for k in self.reducedImageSet]):
        self.msgPrint ( "{:12.3f}".format(bgROI[i,j]) )
      self.msgPrint (os.linesep)      
    if self.dataType == "PD-SNR" or self.dataType == "LCTherm": #raw data is a 1d array signal vs ROI.PD
      for k in range(len(self.reducedImageSet)):
        self.rdPlot.plot(self.rdx, rd[:,k],pen=self.setPlotColor(k),symbolBrush=self.setPlotColor(0), symbolPen='w', name=self.currentROIs.ROIs[0].Name)
    if self.dataType == "Dif" and self.ADCmap: #raw data is a 1d array signal vs ROI.Index
        for k in range(len(self.reducedImageSet)):
          self.rdPlot.plot(self.rdx, rd[:,k],pen=self.setPlotColor(0),symbolBrush=self.setPlotColor(0), symbolPen='w', name=self.currentROIs.ROIs[0].Name)    
    self.rdy = rd   #returns a numpy array of raw data (ave signal in an ROI)
    self.std=std    #returns a numpy array of standard deviation of signal in ROI
    self.background=float(self.ui.lblBackGround.text()) #set background counts for fits
    if self.currentROIs.showBackgroundROI:    #obtain background from signal free region (SNR ROI)
      self.rdBackground=np.zeros(len(self.reducedImageSet))
      for j, pa in enumerate([self.ds.PA[k] for k in self.reducedImageSet]):
        roi=self.bgROI
        array = roi.getArrayRegion(pa,self.imv.getImageItem())
        self.rdBackground[j]= (np.average(array)-self.ds.ScaleIntercept[self.reducedImageSet[j]])/self.ds.ScaleSlope[self.reducedImageSet[j]] 
      self.background=np.average(self.rdBackground)
      self.ui.lblBackGround.setText(str(self.background))
    if self.currentROIs.showSNRROI:    #obtain noise from signal region by subtraction of two images
      self.rdNoise=np.zeros(len(self.reducedImageSet))
      for j, pa in enumerate([self.ds.PA[k] for k in self.reducedImageSet]):
        roi=self.snrROI
        array = roi.getArrayRegion(pa,self.imv.getImageItem())
        self.rdNoise[j]= (np.average(array)-self.ds.ScaleIntercept[self.reducedImageSet[j]])/self.ds.ScaleSlope[self.reducedImageSet[j]] 
      if self.rdNoise.shape[0]==2:
        self.noise=self.rdNoise[1]-self.rdNoise[0]
        self.msgPrint('noise = ' + str(self.noise))
      else:
        self.msgPrint('did not find 2 images for noise = analysis')
    self.line1 = pg.InfiniteLine(pos=np.average(self.rdx),  movable=True, angle=90, label='x={value:0.3f}', 
                labelOpts={'position':0.1, 'color': (200,200,100), 'fill': (200,200,200,50), 'movable': True})
    self.line2 = pg.InfiniteLine( movable=True, angle=0, pen=(0, 0, 200),  hoverPen=(0,200,0), label='y={value:0.3f}', 
                labelOpts={'color': (200,200,100), 'movable': True, 'fill': (200,200,200,50)})
    self.rdPlot.addItem(self.line1)
    self.rdPlot.addItem(self.line2)
       
  def fitData(self):
      self.background=float(self.ui.lblBackGround.text()) #set background counts for fits
      self.messageLogOn=True
      self.messageLog=''
      if self.dataType == "T1":
        self.resultsLogModeX=False
        self.resultsLogModeY=True  
        if self.ui.tabT1.currentIndex() == 0: #T1 Inversion Recovery
          self.fitT1IRData(self.rdx,self.rdy)   #
        if self.ui.tabT1.currentIndex() == 1: #T1 Variable Flip Angle
          self.fitT1VFAData(self.rdx,self.rdy)   #
        if self.ui.tabT1.currentIndex() == 2: #T1 Variable TR repetition time
          self.fitT1VTRData(self.rdx,self.rdy)   # 
      if self.dataType == "T2":
        self.resultsLogModeX=False
        self.resultsLogModeY=True  
        self.fitT2SEData(self.rdx,self.rdy)
      if self.dataType == "PD-SNR":
        self.resultsLogModeX=False
        self.resultsLogModeY=True  
        self.fitPDSNR()    
      if self.dataType == "Dif":
        if not self.ADCmap:
          self.resultsLogModeX=False
          self.resultsLogModeY=False    
          self.fitDifData(self.rdx,self.rdy) 
        if self.ADCmap:
          self.showADCmap()
      if self.dataType == "LCTherm":
        self.resultsLogModeX=False
        self.resultsLogModeY=False  
        self.fitLCTherm(self.rdx,self.rdy[:,0])  
      self.resultsPlot.setLogMode(x=self.resultsLogModeX,y=self.resultsLogModeY) 
      self.fitLog=self.messageLog
      self.messageLogOn=False
          
  def showADCmap(self):
      self.resultsPlot.clear()
      self.ADCvalues=self.rdy[:,0]/1000
      sr = "ADC map" + "\n"
      for i, roi in enumerate(self.pgROIs):
          self.msgPrint ("{:02d}".format(i+1)+ " " )
          sr = sr + "ROI " + "{:02d}".format(i+1)+os.linesep  #build output report text for detailed fitting info
          self.msgPrint ("{:10.3f}".format(self.ADCvalues[i]) )
      self.resultsPlot.plot(np.arange(len(self.pgROIs))+1, self.ADCvalues,pen=self.setPlotColor(7),symbolBrush=self.setPlotColor(7), symbolPen='w')
      self.report =  self.report + self.ui.txtResults.toPlainText() + sr   #add recent results to the beginning of the report 
        
  def fitDifData(self, bValue, data):   
      """Fits Diffusion data, calls fitting routines in DifModel"""
      self.resy=np.zeros((len(self.pgROIs),len(bValue)))
      self.referenceError=np.zeros(len(self.pgROIs))
      self.fitx =np.arange(100) * np.amax(bValue) * 1.1 /100  #generate bValues for fit plot
      self.fity=np.zeros((len(self.pgROIs),100))
      self.Difresults=np.zeros((len(self.pgROIs),DifModel.initializeDifModel ()))  #bValue fitting results, first index = ROI, second =ADC, third = S0
      self.Difstderror=np.zeros((len(self.pgROIs),DifModel.initializeDifModel ()))  #array for standard error from fits
      self.msgPrint ("Diffusion fit summary; ADC in mm2/s" + os.linesep)
      self.msgPrint ("ROI  ADC*1000  ADCerr(%)        Si  Sierr(%)" )
      if self.showReferenceValues:
        self.msgPrint ("  ADCref  ADCdev(%)" + os.linesep)
      else:
        self.msgPrint (os.linesep)
      sr = "Diffusion fitting details" + "\n"
      for i, roi in enumerate(self.pgROIs):
          params=DifModel.initializeDifModel (i,bValue, data[i,:],self.currentROIs.ROIs[i], self.useROIValues, self.background)
          pdicti=params[0] #parameter dictionary
          plist=params[1] #parameter list   
          out = lmfit.minimize(DifModel.DifModel,pdicti,args=(bValue,data[i,:]))
          pdict=out.params
          self.resy[i:]=-out.residual #seems to be defined negative to what is logical
          self.fity[i:]= DifModel.DifModel(pdict, self.fitx, np.zeros(len(self.fitx)))
          self.msgPrint ("{:02d}".format(i+1)+ " " )
          sr = sr + "ROI " + "{:02d}".format(i+1)+os.linesep  #build output report text for detailed fitting info
          for p in plist:
            self.Difresults[i,plist.index(p)]=pdict[p].value #populate output array
            self.Difstderror[i,plist.index(p)]=pdict[p].stderr #populate output array
            if pdict[p].value>0:  #calculate standard error in percent
              se=100*pdict[p].stderr/pdict[p].value
            else:
              se=0
            if p == "ADC":  #multiplier to set ADCs in 10^3mm^2/s
              m=1000.
            else:
              m=1.
            self.msgPrint ("{:10.3f}".format(pdict[p].value*m) + "  " + "{:6.2f}".format(se)+ "  ")
          if self.showReferenceValues:
            self.msgPrint ("{:8.3f}".format(1000*self.currentROIs.ROIs[i].ADC)+ "  " + "{:8.2f}".format((self.Difresults[i,0]-self.currentROIs.ROIs[i].ADC)/self.currentROIs.ROIs[i].ADC*100))
          self.msgPrint (os.linesep)
          sr += lmfit.fit_report(pdict)+os.linesep   #add complete fitting report to output report string
      self.resultsPlot.clear()
      err = pg.ErrorBarItem(x=np.arange(len(self.pgROIs))+1, y=1000*self.Difresults[:,0], top=1000*self.Difstderror[:,0], bottom=1000*self.Difstderror[:,0], beam=0.5)
      self.resultsPlot.addItem(err)
      self.resultsPlot.plot(np.arange(len(self.pgROIs))+1, 1000*self.Difresults[:,0],pen=self.setPlotColor(i),symbolBrush=self.setPlotColor(i), symbolPen='w')
      self.report =  self.report + self.ui.txtResults.toPlainText() + sr   #add recent results to the beginning of the report 
                  
  def fitT2SEData(self, TE, data):   
      """Fits T2-SE data, calls fitting routines in T1SE"""
      #self.currentROIs = self.Phantom.ROIsets[1]
      self.resy=np.zeros((len(self.pgROIs),len(TE)))
      self.referenceError=np.zeros(len(self.pgROIs))
      self.fitx =np.arange(100) * np.amax(TE) * 1.1 /100  #generate TEs for fit plot
      self.fity=np.zeros((len(self.pgROIs),100))
      self.T2results=np.zeros((len(self.pgROIs),T2SE.initializeT2SE ()))  #TE fitting results, first index = ROI, second index = parameter referred to in T1Params
      self.msgPrint ("T2-SE fit summary" + os.linesep)
      self.msgPrint ("!!! not fitting points below noise floor =  " +  "{:6.2f}".format(self.noisefactor*self.background) + "\n" )      
      self.msgPrint ("ROI  T2(ms)  T2err(%)        Si  Sierr(%)         B   Berr(%) T2ref(ms)  T2dev(%)    nFitPoints" + os.linesep)
      sr = "T2-SE fitting details" + "\n"
      for i, roi in enumerate(self.pgROIs):
          params=T2SE.initializeT2SE (i,TE, data[i,:],self.currentROIs.ROIs[i], self.useROIValues, self.background)
          pdicti=params[0] #parameter dictionary
          plist=params[1] #parameter list
          d= data[i,:] > self.noisefactor*self.background  
          nfitpoints = len(data[i,d])
          try: 
            out = lmfit.minimize(T2SE.T2SE,pdicti,args=(TE[d],data[i,d]))
            pdict=out.params
            s=0
            for j, ac in enumerate(d):
                if ac:
                    self.resy[i,j]= out.residual[s]
                    s +=1
                else:
                    self.resy[i,j]=np.NaN
            self.fity[i:]= T2SE.T2SE(pdict, self.fitx, np.zeros(len(self.fitx)))
            self.msgPrint ("{:02d}".format(i+1)+ " " )
            sr = sr + "ROI " + "{:02d}".format(i+1)+os.linesep  #build output report text for detailed fitting info
            for p in plist:
              self.T2results[i,plist.index(p)]=pdict[p].value #populate output array
              if pdict[p].value>0:  #calculate standard error
                se=100*pdict[p].stderr/pdict[p].value
              else:
                se=0
              self.msgPrint ("{:10.2f}".format(pdict[p].value) + "  " + "{:6.2f}".format(se)+ "  ")
            self.msgPrint ("{:8.2f}".format(self.currentROIs.ROIs[i].T2)+ "  " + "{:8.2f}".format((self.T2results[i,0]-self.currentROIs.ROIs[i].T2)/self.currentROIs.ROIs[i].T2*100) + "      " + str(nfitpoints))
            self.msgPrint (os.linesep)
            sr += lmfit.fit_report(pdict)+os.linesep   #add complete fitting report to output report string
          except:
              self.msgPrint('Cannot fit ROI# ' +str(i))
      self.resultsPlot.clear()
      self.resultsX=np.arange(len(self.pgROIs))+1
      self.resultsY=self.T2results[:,0]
      self.resultsPlot.plot(self.resultsX,self.resultsY ,pen=self.setPlotColor(0),symbolBrush=self.setPlotColor(0), symbolPen='w')
      self.report =self.report +  self.ui.txtResults.toPlainText() + sr   #add recent results to the beginning of the report 
          
  def fitT1VFAData(self,FA,data):
      """Fits T1-VFA data, calls fitting routines in T1VFA"""
      self.T1results=np.zeros((len(self.pgROIs),T1VFA.initializeT1VFA ()))  #T1 fitting results, first index = ROI, second index = parameter referred to in T1Params
      self.resy=np.zeros((len(self.pgROIs),len(FA)))
      self.referenceError=np.zeros(len(self.pgROIs))
      self.fitx =np.arange(100) * np.amax(FA) * 1.1 /100  #generate flip angles for fit plot
      self.fity=np.zeros((len(self.pgROIs),100))
      self.msgPrint ("T1-VFA fit summary" + os.linesep)
      self.msgPrint ("ROI  T1(ms)    T1err(%)    S90  S90err(%)     TR(ms)  TRdev    TE(ms)  TEdev    B  Berr(%)        dFA    dFAerr")
      if self.showReferenceValues:
        self.msgPrint ("  T1ref(ms) T1dev(%)" + os.linesep)
      else:
        self.msgPrint (os.linesep)
      sr = "T1-VFA fitting details" + "\n"
      for i, roi in enumerate(self.pgROIs):
          params=T1VFA.initializeT1VFA (i,FA, data[i,:],TR=np.array([self.ds.TR[j] for j in self.reducedImageSet])[0],TE=np.array([self.ds.TE[j] for j in self.reducedImageSet])[0] ) #note VFA model needs TR
          pdicti=params[0] #parameter dictionary
          plist=params[1] #parameter list   
          out = lmfit.minimize(T1VFA.T1VFA,pdicti,args=(FA*np.pi/180.0,data[i,:]))
          pdict=out.params
          self.resy[i:]=out.residual
          self.fity[i:]= T1VFA.T1VFA(pdict, self.fitx*np.pi/180.0, np.zeros(len(self.fitx)))
          self.msgPrint ("{:02d}".format(i+1))
          sr= sr + ("ROI " + "{:02d}".format(i+1)+os.linesep)  #build output report text
          for p in plist:
            self.T1results[i,plist.index(p)]=pdict[p].value #populate output array
            if pdict[p].value>0:
              se=100*pdict[p].stderr/pdict[p].value
            else:
              se=0
            self.msgPrint ("{:10.3f}".format(pdict[p].value) + "  " + "{:6.2f}".format(se)+ "  ")
          self.referenceError[i]=(self.T1results[i,0]-self.currentROIs.ROIs[i].T1)/self.currentROIs.ROIs[i].T1*100
          self.msgPrint ("{:8.2f}".format(self.currentROIs.ROIs[i].T1)+ "  " + "{:8.2f}".format(self.referenceError[i]))
          self.msgPrint (os.linesep)
          sr += lmfit.fit_report(pdict)+os.linesep
      self.resultsPlot.clear()
      self.resultsX=np.arange(len(self.pgROIs))+1
      self.resultsY=self.T1results[:,0]
      self.resultsPlot.plot(self.resultsX,self.resultsY ,pen=self.setPlotColor(0),symbolBrush=self.setPlotColor(0), symbolPen='w')
      self.report =self.report +  self.ui.txtResults.toPlainText() + sr
 
  def fitT1VTRData(self,TR,data):
      """Fits T1-VTR data, calls fitting routines in T1VTR"""
      self.T1results=np.zeros((len(self.pgROIs),T1VTR.initializeT1VTR ()))  #T1 fitting results, first index = ROI, second index = parameter referred to in T1Params
      self.resy=np.zeros((len(self.pgROIs),len(TR)))
      self.referenceError=np.zeros(len(self.pgROIs))
      self.fitx =np.arange(100) * np.amax(TR) * 1.1 /100  #generate flip angles for fit plot
      self.fity=np.zeros((len(self.pgROIs),100))
      self.msgPrint ("T1-VTR fit summary" + os.linesep)
      self.msgPrint ("ROI  T1(ms)T1err(%)    S0  S0err(%)    TE(ms)  TEerr        FA  FAerr(%)  T1ref(ms) T1dev(%)" + os.linesep)
      sr = "T1-VTR fitting details" + "\n"
      for i, roi in enumerate(self.pgROIs):
          params=T1VTR.initializeT1VTR (i, TR, data[i,:], self.currentROIs.ROIs[i], self.ds, self.reducedImageSet)  #note VTR model needs FA, TE
          pdicti=params[0] #parameter dictionary
          plist=params[1] #parameter list   
          out = lmfit.minimize(T1VTR.T1VTR,pdicti,args=(TR,data[i,:]))
          pdict=out.params
          self.resy[i:]=out.residual
          self.fity[i:]= T1VTR.T1VTR(pdict, self.fitx, np.zeros(len(self.fitx)))
          self.msgPrint ("{:02d}".format(i+1))
          sr= sr + ("ROI " + "{:02d}".format(i+1)+os.linesep)  #build output report text
          for p in plist:
            self.T1results[i,plist.index(p)]=pdict[p].value #populate output array
            if pdict[p].value>0:
              se=100*pdict[p].stderr/pdict[p].value
            else:
              se=0
            self.msgPrint ("{:10.2f}".format(pdict[p].value) + "  " + "{:6.2f}".format(se)+ "  ")
          self.referenceError[i]=(self.T1results[i,0]-self.currentROIs.ROIs[i].T1)/self.currentROIs.ROIs[i].T1*100
          self.msgPrint ("{:8.2f}".format(self.currentROIs.ROIs[i].T1)+ "  " + "{:8.2f}".format(self.referenceError[i]))
          self.msgPrint (os.linesep)
          sr += lmfit.fit_report(pdict)+os.linesep
      self.resultsPlot.clear()
      self.resultsX=np.arange(len(self.pgROIs))+1
      self.resultsY=self.T1results[:,0]
      self.resultsPlot.plot(self.resultsX,self.resultsY ,pen=self.setPlotColor(0),symbolBrush=self.setPlotColor(0), symbolPen='w')
      self.ui.rbResultsY.setChecked()
      self.report =self.report +  self.ui.txtResults.toPlainText() + sr
      
  def fitT1IRData(self,TI,data):
      """Fits T1-IR data, calls fitting routines in T1IRabs"""
      self.resy=np.zeros((len(self.pgROIs),len(TI)))
      self.referenceError=np.zeros(len(self.pgROIs))
      self.fitx =np.arange(100) * np.amax(TI) * 1.1 /100  #generate TIs for fit plot
      self.fity=np.zeros((len(self.pgROIs),100))
      self.T1results=np.zeros((len(self.pgROIs),T1IRabs.initializeT1IRabs ()))  #T1 fitting results, first index = ROI, second index = parameter referred to in T1Params
      self.msgPrint ("T1-IR fit summary" + os.linesep)
      self.msgPrint ("ROI   T1(ms) T1err(%)        Si  Sierr(%)       B    Berr(%)  T1ref(ms)   T1dev(%)" + os.linesep)
      sr = "T1-IR fitting details" + "\n"
      if self.ui.rbDelta.isChecked():
          delta=float(self.ui.leDelta.text())
      else:
          delta=None
      for i, roi in enumerate(self.pgROIs):
          params=T1IRabs.initializeT1IRabs (i,TI, data[i,:],self.currentROIs.ROIs[i], self.useROIValues, delta=delta)
          pdicti=params[0] #parameter dictionary
          plist=params[1] #parameter list
          out = lmfit.minimize(T1IRabs.T1IRabs,pdicti,args=(TI,data[i,:]))
          pdict=out.params
          self.resy[i:]=out.residual
          self.fity[i:]= T1IRabs.T1IRabs(pdict, self.fitx, np.zeros(len(self.fitx)))
          self.msgPrint ("{:02d}".format(i+1)+ " " )
          sr = sr + "ROI " + "{:02d}".format(i+1)+os.linesep  #build output report text for detailed fitting info
          for p in plist:
            self.T1results[i,plist.index(p)]=pdict[p].value #populate output array
            self.msgPrint ("{:10.2f}".format(pdict[p].value) + "  " + "{:6.2f}".format(100*pdict[p].stderr/pdict[p].value)+ "  ")
          if self.showReferenceValues:
            self.referenceError[i]=(self.T1results[i,0]-self.currentROIs.ROIs[i].T1)/self.currentROIs.ROIs[i].T1*100
            self.msgPrint ("{:8.2f}".format(self.currentROIs.ROIs[i].T1)+ "  " + "{:8.2f}".format(self.referenceError[i]))
          self.msgPrint (os.linesep)
          sr += lmfit.fit_report(pdict)+os.linesep   #add complete fitting report to output report string
      self.resultsPlot.clear()
      self.resultsX=np.arange(len(self.pgROIs))+1
      self.resultsY=self.T1results[:,0]
      self.resultsPlot.plot(self.resultsX,self.resultsY ,pen=self.setPlotColor(0),symbolBrush=self.setPlotColor(0), symbolPen='w')
      self.ui.rbResultsY.setChecked(True)
      self.report =self.report +  self.ui.txtResults.toPlainText() + sr   #add recent results to the beginning of the report 


  def makeT1IRMap(self,TI,data3d, baseline=0, max=10000.0):
      """Makes a T1-IR map, calls fitting routines in T1IRmap"""
      t1map=np.zeros((data3d.shape[1],data3d.shape[2]))
      for i in range(data3d.shape[1]):
        for j in range(data3d.shape[2]):
          if np.amax(data3d[:,i,j]) >baseline:
            params=T1IRmap.init (TI, data3d[:,i,j])
            pdicti=params[0] #parameter dictionary
            plist=params[1] #parameter list
            out = lmfit.minimize(T1IRmap.objFunction,pdicti,args=(TI,data3d[:,i,j]))
            t1map[i,j]= out.params['T1'].value
          else:
            t1map[i,j]=np.nan
      t1map[t1map > max] = np.nan    #set large values to zero
      return t1map

  def makeT2SEMap(self,TE,data3d, baseline=0, max=10000.0):
      """Makes a T2-SE map, calls fitting routines in T2SEmap"""
      t2map=np.zeros((data3d.shape[1],data3d.shape[2]))
      for i in range(data3d.shape[1]):
        for j in range(data3d.shape[2]):
          if np.amax(data3d[:,i,j]) >baseline:
            params=T2SEmap.init (TE, data3d[:,i,j])
            pdicti=params[0] #parameter dictionary
            plist=params[1] #parameter list
            out = lmfit.minimize(T2SEmap.objFunction,pdicti,args=(TE,data3d[:,i,j]))
            t2map[i,j]= out.params['T2'].value
          else:
            t2map[i,j]=np.nan
      t2map[t2map > max] = np.nan    #set large values to zero
      return t2map

  def makeADCMap(self,bvalue,data3d, baseline=0, max=10000.0):
      """Makes a diffusion ADC map, calls fitting routines in DifMap"""
      map=np.zeros((data3d.shape[1],data3d.shape[2]))
      for i in range(data3d.shape[1]):
        for j in range(data3d.shape[2]):
          if np.amax(data3d[:,i,j]) >baseline:
            params=DifMap.init (bvalue, data3d[:,i,j])
            pdicti=params[0] #parameter dictionary
            plist=params[1] #parameter list
            out = lmfit.minimize(DifMap.objFunction,pdicti,args=(bvalue,data3d[:,i,j]))
            map[i,j]= out.params['ADC'].value
          else:
            map[i,j]=np.nan
      map[map > max] = np.nan    #set large values to zero
      return map 
       
  def fitPDSNR(self):
    self.msgPrint('PD-SNR')
    
  def fitLCTherm(self,Tc,data):
      """Fits LCTherm data"""
      self.fitx =np.amin(Tc)+np.arange(100) * (np.amax(Tc)-np.amin(Tc)) /100  #generate Ts for fit plot
      self.fity=np.zeros(100)
      self.msgPrint ("LC Thermometer fit summary" + os.linesep)
      params=LCPowerLawFit.initialize (Tc,data)
      pdicti=params[0] #parameter dictionary
      plist=params[1] #parameter list
      out = lmfit.minimize(LCPowerLawFit.model,pdicti,args=(Tc,data))
      pdict=out.params
      self.msgPrint('T={:4.2f}, gamma={:4.3f}'.format(pdict['T'].value,pdict['gamma'].value))
      self.Temperature=float(pdict['T'].value)

      self.ui.lcdTemp.display(self.Temperature)
      self.fity= LCPowerLawFit.model(pdict, self.fitx, np.zeros(len(self.fitx)))
      self.resultsPlot.clear()
      self.resultsPlot.plot(Tc,data ,pen=self.setPlotColor(0),symbolBrush=self.setPlotColor(0), symbolPen='w')
      self.resultsPlot.plot(self.fitx,self.fity ,pen=self.setPlotColor(1))

  def replotResults(self):
      self.resultsPlot.clear()
      if self.ui.rbResultsY.isChecked():
          self.setDataType(self.dataType)
          self.resultsPlot.setLogMode(x=self.resultsLogModeX,y=self.resultsLogModeY)
          self.resultsPlot.plot(self.resultsX,self.resultsY ,pen=self.setPlotColor(0),symbolBrush=self.setPlotColor(0), symbolPen='w')
      if self.ui.rbResultsError.isChecked():
          self.resultsPlot.setLogMode(x=False,y=False)
          self.resultsPlot.setLabel('left', 'Deviation from Reference (%)') 
          self.resultsPlot.plot(self.resultsX,self.referenceError ,pen=self.setPlotColor(0),symbolBrush=self.setPlotColor(0), symbolPen='w')
           
  def viewFits(self):
        """Opens FitPlot window to display quality of fits"""
        if hasattr(self,"viewFitsWindow")==False:
          self.viewFitsWindow=FitPlots.FitPlots(self.rdx, self.rdy, self.fitx, self.fity, self.resy,header=self.dataHeader)
          self.viewFitsWindow.setWindowTitle("View Fits:" + self.dataType + "  analysis")
        else:
          self.viewFitsWindow.__init__(self.rdx, self.rdy, self.fitx, self.fity, self.resy,header=self.dataHeader)
          self.viewFitsWindow.setWindowTitle("View Fits:" + self.dataType + "  analysis")
        self.viewFitsWindow.show()

  def makeMap(self):
    self.map=mapWindow(self)   #define separate window object  for the map
    #self.map.win      #image stack window
    #self.map.imv    #image stack pyqtgraph image view object
    self.map.win.setWindowTitle('Generating ' + self.dataType + ' Map ...........')
    self.map.win.show()
    self.bAllSlices=self.ui.rbAllSlices.isChecked()
    if self.bAllSlices == True:   #analyze all slices together
      self.reducedImageSet= range(1,len(self.ds.PA))
      self.msgPrint ("Fitting all Slice locations(mm)=" + str(self.ds.SliceLocation[1:]) + "\n")
    else:   # only analyze images which are at the current slice location
      currentSL = self.ds.SliceLocation[self.nCurrentImage]   
      self.reducedImageSet= [i for i, val in enumerate(self.ds.SliceLocation) if np.isclose(val, currentSL , rtol=1e-04, atol=1e-03)] #[0]] something changed, not returning an array now
      self.msgPrint ("Fitting slice location(mm)=" + "{:6.1f}".format(self.ds.SliceLocation[self.nCurrentImage]) + "\n")   
#     text, ok = QInputDialog.getText(self, 'Image baseline', 'Noise level(voxels<noise level will not be fit)')
#     if ok:
#         baseline=float(text)
#     else:
#         baseline=0
    max=np.zeros(len(self.reducedImageSet))
    for i, image in enumerate(self.reducedImageSet):
        max[i]= np.amax(self.ds.PA[image])
    baseline = 0.05*np.amax(max)    
    mapMaking = threading.Thread(target=self.mapMaking(baseline=baseline,dataType=self.dataType))
    mapMaking.start()
    
  def plotMap(self):
    self.map.imv.setImage(self.mapArray)
    self.map.win.setWindowTitle(self.mapTitle)

  def addMaptoStack(self):
    self.ds.addImage(self.mapArray)
    self.nImages += 1
    self.ui.lblnImages.setText(str(self.nImages))
    self.ui.vsImage.setMinimum(1)       #set slider to go from 1 to the number of images
    self.ui.vsImage.setMaximum(self.nImages)
    self.nCurrentImage = self.nImages
    self.ui.vsImage.setValue(self.nCurrentImage)
    self.displayCurrentImage()
        
  def mapMaking(self, baseline=0, dataType=''):
    '''makes T1, T2, ADC maps from current image stack, will not fit any voxels with data below baseline'''
    # T1 data
    if dataType == "T1":
      if self.ui.tabT1.currentIndex() == 0: #T1 Inversion Recovery
        ti = np.array([self.ds.TI[i] for i in self.reducedImageSet])
        data3d = np.array([self.ds.PA[i] for i in self.reducedImageSet])
        self.mapTitle='T1 IR Map'
        self.mapArray=self.makeT1IRMap(ti,data3d, baseline=baseline)
      if self.ui.tabT1.currentIndex() == 1: #T1 VFA
        fa = np.array([self.ds.FA[i] for i in self.reducedImageSet])
        data3d = np.array([self.ds.PA[i] for i in self.reducedImageSet])
        self.mapTitle='T1 FA Map'
        self.mapArray=self.makeT1IRMap(fa,data3d, baseline=baseline)
      if self.ui.tabT1.currentIndex() == 2: #T1 VTR
        self.rdx = np.array([self.ds.TR[j] for j in self.reducedImageSet])
        self.msgPrint ( "T1VTR map")     
    # T2 data
    if dataType == "T2":
      te = np.array([self.ds.TE[i] for i in self.reducedImageSet])
      data3d = np.array([self.ds.PA[i] for i in self.reducedImageSet])
      self.mapTitle='TE SE Map'
      self.msgPrint ( "T2 SE map")
      self.mapArray=self.makeT2SEMap(te,data3d, baseline=baseline)
    # Diffusion data
    if dataType == "Dif":
      bvalue = np.array([self.ds.bValue[i] for i in self.reducedImageSet])
      data3d = np.array([self.ds.PA[i] for i in self.reducedImageSet])
      self.mapTitle='ADC Map'
      self.msgPrint ( "ADC map")
      self.mapArray=self.makeADCMap(bvalue,data3d, baseline=baseline)
    self.msgPrint (os.linesep)
    self.plotMap()

  def riAnalysis(self):
    '''make sub array with resolution inset, with all non-close voxels set to zero'''
       #make a sub array with dimension in mm of riview
    lroi=self.currentROIs.ROIs[-1]  #last roi with rectangular bounding box
    riview=lroi.dx/2
    fovx=self.ds.FoVX[self.nCurrentImage]
    fovy=self.ds.FoVY[self.nCurrentImage]
    dx=self.ds.PixelSpacingX[self.nCurrentImage]
    dy=self.ds.PixelSpacingY[self.nCurrentImage]
    ic=int((lroi.Xcenter+fovx/2)/dx)
    jc=int((lroi.Zcenter+fovy/2)/dy)
    il=ic-int(riview/dx)
    ir=ic+int(riview/dx)
    jl=jc-int(riview/dy)
    jr=jc+int(riview/dy)
    im=self.ds.PA[self.nCurrentImage][il:ir,jl:jr].astype(float)
    fovx=im.shape[0]*self.ds.PixelSpacingX[self.nCurrentImage]    #field of view of cropped image
    fovy=im.shape[1]*self.ds.PixelSpacingY[self.nCurrentImage]
    self.riWin=Resolution.resolutionWindow()    #create an image widow for syntheitc images
    self.riWin.win.show()
    self.riWin.addImage(im,self.currentROIs,fovx=fovx,fovy=fovy, ofovx=self.ds.FoVX[self.nCurrentImage],ofovy=self.ds.FoVY[self.nCurrentImage])


  def sliceProfileSetup(self):
    '''slice profile analysis'''
    self.SystemPhantom()    #load system phantom and select fiduial ROIs
    self.ui.hsROISet.setValue(5)
    self.changeROISet
    self.dataType = "SP"
    self.setDataType(self.dataType)
    #self.ui.tabSliceProfile.setCurrentIndex(0)
    self.displayCurrentImage()
    
  def sliceProfileAnalysis(self):
    nav=11
    nSGorder=5
    pens=[pg.mkPen(self.currentROIs.ROIs[0].color, width=2),pg.mkPen(self.currentROIs.ROIs[1].color, width=2)]
    self.rdPlot.clear()
    self.rdPlot.setTitle('Ramp Profile')
    self.rdPlot.setLabel('bottom', "Distance Along Slice(mm)")
    self.rdPlot.setLabel('left', "Ave signal")
    self.resultsPlot.clear()
    self.resultsPlot.setLabel('bottom', "Distance Perpenduclar to Slice(mm)")
    self.resultsPlot.setLabel('left', "Ave signal")
    self.resultsPlot.setTitle('Slice Profile')
    alpha=SystemPhantom.SPangle
    dx=self.ds.PixelSpacingX[self.nCurrentImage]    #assume in plane pixels are isotropic
    self.spWidth=np.zeros(2)
    self.sliceProfiles=[]   #make a list of slice profile arrays
    nav=self.ui.sbSliceProfilePoints.value()
    nSGorder=self.ui.sbSliceProfileFittingOrder.value()
    for i, roi in enumerate(self.pgROIs):
      roi.setPen(pens[i])
      arr=roi.getArrayRegion(self.ds.PA[self.nCurrentImage], self.imv.getImageItem())
      ls=np.average(arr,axis=0)
      ls=(ls-np.amin(ls))   #noramlize line scan so it goes between 0 and 1
      ls/=np.amax(ls)
      idx = (np.abs(ls - 0.5)).argmin()   #find element whose value is closest to 0.5
      x=(np.arange(ls.shape[0])-idx)*dx
      self.rdPlot.plot(x,ls,pen=pens[i])    #plot linescans
      sp=signal.savgol_filter(ls,nav,nSGorder,deriv=1, mode='nearest')    #take first derivative with SG averaging
      if np.average(sp) <0: #make sure both slice profiles are positive
        sp=-sp
      spmax=np.mean(sp[(sp>0.95*np.amax(sp))])   #Take the peak height as the average of all points greater than 90% of max
      self.spWidth[i]=dx*np.tan(alpha)/spmax    #since the integral of sp = 1, the effective width will be 1/height
      y=x*np.tan(alpha)
      self.sliceProfiles.append(np.vstack((y,sp)))
      if self.ui.rbShowRawProfile.isChecked():
        self.resultsPlot.plot(y,sp,pen=pens[i])
    wDif=self.spWidth[1]-self.spWidth[0]
    wSum=self.spWidth[1]+self.spWidth[0]
    theta=0.5*np.arcsin(wDif*np.sin(2*alpha/wSum))    #theta is the tilt of the phantom relative to the base of the wedges, should be small if phantom installed as prescribed
    self.msgPrint('Raw slice widths(mm)= {:.3f} , {:.3f}; corrected slice width(mm)={:.3f}, Offset angle (deg)={:.3f} \n'.format(self.spWidth[0],self.spWidth[1],np.mean(self.spWidth),theta*180/np.pi ))
    if self.ui.rbShowCorrectedProfile.isChecked(): #display corrected slice profiles  
      for i in range(len(self.sliceProfiles)):    #find corrected y values
        y=self.sliceProfiles[i][0,:]*np.tan(alpha+theta*(-2*i+1))/np.tan(alpha)
        if i==1:  #reverse first profile so -y means towards top of phantom (superior), +y towards bottom (inferior) 
          y=-y
        sp=self.sliceProfiles[i][1,:]*np.tan(alpha)/np.tan(alpha+theta*(-2*i+1))  #correct amplitude for stretching
        self.resultsPlot.plot(y,sp,pen=pens[i])
                     
  def fiducialAnalysis(self):
    '''setup for fiducial array analysis'''
    self.SystemPhantom()    #load system phantom and select fiduial ROIs
    self.ui.hsROISet.setValue(3)
    self.changeROISet
    #self.sortOnSliceLocReverse()    #need bottom of stack to be positive y, plate 1
    self.nCurrentImage=int(len(self.ds.PA)/2)   #set image to middle of stack
    self.ui.vsImage.setValue(self.nCurrentImage)
    self.displayCurrentImage()
    
  def openFiducialAnalysisWindow(self):
    '''make new window for fiducial array analysis'''
    try:
      self.fidWin=Fiducial.fiducialWindow(pw=self)    #create an image widow for synthetic images
      self.fidWin.win.show()
      self.fidWin.addImage(self.ds,self.currentROIs, ny=self.nCurrentImage)
    except:
      self.msgPrint('Cannot do fiducial analysis with current image stack')
               
  
  def saveData(self):
      '''output data to fileName.csv and fit function to filenameFit.csv'''
      f = QFileDialog.getSaveFileName(parent=None, caption="Report File Name", directory = '', selectedFilter = ".csv")
      if not f:  #if cancel is pressed return
        return None     
      if type(f)==tuple:    #passes  string with PyQt4 and a tuple with PyQt5
        fileName=f[0]
      else:
        fileName=f
      file=str(fileName)
      idx=file.index('.')
      avefile=file[:idx]+'AveSignal'+file[idx:]
      fitfile=file[:idx]+'Fit'+file[idx:]
      resfile=file[:idx]+'Res'+file[idx:]
      stdfile=file[:idx]+'STD'+file[idx:]
      fitParamfile=file[:idx]+'FitParam'+file[idx:]
      header=self.dataHeader+'time(ms)'
      for i in range(self.rdy.shape[0]):
        header+=',ROI' +str(i+1) + 'Ave'      
      a1=self.rdx.reshape(self.rdx.shape[0],-1)
      a2=np.transpose(self.rdy) 
      np.savetxt(avefile,np.column_stack((a1,a2)) ,fmt='%.6e', header=header, delimiter=',')
    #STD data
      header=self.dataHeader+'time(ms)'
      for i in range(self.rdy.shape[0]):
        header+=',ROI' +str(i+1) + 'STD'      
      a1=self.rdx.reshape(self.rdx.shape[0],-1)
      a2=np.transpose(self.std) 
      np.savetxt(stdfile,np.column_stack((a1,a2)) ,fmt='%.6e', header=header, delimiter=',')
    #Fit data
      header=self.dataHeader+'time(ms)'
      for i in range(self.rdy.shape[0]):
        header+=',ROI' +str(i+1) + 'fit'       
      a1=self.fitx.reshape(self.fitx.shape[0],-1)
      a2=np.transpose(self.fity) 
      np.savetxt(fitfile,np.column_stack((a1,a2)) ,fmt='%.6e', header=header, delimiter=',')
    #Residual data
      header=self.dataHeader+'time(ms)'
      for i in range(self.rdy.shape[0]):
        header+=',ROI' +str(i+1) + 'Residuals'       
      a1=self.rdx.reshape(self.rdx.shape[0],-1)
      a2=np.transpose(self.resy) 
      np.savetxt(resfile,np.column_stack((a1,a2)) ,fmt='%.6e', header=header, delimiter=',')        
    #Fit parameters
      f= open(fitParamfile, 'w')
      f.write(self.fitLog)
      f.close()

  def saveROIarrays(self):
    #save ROI data arrays
      f = QFileDialog.getSaveFileName(parent=None, caption="BAse File Name for ROI Arrays", directory = '', selectedFilter = ".csv")
      if not f:  #if cancel is pressed return
        return None     
      if type(f)==tuple:    #passes  string with PyQt4 and a tuple with PyQt5
        fileName=f[0]
      else:
        fileName=f
      file=str(fileName)
      idx=file.index('.')

      for i, roi in enumerate(self.pgROIs):
        for j, pa in enumerate([self.ds.PA[k] for k in self.reducedImageSet]):
          fext="ROI-" +"{:02d}".format(i+1)+"Im" + "{:02d}".format(j+1)+ self.dataType+ '_'+ "{:.2f}".format(self.rdx[j])
          array = roi.getArrayRegion(pa,self.imv.getImageItem())   
          header=self.dataHeader+ ' ROI' +str(i+1) + " Image" + "{:02d}".format(j+1)      
          lfile=file[:idx]+fext+file[idx:]
          np.savetxt(lfile,array.filled(fill_value=0.0),fmt='%.6e', header=header, delimiter=',')
          
  def viewReport(self):
        if hasattr(self,"reportwindow")==False:
          self.reportwindow=messageWindow()
        self.reportwindow.win.setWindowTitle("Report:" + self.dataType + "  analysis")
        self.reportwindow.text.setText(self.report)        
        self.reportwindow.win.show()

        
  def saveReport(self):
      f = QFileDialog.getSaveFileName(parent=None, caption="Report File Name", directory = '', selectedFilter = ".dat")
      if not f:  #if cancel is pressed return
        return None     
      if type(f)==tuple:    #passes  string with PyQt4 and a tuple with PyQt5
        fileName=f[0]
      else:
        fileName=f
      file= open(fileName, 'w')
      self.Phantom.ROIsets[0]=self.currentROIs
      file.write(self.report)
      file.close()
      
  def clearReport(self):
      self.report = "" 
      if hasattr(self,"reportwindow")==True:
        self.reportwindow.text.clear()   
             
  def setPlotColor(self,i):
      color = (int(np.base_repr(i, base=3, padding=3)[-1])*127, int(np.base_repr(i, base=3, padding=3)[-2])*127,int(np.base_repr(i, base=3, padding=3)[-3])*127)
      return color

  def rgb_to_hex(self,rgb):
    return '#%02x%02x%02x' % rgb
   
  def T1Analysis(self):
      self.dataType = "T1"
      self.setDataType(self.dataType)
          
  def T2Analysis(self):
      self.dataType = "T2"
      self.setDataType(self.dataType)
         
  def PDSNRAnalysis(self):
      self.dataType = "PD-SNR"
      self.setDataType(self.dataType)
      
  def LCThermometerAnalysis(self):
      self.dataType = "LCTherm"
      self.setDataType(self.dataType)      
      
  def RIAnalysis(self):
      self.dataType = "RI"
      self.setDataType(self.dataType)
#       self.RIwin = Resolution.resolutionWindow(ds=self.ds)
#       self.RIwin.win.show()
       
  def diffusionAnalysis(self):
      self.dataType = "Dif"
      self.setDataType(self.dataType)
           
  def setDataType(self,dataType):
    if dataType == "T1":
      self.ui.stackedModels.setCurrentIndex(0)
      self.setWindowTitle(self.wTitle + 'T1 Analysis')
      self.rdPlot.setTitle("T1 raw data")
      self.rdPlot.setLabel('left', "Ave. Counts", units='A')
      self.rdPlot.setLabel('bottom', "TI(ms)")
      self.resultsPlot.setLabel('bottom', "ROI")
      self.resultsPlot.setLabel('left', "T1(ms)")
      self.resultsPlot.setTitle("T1 Results")
      self.roiPen = pg.mkPen('g', width=3)
      self.ROIset = "T1Array"   #determines which ROI set to use via ROIsetdict
    if dataType == "T2":
      self.ui.stackedModels.setCurrentIndex(1)
      self.setWindowTitle(self.wTitle + 'T2 Analysis')
      self.rdPlot.setTitle("T2 raw data")
      self.rdPlot.setLabel('left', "Ave. Counts", units='A')
      self.rdPlot.setLabel('bottom', "TE")
      self.resultsPlot.setLabel('bottom', "ROI")
      self.resultsPlot.setLabel('left', "T2(ms)")
      self.resultsPlot.setTitle("T2 Results")
      self.roiPen = pg.mkPen('r', width=3)
      self.ROIset = "T2Array"   #determines which ROI set to use via ROIsetdict 
    if dataType == "PD-SNR":
      self.ui.stackedModels.setCurrentIndex(2)
      self.setWindowTitle(self.wTitle + 'Proton Density/ SNR Analysis')
      self.rdPlot.setTitle("PD raw data")
      self.rdPlot.setLabel('left', "Ave. Counts", units='A')
      self.rdPlot.setLabel('bottom', "PD(%)E")
      self.resultsPlot.setTitle("PD Results")
      self.roiPen = pg.mkPen('y', width=3)
      self.ROIset = "PDArray"   #determines which ROI set to use via ROIsetdict 
      self.InitialROIs = self.Phantom.ROIsets[2]     #Current ROIs are initial ROIs with global rotation and translation
      self.resetROIs()
      self.useROIValues = True
    if dataType == "Dif":
      self.ui.stackedModels.setCurrentIndex(3)
      self.setWindowTitle(self.wTitle + 'Diffusion Analysis')
      self.rdPlot.setTitle("Diffusion raw data")
      self.rdPlot.setLabel('left', "Ave. Counts", units='A')
      self.rdPlot.setLabel('bottom', "b(s/mm^2)")
      self.resultsPlot.setLabel('bottom', "ROI")
      self.resultsPlot.setLabel('left', "ADC*1000(mm^2/s)")
      self.resultsPlot.setTitle("Diffusion Results")
      self.roiPen = pg.mkPen('b', width=3)
    if dataType == "SP":    #sliceprofile
      self.ui.stackedModels.setCurrentIndex(4)
      self.ROIset = "SP"   #determines which ROI set to use via ROIsetdict
    if dataType == "RI":    #resolution inset
      self.ui.stackedModels.setCurrentIndex(5)
      self.ROIset = "RI"   #determines which ROI set to use via ROIsetdict
    if dataType == "LCTherm":
      self.setWindowTitle(self.wTitle + 'Liquid Crytal Thermometer Analysis')
      self.rdPlot.setTitle("LCTherm raw data")
      self.rdPlot.setLabel('left', "Ave. Counts", units='A')
      self.rdPlot.setLabel('bottom', "Temperature(C)")
      self.resultsPlot.setTitle("Temperature Results")
      self.roiPen = pg.mkPen('y', width=3)
      self.ROIset = "LCTherm"   #determines which ROI set to use via ROIsetdict
      self.ui.stackedModels.setCurrentIndex(6) 

    if dataType == "":
      self.setWindowTitle(self.wTitle)
      self.rdPlot.setLabel('left', "Counts", units='A')
      self.rdPlot.setLabel('bottom', "X Axis")
      self.rdPlot.setTitle("Raw Data")
      self.rdPlot.clear()
      self.resultsPlot.clear()
      self.resultsPlot.setLabel('bottom', "X Axis", units='s')
      self.resultsPlot.setTitle("Results")
      self.roiPen = pg.mkPen('p', width=3)      

  def TEChanged(self):
      if self.TEChangedFlag == True:
        try:
          self.ds.TE[self.nCurrentImage]=float(self.ui.lblTE.toPlainText())
        except:
          pass
      self.TEChangedFlag == True
      
  def TRChanged(self):
      if self.TRChangedFlag == True:
        try:
          self.ds.TR[self.nCurrentImage]=float(self.ui.lblTR.toPlainText())
        except:
          pass
      self.TRChangedFlag == True

  def getDouble(self, text1='Qinput dialog box', text2='inputtext'):
      scale,ok = QInputDialog.getDouble(self, "Resolution","enter standard deviation (mm)", 0, 0, 10)
      if OK:
        return double
      else:
        return None
              
  def msgPrint (self, s):
          self.ui.txtResults.insertPlainText(s)
          if self.messageLogOn== True:
            self.messageLog+=s

class fCircleROI(pg.EllipseROI):   #originally was pg.EllipseROI
    """Defines a circular ROI using pyqtgraph's EllipseROI"""
    def __init__(self, callingform, pos, size, label,   **args):   #passing in calling form, could be a problem
        pg.ROI.__init__(self, pos, size, **args)
        self.aspectLocked = True
        self.PhantomViewerInstance=callingform
        self.Index = 0
        self.label = pg.TextItem(label, callingform.lblColor, anchor = (0,0))
        self.label.setPos(pos[0],pos[1])
        self.label.setFont(callingform.lblFont)
        
        self.path=None    #added when updating to pyqtgraph 0.11 Do not know why it is not there

    def getArrayRegion(self, arr, img=None):
        arr = pg.ROI.getArrayRegion(self, arr, img)
        if arr is None or arr.shape[0] == 0 or arr.shape[1] == 0:
            return None
        w = arr.shape[0]
        h = arr.shape[1]
        ## generate an ellipsoidal mask
        mask = np.fromfunction(lambda x,y: (((x+0.5)/(w/2.)-1)**2+ ((y+0.5)/(h/2.)-1)**2)**0.5 > 1, (w, h))
        maskedArr = ma.masked_array(arr, mask)
        return maskedArr

    def mouseDragEvent(self, ev):       #Dragging ROI event: translates ROIs
        if ev.isStart():  
            if ev.button() == Qt.LeftButton:
                self.setSelected(True)
                if self.translatable:
                    self.isMoving = True
                    self.preMoveState = self.getState()
                    self.cursorOffset = self.pos() - self.mapToParent(ev.buttonDownPos())
                    self.sigRegionChangeStarted.emit(self)
                    ev.accept()
                else:
                    ev.ignore()
        elif ev.isFinish():
            if self.translatable:
                if self.isMoving:
                    self.stateChangeFinished()
                self.isMoving = False
            return
        if self.translatable and self.isMoving and ev.buttons() == Qt.LeftButton:
            snap = True if (ev.modifiers() & Qt.ControlModifier) else None
            newPos = self.mapToParent(ev.pos()) + self.cursorOffset
#            self.translate(newPos - self.pos(), snap=snap, finish=False)
            self.PhantomViewerInstance.translateROIs(newPos - self.pos(),snap, self.Index)

    def setMouseHover(self, hover):
        '''Inform the ROI that the mouse is(not) hovering over it'''
        if self.mouseHovering == hover:
            return
        self.mouseHovering = hover
        if hover:
            self.currentPen = fn.mkPen(255, 255, 0)
            array = self.getArrayRegion(self.PhantomViewerInstance.ds.PA[self.PhantomViewerInstance.nCurrentImage],self.PhantomViewerInstance.imv.getImageItem())
            mean=array.mean()
            std=array.std()
            self.PhantomViewerInstance.currentROI=self.PhantomViewerInstance.currentROIs.ROIs[self.Index-1]           
            self.PhantomViewerInstance.currentROI.SignalAve= mean
            self.PhantomViewerInstance.currentROI.SignalRMS= std
            self.PhantomViewerInstance.currentROI.array= array
            self.PhantomViewerInstance.showROIInfo()
        else:
            self.currentPen = self.pen
            #self.PhantomViewerInstance.currentROI=-1
        self.update() 
        
class fRectROI(pg.RectROI):
    """Defines a rectangular ROI using pyqtgraph's RectROI"""
    def __init__(self, callingform, pos, size, label,   **args):   
        pg.ROI.__init__(self, pos, size, **args)
        self.aspectLocked = False
        self.PhantomViewerInstance=callingform
        self.Index = 0
        self.label = pg.TextItem(label, callingform.lblColor, anchor = (0,0))
        self.label.setPos(pos[0],pos[1])
        self.label.setFont(callingform.lblFont)

    def mouseDragEvent(self, ev):       #Dragging ROI event: translates ROIs
        if ev.isStart():  
            if ev.button() == Qt.LeftButton:
                self.setSelected(True)
                if self.translatable:
                    self.isMoving = True
                    self.preMoveState = self.getState()
                    self.cursorOffset = self.pos() - self.mapToParent(ev.buttonDownPos())
                    self.sigRegionChangeStarted.emit(self)
                    ev.accept()
                else:
                    ev.ignore()
        elif ev.isFinish():
            if self.translatable:
                if self.isMoving:
                    self.stateChangeFinished()
                self.isMoving = False
            return
        if self.translatable and self.isMoving and ev.buttons() == Qt.LeftButton:
            snap = True if (ev.modifiers() & Qt.ControlModifier) else None
            newPos = self.mapToParent(ev.pos()) + self.cursorOffset
#            self.translate(newPos - self.pos(), snap=snap, finish=False)
            self.PhantomViewerInstance.translateROIs(newPos - self.pos(),snap, self.Index)

    def setMouseHover(self, hover):
        '''Inform the ROI that the mouse is(not) hovering over it'''
        if self.mouseHovering == hover:
            return
        self.mouseHovering = hover
        if hover:
            self.currentPen = fn.mkPen(255, 255, 0)
            array = self.getArrayRegion(self.PhantomViewerInstance.ds.PA[self.PhantomViewerInstance.nCurrentImage],self.PhantomViewerInstance.imv.getImageItem())
            mean=array.mean()
            std=array.std()
            self.PhantomViewerInstance.currentROI=self.PhantomViewerInstance.currentROIs.ROIs[self.Index-1]           
            self.PhantomViewerInstance.currentROI.SignalAve= mean
            self.PhantomViewerInstance.currentROI.SignalRMS= std
            self.PhantomViewerInstance.showROIInfo()
        else:
            self.currentPen = self.pen
            #self.PhantomViewerInstance.currentROI=-1
        self.update() 
        
class messageWindow():
  '''General text window'''
  def __init__(self, wtitle='Message window', image=None, parent = None):
    '''Defines message window, contains output data:
    self.win, self.text are the QMainWIndow and QtextEdit, respectively'''    
    self.win = QMainWindow()
    #self.win.setWindowFlags(Qt.WindowStaysOnTopHint)        #having the window always on top can be a problem
    self.win.resize(1000,500)
    self.text=QTextEdit()
    self.win.setCentralWidget(self.text)
    self.win.setWindowTitle(wtitle)
    self.menu = self.win.menuBar()
    
    self.fileMenu = self.menu.addMenu('&file')
    self.actionSaveFile = QAction('Save to file', self.win)
    self.actionSaveFile.triggered.connect(self.saveFile)
    self.fileMenu.addAction(self.actionSaveFile)
    
  def saveFile(self):
    f = QFileDialog.getSaveFileName(self.win, "File Name", '', "txt(*.txt)")
    if not f:  #if cancel is pressed return
      return None     
    if type(f)==tuple:    #passes  string with PyQt4 and a tuple with PyQt5
      fileName=f[0]
    else:
      fileName=f
    file= open(fileName, 'w')
    file.write(self.text.toPlainText())
    file.close()
        
class imageStackWindow(PhantomViewer):
  def __init__(self, pw, parent = None):
    '''Define image stack windows and menus, pw is the parent PhantomViewer window'''    
    super(PhantomViewer, self).__init__()
    self.win = QMainWindow()
#    self.win.setWindowFlags(Qt.WindowStaysOnTopHint)        #having the window always on top can be a problem
    self.win.resize(800,600)
    self.pw=pw    #parent window to allow access to parent methods
    self.imv = pg.ImageView( view = pg.PlotItem())
    self.win.setCentralWidget(self.imv)
    self.win.setWindowTitle('Image Stack')
    point = self.rect().topRight()
    self.win.move(point + QPoint(self.width()/2, 0))
    self.win.statusBar = QStatusBar()
    self.win.setStatusBar(self.win.statusBar)
    self.win.statusBar.show()

    self.menu = self.win.menuBar()
#    self.menu.setToolTipsVisible(True)
    self.imageMenu = self.menu.addMenu('&ImageFiles')
    self.actionSelectImages = QAction('Select/Add Images', self)
    self.actionSelectImages.setShortcut('Ctrl+S')
    self.actionSelectImages.setStatusTip('Opens DICOM, DicomDir, tif, fdf' )
    self.actionSelectImages.triggered.connect(pw.openFile)
    self.imageMenu.addAction(self.actionSelectImages) 
    self.actionClear_All_Images = QAction('Clear All Images', self)
    self.actionClear_All_Images.setShortcut('Ctrl+C')
    self.actionClear_All_Images.setStatusTip('Clear All Images')
    self.actionClear_All_Images.triggered.connect(pw.clearImages)
    self.imageMenu.addAction(self.actionClear_All_Images)

    self.actionWrite_Animated_GIF = QAction('Write Animated GIF', self)
    self.actionWrite_Animated_GIF.setStatusTip('Write Animated GIF')
    self.actionWrite_Animated_GIF.triggered.connect(pw.writeAnimatedGIF)
    self.imageMenu.addAction(self.actionWrite_Animated_GIF)
    
    self.actionWriteDICOM = QAction('Write DICOM', self)
    self.actionWriteDICOM.setStatusTip('Write DICOM')
    self.actionWriteDICOM.triggered.connect(pw.writeDicomFiles)
    self.imageMenu.addAction(self.actionWriteDICOM)
    
    self.imageStackMenu = self.menu.addMenu('&ImageStacks')
    self.actionSort_on_Slice_Loc = QAction('Sort on Slice Loc (neg at bot)', self)
    self.actionSort_on_Slice_Loc.triggered.connect(pw.sortOnSliceLoc)
    self.actionSort_on_Slice_Loc.setStatusTip('Sort so negative slice at bottom')
    self.imageStackMenu.addAction(self.actionSort_on_Slice_Loc)
    self.actionSort_on_Slice_Loc_Reverse = QAction('Sort on Slice Loc (neg at top)', self)
    self.actionSort_on_Slice_Loc_Reverse.triggered.connect(pw.sortOnSliceLocReverse)
    self.actionSort_on_Slice_Loc_Reverse.setStatusTip('Sort so negative slice at top')
    self.imageStackMenu.addAction(self.actionSort_on_Slice_Loc_Reverse)        
    self.actionDuplicateImageStack = QAction('Duplicate Image Stack', self)
    self.actionDuplicateImageStack.triggered.connect(pw.duplicateImageStack)
    self.imageStackMenu.addAction(self.actionDuplicateImageStack)
    self.actionChangeImageStack = QAction('Change Image Stack', self)
    self.actionChangeImageStack.triggered.connect(pw.changeImageStack)
    self.imageStackMenu.addAction(self.actionChangeImageStack) 
    self.actionRotate90 = QAction('Rotate Images 90 deg', self)
    self.actionRotate90.triggered.connect(pw.rotateIm90)
    self.imageStackMenu.addAction(self.actionRotate90) 
    self.actionFlip0 = QAction('Flip Images about vertical', self)
    self.actionFlip0.triggered.connect(pw.flipIm0)
    self.imageStackMenu.addAction(self.actionFlip0) 
    self.actionFlip1 = QAction('Flip Images about horizontal', self)
    self.actionFlip1.triggered.connect(pw.flipIm1)
    self.imageStackMenu.addAction(self.actionFlip1)         
                   
    self.phantomMenu = self.menu.addMenu('&Phantoms')
    self.actionOpenPhantomFile = QAction('Open phantom file', self)
    self.actionOpenPhantomFile.setStatusTip('Open phantom file')
    self.actionOpenPhantomFile.triggered.connect(pw.openPhantomFile)
    self.phantomMenu.addAction(self.actionOpenPhantomFile)
    self.actionSavePhantomFile = QAction('Save phantom file', self)
    self.actionSavePhantomFile.setStatusTip('Save phantom file')
    self.actionSavePhantomFile.triggered.connect(pw.savePhantomFile)
    self.phantomMenu.addAction(self.actionSavePhantomFile)  
    
    self.actionSystemPhantom = QAction('System Phantom', self)
    self.actionSystemPhantom.setStatusTip('NIST/ISMRM System Phantom')
    self.actionSystemPhantom.triggered.connect(pw.SystemPhantom)
    self.phantomMenu.addAction(self.actionSystemPhantom) 
    self.actionDiffusionPhantom = QAction('Diffusion Phantom', self)
    self.actionDiffusionPhantom.setStatusTip('NIST/RSNA/NCI Diffusion Phantom')
    self.actionDiffusionPhantom.triggered.connect(pw.DiffusionPhantom)
    self.phantomMenu.addAction(self.actionDiffusionPhantom) 
    self.actionBreastPhantom = QAction('Breast Phantom', self)
    self.actionBreastPhantom.setStatusTip('NIST/UCSF/NCI Breast Phantom')
    self.actionBreastPhantom.triggered.connect(pw.BreastPhantom)
    self.phantomMenu.addAction(self.actionBreastPhantom) 
    self.actionNISThcpPhantom = QAction('NIST hcp Phantom', self)
    self.actionNISThcpPhantom.setStatusTip('NIST hcp Phantom')
    self.actionNISThcpPhantom.triggered.connect(pw.NISThcpPhantom)
    self.phantomMenu.addAction(self.actionNISThcpPhantom)
    self.actionNISThcpCoronalPhantom = QAction('NIST hcp Coronal Phantom', self)
    self.actionNISThcpCoronalPhantom.setStatusTip('NIST hcp Coronal Phantom')
    self.actionNISThcpCoronalPhantom.triggered.connect(pw.NISThcpCoronalPhantom)
    self.phantomMenu.addAction(self.actionNISThcpCoronalPhantom) 
    self.actionShowPhantomProperties = QAction('Current phantom info', self)
    self.actionShowPhantomProperties.setStatusTip('Current phantom info')
    self.actionShowPhantomProperties.triggered.connect(pw.showPhantomProperties)
    self.actionShowPhantomProperties.setShortcut('Ctrl+P')
    self.phantomMenu.addAction(self.actionShowPhantomProperties)   
    
    self.roiMenu = self.menu.addMenu('&ROIs')
    self.actionRoiLabel = QAction('Show/hide ROI labels', self)
    self.actionRoiLabel.triggered.connect(pw.toggleRoiLabels)
    self.roiMenu.addAction(self.actionRoiLabel)      
    
#Analysis    
    self.analysisMenu = self.menu.addMenu('&Analysis')
    self.actionRI = QAction('Resolution Inset', self)
    self.actionRI.triggered.connect(pw.riAnalysis)
    self.analysisMenu.addAction(self.actionRI)
    
    self.actionFiducial = QAction('Fiducial Array Setup', self)
    self.actionFiducial.triggered.connect(pw.fiducialAnalysis)
    self.actionFiducial.setStatusTip('Make sure bottom of phantom at bottom of stack, coarse align, click Fiducual Array Analysis')
    self.analysisMenu.addAction(self.actionFiducial) 
          
    self.actionFiducialAnalysis = QAction('Fiducial Array Analysis', self)
    self.actionFiducialAnalysis.triggered.connect(pw.openFiducialAnalysisWindow)
    self.actionFiducialAnalysis.setStatusTip('Open analysis window, ')
    self.analysisMenu.addAction(self.actionFiducialAnalysis)
     
    self.actionSliceProfileSetup = QAction('Slice Profile Setup', self)
    self.actionSliceProfileSetup.triggered.connect(pw.sliceProfileSetup)
    self.actionSliceProfileSetup.setStatusTip('Sets up slice profile anlaysis, Align rectangular ROIs, click Slice Profile Analysis ')
    self.analysisMenu.addAction(self.actionSliceProfileSetup) 
    
    self.actionSliceProfileAnalysis = QAction('Slice Profile Analysis', self)
    self.actionSliceProfileAnalysis.triggered.connect(pw.sliceProfileAnalysis)
    self.actionSliceProfileAnalysis.setStatusTip('Sets up slice profile analysis, Align rectangular ROIs, click Slice Profile Analysis ')
    self.analysisMenu.addAction(self.actionSliceProfileAnalysis)           
    
    self.actionInterpolate = QAction('Scale image', self)
    self.actionInterpolate.setStatusTip('Resizes image with (n,m) rows and columns to int((scale*n,scale*m)); scale > 1 interpolates, scale<1 downsamples')
    self.actionInterpolate.triggered.connect(pw.interpolate)
    self.analysisMenu.addAction(self.actionInterpolate)  
    self.actionCropToROI = QAction('Crop To ROI', self)
    self.actionCropToROI.setStatusTip('Draw ROI')
    self.actionCropToROI.triggered.connect(pw.cropToROI)
    self.analysisMenu.addAction(self.actionCropToROI)
    self.actionFindBackgroundAroundROIs = QAction('FindBackgroundAroundROIs', self)
    self.actionFindBackgroundAroundROIs.setStatusTip('Find Background Around ROIs')
    self.actionFindBackgroundAroundROIs.triggered.connect(pw.findBackgroundAroundROIs)
    self.analysisMenu.addAction(self.actionFindBackgroundAroundROIs)
    
    self.windowMenu = self.menu.addMenu('&Windows')
    self.actionWinAlwaysOnTop = QAction('Window always on top', self)
    self.actionWinAlwaysOnTop.triggered.connect(self.winOnTop)
    self.windowMenu.addAction(self.actionWinAlwaysOnTop)
    self.actionWinNotOnTop = QAction('Window not always on top', self)
    self.actionWinNotOnTop.triggered.connect(self.winNotOnTop)
    self.windowMenu.addAction(self.actionWinNotOnTop)          
    
    self.imv.getView().setLabel('bottom',"H","mm")
    self.imv.getView().setLabel('left',"V","mm")
    
  def winOnTop(self):
        self.win.setWindowFlag(Qt.WindowStaysOnTopHint, True)
        self.pw.imswin.show()
        
  def winNotOnTop(self):
        self.win.setWindowFlag(Qt.WindowStaysOnTopHint, False)
        self.pw.imswin.show()

class mapWindow(PhantomViewer):
  def __init__(self, pw, parent = None):
    '''Define image stack windows and menus, pw is the parent PhantomViewer window'''    
    super(PhantomViewer, self).__init__()
    self.win = QMainWindow()
    #self.win.setWindowFlags(Qt.WindowStaysOnTopHint)        #having the window always on top can be a problem
    self.win.resize(800,600)
    self.imv = pg.ImageView( view = pg.PlotItem())
    self.win.setCentralWidget(self.imv)
    self.win.setWindowTitle('Map')
    point = self.rect().topRight()
    self.win.move(point + QPoint(self.width()/2, 0)) 
    self.menu = self.win.menuBar()
    self.imageMenu = self.menu.addMenu('&Images')
    self.actionAddToStack = QAction('Add Map to stack', self)
    self.actionAddToStack.setStatusTip('AddToStack')
    self.actionAddToStack.triggered.connect(pw.addMaptoStack)
    self.imageMenu.addAction(self.actionAddToStack)

class phantomImageWindow(PhantomViewer):
  def __init__(self, fname, parent = None):
    '''Defines image stack windows and menus, rv is the parent PhantomViewer window'''    
    super(PhantomViewer, self).__init__()
    self.win = QMainWindow()
    self.win.setWindowFlags(Qt.WindowStaysOnTopHint)
    self.win.resize(800,600)
    self.win.setWindowTitle('Current Phantom')
    self.setStyleSheet("background-image: url(" + str(fname) + "); background-repeat: no-repeat; background-position: center;")
    #self.win.setCentralWidget(self.imv)

#Useful for debugging Qt applications where the app closes without giving error message
sys._excepthook = sys.excepthook 
def exception_hook(exctype, value, traceback):
    print("Missed Exception:", exctype, value, traceback)
    sys._excepthook(exctype, value, traceback) 
    sys.exit(1) 
#*******

sys.excepthook = exception_hook                               
if __name__ == '__main__':
    app = QApplication(sys.argv)
    main = PhantomViewer("T1")
    main.show()
    sys.exit(app.exec_())

    