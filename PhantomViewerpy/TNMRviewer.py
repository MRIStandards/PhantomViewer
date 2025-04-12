'''
Created on Nov 7, 2023

@author: russek

Classes to display 4d NMR/MRI data from Tecmag TNMR

Conventions:

all stored parameter units are SI (T, m, s, A), ***displayed units may be mm, ms, mT
Assume tnt data is a 4 dimensional array with RO, slice, phase, parameter), switch RO and phase so the internal working array is slice, RO, phase, parameter.

If the following strings are found in the filename
    'TRISCOUT' then the protocol is interpreted as a TRISCOUT, assumes 3 images in each readout that are separated as 
    'T2SE' then the protocol is interpreted as a T2SE, looks for a table 'teDelay', calculates TE array as self.TEArray=2*(self.teDelay+self.TE0)
    'T1IRSE' then the protocol is interpreted as a T1IRSE, looks for a table 'tiDelay', calculates TI array as self.TEArray=self.tiDelay+self.TI0
    'GE_FLASH' then the protocol is interpreted as a GE_FLASH
    'GEIR_FLASH' then the protocol is interpreted as a GEIR_FLASH

'''
VersionID='TNMRviewer V12-1-2023'
import sys, os, time, datetime, copy        #useful system libraries
import numpy as np
from PyQt5.Qt import PYQT_VERSION_STR
from PyQt5.QtCore import  Qt, QPoint, QTimer
from PyQt5.QtGui import  QFont, QColor, QPainter, QPixmap, QTextOption, QScreen, QPen, QTextCursor
from PyQt5.QtWidgets import * #QApplication, QMainWindow,  QWidget, QProgressDialog, QInputDialog, QColorDialog, QLineEdit, QFileDialog, QAction, QTextEdit, QToolTip, QStatusBar, QMenuBar, QMessageBox, QScrollBar
from pyqtgraph.graphicsItems.ScatterPlotItem import ScatterPlotItem
QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True) #enable highdpi scaling
QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True) #use highdpi icons  
import pyqtgraph as pg
import pyqtgraph.opengl as gl
from processTNT import TNTfile
from scipy import constants
from scipy.ndimage import zoom   #used for interpolation of images
import pydicom    #pydicom is used to import DICOM images  pydicom.UID
from pydicom.dataset import Dataset, FileDataset
try:
    import ImageList  #class to make an image list from a stack of image files, used on PhantomViewer
except:
    pass
try:
  from PIL import Image  #imports tif gif
except:
  pass     
try:
  import imageio
except:
   print ('Can not import imageio needed for animated GIFs, try pip install imageio')

class TNMRviewer(QMainWindow):
      def __init__(self, daughter=None, parent = None):
        '''Define image window to display and reconstruct raw MRI data'''    
        super(TNMRviewer, self).__init__(parent)
        #self.win.setWindowFlags(Qt.WindowStaysOnTopHint)
        self.resize(1000,800)
        if daughter!=None:
            self.daughter=daughter
        self.imv = pg.ImageView( view = pg.PlotItem(),discreteTimeLine=True)
        # Creating a grid layout, adding slider on the left for 4th dimension, image view widget on the right
        layout = QGridLayout()
        self.paramSlider = QSlider(self)
        self.paramSlider.setStyleSheet("background : lightgray;")
        self.paramSlider.setSingleStep(1)
        self.paramSlider.setTickInterval(1)
        self.paramSlider.setTickPosition(3)
        self.lblParamId=QLabel()
        self.lblParamId.setText('Param')
        self.lblParamIndex=QLabel()
        self.lblParamIndex.setText('0')
        widget = QWidget()
        widget.setLayout(layout)
 
        # adding label in the layout
        layout.addWidget(self.paramSlider, 0, 0,)
        layout.addWidget(self.lblParamId, 2, 0)
        layout.addWidget(self.lblParamIndex, 1, 0)
        layout.setRowMinimumHeight(0,500)
        # plot window goes on right side, spanning 3 rows
        layout.addWidget(self.imv, 0, 1, 0, 3)
 
        # setting this widget as central widget of the main window
        self.setCentralWidget(widget)
        self.paramSlider.valueChanged.connect(self.paramSliderChange)
        self.imv.sigTimeChanged.connect(self.sliceSliderChange)
        self.imv.timeLine.setPen(width=10, color='lightblue')
 
        #self.setCentralWidget(self.imv)
        self.programID='TNMRviewer: '
        self.setWindowTitle(self.programID)
        point = self.rect().topRight()
        self.move(point + QPoint(int(self.width()/2), 0)) 
        self.statusBar = QStatusBar()
        self.setStatusBar(self.statusBar)
        self.menu = self.menuBar()
        self.view3DColor = QColor(255, 255 ,255 , alpha=10)
        self.view3DBackground = QColor(155, 155 ,255 , alpha=10)
        self.view3DTransparency = 2   #set transparency scaling for 3dview, 1 = transparency set by voxel value
        self.view3Dinvert = False    #flag to invert contrast in 3d image
        self.InstitutionName='NIST MIG'
        self.FoVX=160       #field of view in mm
        self.FoVY=160
        self.xscale=1       #=1 if image displays voxels, =FoV/#voxels if distance display, note xscale is in mm/voxel
        self.yscale=1 
        self.pIndex=2       #Set phase dimension (usually slice=0, readout=1, phase=2, parameter=3)
        self.roIndex=1      #Set readout dimension (usually slice=0, readout=1, phase=2, parameter=3)
        self.paramIndex=0
        self.sliceIndex=0
        self.imageOrientation=('Slice','Readout','Phase','Parameter')   #tuple to indicate input image orientation: default  first index = slice, second=Readout, third=Phase, fourth=Parameter eg TI, TE, b-value
        self.Or2Index={'Slice':0, 'Readout':1, 'Phase':2, 'Parameter':3}  #dictionary relating axes to to index
        self.Index2Or={0:'Slice', 1:'Readout', 2:'Phase', 3:'Parameter'}  #dictionary relating axes to to index
        self.plotOrientation=(0,1,2,3)  #tuple to determine plot dimension
        self.hiddenIndex='{}={}'.format('Parameter',0) 
        self.TE0=0.006501  # 1/2 minimum TE time in (s) TE=2(teDelay+TE0) 
        self.GEMSTE0=0.0072  #  minimum GEMS TE time in (s) TE=teDelay+GEMSTE0 for standard GEMS it is 7.2ms 
        self.TI0=0.004505  #  minimum TI time in (s) TI=tiDelay+TI0
        self.messages=''#messAages generated from image processing
        self.messageTextBox=self.InfoWindow()     #generate a message textbox
        self.messageTextBox.setWindowTitle('Messages')
        self.HorizontalLabel=''
        self.VerticalLabel=''
        self.HorizontalUnits='voxel index'
        self.VerticalUnits='voxel index'
        self.fileName=''

            #Physical Constants
        self.GammaPMHzperT=constants.physical_constants["proton gyromag. ratio over 2 pi"][0]   #in MHz/T 
        self.GammaPHzperT=1E6*constants.physical_constants["proton gyromag. ratio over 2 pi"][0]   #in Hz/T 
        self.GammaWaterProtonRadperT=constants.physical_constants["proton gyromag. ratio"][0]*(1-constants.physical_constants["proton mag. shielding correction"][0])                                                                                
        self.Gamma=self.GammaWaterProtonRadperT
        self.Gammaf=self.Gamma/2/np.pi
        
        #Gradient calibrations
        self.GxCalMRI=1.0E-3     #T/m/A
        self.GyCalMRI=1.0E-3     #T/m/A
        self.GzCalMRI=1.0E-3     #mT/m/A
        self.GxCal=self.GxCalMRI
        self.GyCal=self.GyCalMRI
        self.GzCal=self.GzCalMRI

        self.AxmaxMRI=85.0     #Maximum X gradient amp current when TNMR Amplitude is 100
        self.AymaxMRI=85.0    #Maximum Y gradient amp current when TNMR Amplitude is 100
        self.AzmaxMRI=85.0    #Maximum Z gradient amp current when TNMR Amplitude is 100
        self.Axmax=self.AxmaxMRI     #Maximum X gradient amp current when TNMR Amplitude is 100
        self.Aymax=self.AymaxMRI    #Maximum Y gradient amp current when TNMR Amplitude is 100
        self.Azmax=self.AzmaxMRI    #Maximum Z gradient amp current when TNMR Amplitude is 100

        self.GrCal=self.GxCal   #gradient calibration in slice direction
        self.GpCal=self.GyCal
        self.GsCal=self.GzCal
        self.Armax=self.Axmax    #Maximum X gradient amp current when TNMR Amplitude is 100
        self.Apmax=self.Aymax    #Maximum Y gradient amp current when TNMR Amplitude is 100
        self.Asmax=self.Azmax    #Maximum Z gradient amp current when TNMR Amplitude is 100 
        self.Gsmax=self.GsCal*self.Asmax        #maximum gradient in slice direction
        self.Grmax=self.GrCal*self.Armax        #maximum gradient in slice direction
        self.Gpmax=self.GpCal*self.Apmax        #maximum gradient in slice direction
        
        self.Gs=19      #Dac multiplier
        self.Gslice=self.Gs*self.Gsmax/100      

        self.fileMenu = self.menu.addMenu('&File')
        self.actionOpenTNTFile = QAction('Open .tnt file', self)
        self.actionOpenTNTFile.setStatusTip('Opens Tecmag file, assumed to be raw spin echo data')
        self.fileMenu.addAction(self.actionOpenTNTFile)
        self.actionOpenTNTFile.triggered.connect(self.openTNTFile)
        
        self.actionShowTNTComment = QAction('Show file comment', self)
        self.actionShowTNTComment.setStatusTip('Shows comment in tnt file and allows modification')
        self.fileMenu.addAction(self.actionShowTNTComment)
        self.actionShowTNTComment.triggered.connect(self.showTNTComment)
        
        self.actionShowTNTHeader = QAction('Show file header', self)
        self.actionShowTNTHeader.setStatusTip('Shows tnt file header and allows modification')
        self.fileMenu.addAction(self.actionShowTNTHeader)
        self.actionShowTNTHeader.triggered.connect(self.printFullTNTHeader)
        
        self.actionShowMessages = QAction('Show messages', self)
        self.actionShowMessages.setStatusTip('Shows messages')
        self.fileMenu.addAction(self.actionShowMessages)
        self.actionShowMessages.triggered.connect(self.showMessages)   
        
        self.actionWriteDICOM = QAction('Save as DICOM', self)
        self.fileMenu.addAction(self.actionWriteDICOM)
        self.actionWriteDICOM.triggered.connect(self.writeDicomFiles)     
        
        self.actionWriteAnimatedGIF = QAction('Save as animated GIF', self)
        self.fileMenu.addAction(self.actionWriteAnimatedGIF)
        self.actionWriteAnimatedGIF.triggered.connect(self.writeAnimatedGIF)     

        if daughter!=None:        
            self.actionMakeDaughter = QAction('Show Daughter Image Viewer', self)
            self.fileMenu.addAction(self.actionMakeDaughter)
            self.actionMakeDaughter.triggered.connect(self.makeDaughter) 
        
        self.imageMenu = self.menu.addMenu('&Images')
 
        self.actionplotMag = QAction('Plot Raw Magnitude', self)
        self.actionplotMag.setStatusTip('Plot raw magnitude images')
        self.imageMenu.addAction(self.actionplotMag)
        self.actionplotMag.triggered.connect(lambda: self.plotMag(sliceindex=self.sliceIndex))
        self.actionplotPhase = QAction('Plot Raw Phase', self)
        self.actionplotPhase.setStatusTip('Plot raw phase images')
        self.imageMenu.addAction(self.actionplotPhase)
        self.actionplotPhase.triggered.connect(lambda: self.plotPhase(sliceindex=self.sliceIndex))
        self.actionplotFFTMag = QAction('Plot Reconstructed Magnitude', self)
        self.actionplotFFTMag.setStatusTip('Plot FFT magnitude images')
        self.imageMenu.addAction(self.actionplotFFTMag)
        self.actionplotFFTMag.triggered.connect(lambda: self.plotFFTMag(sliceindex=self.sliceIndex))
        self.actionplotFFTPhase = QAction('Plot Reconstructed Phase', self)
        self.actionplotFFTPhase.setStatusTip('Plot FFT phase images')
        self.imageMenu.addAction(self.actionplotFFTPhase)
        self.actionplotFFTPhase.triggered.connect(lambda: self.plotFFTPhase(sliceindex=self.sliceIndex))
         
        self.actionSetImageOrientation = QAction('Set image orientation', self)
        self.imageMenu.addAction(self.actionSetImageOrientation)
        self.actionSetImageOrientation.triggered.connect(self.setImageOrientation)
        
        self.actionFlipROandPhase = QAction('Rotate image 180deg', self)
        self.imageMenu.addAction(self.actionFlipROandPhase)
        self.actionFlipROandPhase.triggered.connect(self.flipROandPhase)


#         self.actionExportImage = QAction('Export Image', self)
#         self.imageMenu.addAction(self.actionExportImage)
#         self.actionExportImage.triggered.connect(self.exportImages)
        
        self.imageMenu = self.menu.addMenu('&Processing')
                
        self.actionZeroPad = QAction('ZeroPad to equalize dimensions in FFT Recon', self)
        self.imageMenu.addAction(self.actionZeroPad)
        self.actionZeroPad.triggered.connect(self.zeroPad)

        self.actionResizeImage = QAction('Resize Image', self)
        self.imageMenu.addAction(self.actionResizeImage)
        self.actionResizeImage.triggered.connect(self.interpolate)
        self.actionNormalizeImage = QAction('normalize Images', self)
        self.imageMenu.addAction(self.actionNormalizeImage)
        self.actionNormalizeImage.triggered.connect(self.normalizeImages)
        self.actionSetFoV = QAction('Set field of view (FoV)', self)
        self.imageMenu.addAction(self.actionSetFoV)
        self.actionSetFoV.triggered.connect(self.setFoV)
        self.actionToggleScaledView = QAction('Toggle axes from Voxel index <-> distance(mm)', self)
        self.imageMenu.addAction(self.actionToggleScaledView)
        self.actionToggleScaledView.triggered.connect(self.toggleScaledVoxelView)
        
        self.actionZeroNoise = QAction('Zero Noise', self)
        self.imageMenu.addAction(self.actionZeroNoise)
        self.actionZeroNoise.triggered.connect(self.zeroNoise)
        
        self.actionZeroDataBorder = QAction('Zero border of kspace data', self)
        self.imageMenu.addAction(self.actionZeroDataBorder)
        self.actionZeroDataBorder.triggered.connect(self.zeroDataBorder)
        
        self.actionZeroKspacePoint = QAction('Enable Zero kspace point on mouse double click', self)
        self.imageMenu.addAction(self.actionZeroKspacePoint)
        self.actionZeroKspacePoint.triggered.connect(self.zeroKspacePoint)
                
        self.imageMenu = self.menu.addMenu('&3D Images')
        self.action3DImage = QAction('Plot 3D Reconstructed Image', self)
        self.imageMenu.addAction(self.action3DImage)
        self.action3DImage.triggered.connect(self.view3d)
        self.action3DImageInv = QAction('Plot inverted 3D Reconstructed Image', self)
        self.imageMenu.addAction(self.action3DImageInv)
        self.action3DImageInv.triggered.connect(lambda: self.view3d(invert=True))
        self.action3DImageColor = QAction('Set 3D Plot Color', self)
        self.imageMenu.addAction(self.action3DImageColor)
        self.action3DImageColor.triggered.connect(self.set3DColor)
        self.action3DImageTransparency = QAction('Set 3D Plot Transparency', self)
        self.imageMenu.addAction(self.action3DImageTransparency)
        self.action3DImageTransparency.triggered.connect(self.set3DTransparency)
        self.scaledData=False       #flag to indicate if scaled or raw data should be shown
        
        self.imageMenu = self.menu.addMenu('&NIST_7T_MRI')
        self.actionProcessScout = QAction('Process Scout', self)
        self.imageMenu.addAction(self.actionProcessScout)
        self.actionProcessScout.triggered.connect(self.processScout) 
        self.actionProcessGEFlash = QAction('Process GE_Flash or SE_T1', self)
        self.imageMenu.addAction(self.actionProcessGEFlash)
        self.actionProcessGEFlash.triggered.connect(self.processGEFlash)               
        self.actionExtractCh1 = QAction('Extract Ch1', self)
        self.imageMenu.addAction(self.actionExtractCh1)
        self.actionExtractCh1.triggered.connect(self.extractCh1)
        self.actionExtractCh2 = QAction('Extract Ch2', self)
        self.imageMenu.addAction(self.actionExtractCh2)
        self.actionExtractCh2.triggered.connect(self.extractCh2)
        self.actionExtractCh3 = QAction('Extract Ch3', self)
        self.imageMenu.addAction(self.actionExtractCh3)
        self.actionExtractCh3.triggered.connect(self.extractCh3)
        self.actionExtractCh4 = QAction('Extract Ch4', self)
        self.imageMenu.addAction(self.actionExtractCh4)
        self.actionExtractCh4.triggered.connect(self.extractCh4)
        self.actionSeparateImages = QAction('Separate Images', self)
        self.imageMenu.addAction(self.actionSeparateImages)
        self.actionSeparateImages.triggered.connect(self.separateImages)   
        
        self.imageMenu = self.menu.addMenu('&Measure')
        self.actionAddVerticalMarkers = QAction('Add vertical markers', self)
        self.imageMenu.addAction(self.actionAddVerticalMarkers)
        self.actionAddVerticalMarkers.triggered.connect(self.addVerticalMarkers)
        self.actionAddHorizontalMarkers = QAction('Add Horizontal markers', self)
        self.imageMenu.addAction(self.actionAddHorizontalMarkers)
        self.actionAddHorizontalMarkers.triggered.connect(self.addHorizontalMarkers)  
        self.actionRemoveMarkers = QAction('Remove markers', self)
        self.imageMenu.addAction(self.actionRemoveMarkers)
        self.actionRemoveMarkers.triggered.connect(self.removeMarkers)   
        
        self.imageMenu = self.menu.addMenu('&Plot')
        self.actionToggleHistogramLogMode = QAction('Toggle Histogram Log Mode', self)
        self.imageMenu.addAction(self.actionToggleHistogramLogMode)
        self.actionToggleHistogramLogMode.triggered.connect(self.toggleHistogramLogMode)
             
        
        self.imv.getView().setLabel('bottom',"H","voxel index")
        self.imv.getView().setLabel('left',"V","voxel index")
        self.dataAxes={'t':0, 'x':1, 'y':2}
        
        self.p1=pg.mkPen('r', width=1)
        self.p2=pg.mkPen('g', width=1)
        self.p3=pg.mkPen('c', width=1)
        self.p4=pg.mkPen('y', width=1)
        self.infv1 = pg.InfiniteLine(movable=True, angle=90, pen=self.p1,label='x={value:0.2f}', 
                       labelOpts={'position':0.1, 'color': (200,200,100), 'fill': (200,200,200,50), 'movable': True})
        self.infv2 = pg.InfiniteLine(movable=True, angle=90, pen=self.p1, label='x={value:0.2f}', 
                       labelOpts={'position':.1, 'color': (200,100,100), 'fill': (200,200,200,50), 'movable': True})
        self.infh1 = pg.InfiniteLine(movable=True, angle=0, pen=self.p2,  hoverPen=(0,200,0), label='y={value:0.2f}', 
                   labelOpts={'color': (200,0,0), 'movable': True, 'fill': (0, 0, 200, 100)})
        self.infh2 = pg.InfiniteLine(movable=True, angle=0, pen=self.p2,  hoverPen=(0,200,0), label='y={value:0.2f}', 
                   labelOpts={'color': (200,0,0), 'movable': True, 'fill': (0, 0, 200, 100)})
        self.proxy2 = pg.SignalProxy(self.imv.view.scene().sigMouseMoved, rateLimit=60, slot=self.mouseMove)
        #self.proxy2 = pg.SignalProxy(self.imv.view.scene().sigMouseMoved, rateLimit=60, slot=self.mouseMove)
        self.proxy3 = pg.SignalProxy(self.imv.view.scene().sigMouseClicked, rateLimit=60, slot=self.mouseButtonPress)
        self.zeroKspacePointonMouseClick=False  #flag to zero kspace point on left mouse click
        self.markerdH=0.0  #distance between horizontal markers
        self.markerdV=0.0   #distance between horizontal markers
        self.lineIntPoint=(0,0)     #initial position of measurement line, second point set by mouse double click
        self.measureLine = pg.LineROI((0,0),(0,0),0, pen=self.p4)
        self.measureLineLength=0        #measurement line length

        self.infh1.sigPositionChanged.connect(self.updateHMarkerLabel)
        self.infh2.sigPositionChanged.connect(self.updateHMarkerLabel)
        self.infv1.sigPositionChanged.connect(self.updateVMarkerLabel)
        self.infv2.sigPositionChanged.connect(self.updateVMarkerLabel)

      def clear(self):
          self.imv.clear()
          
      def openTNTFile(self):
        '''Open .tnt file, data is extracted as 4d array with readout,slice,phase,parameter 
        reorder to slice, readout, phase, parameter'''
        f = QFileDialog.getOpenFileName(self,'Open .tnt file', '', "Tecmag Files (*.tnt)")        #Note self.FileName is a qString not a string
        if f[0]=='':
          return 'cancel'
        self.fileName=f[0]
        self.tntfile=TNTfile(self.fileName)   #open tecmag tnt file
        self.message('filename=' +self.fileName)
        self.ProtocolName=""
        if self.fileName.find('SEMS')!= -1:
            self.ProtocolName="SEMS"    
        if self.fileName.find('SEMS_IR')!= -1:
            self.ProtocolName="SEMS_IR"  
        if self.fileName.find('GE_FLASH')!= -1:
            self.ProtocolName="GE_FLASH"    
        if self.fileName.find('GEMS')!= -1:
            self.ProtocolName="GEMS"   
        if self.fileName.find('GEMS_IR')!= -1:
            self.ProtocolName="GEMS_IR"    
        if self.fileName.find('SCOUT')!= -1:
            self.ProtocolName="SCOUT"       
        if self.fileName.find('PGSE_Dif')!= -1:
            self.ProtocolName="PGSE_Dif"       
        self.message('ProtocolName=' +self.ProtocolName)
        self.tntData=np.swapaxes(self.tntfile.DATA, 0,1) #input 4 dimensional data array, rearrange to get slice, readout, phase, parameter
        self.nSlice= self.tntData.shape[0]
        self.nReadout= self.tntData.shape[1]
        self.nPhase= self.tntData.shape[2]
        self.nParameter= self.tntData.shape[3]
        self.tfid=self.tntfile.fid_times()      #time points for FID waveforms
        self.freq=self.tntfile.freq_Hz()
        self.aquireTime=self.tntfile.spec_acq_time()
        self.ImagingFrequency=self.tntfile.ob_freq[0]        #frequency points for spectra
        self.MagneticFieldStrength=1E6*self.ImagingFrequency/self.Gammaf
        self.message('Magnetic field strength (T)={:.6f}'.format(self.MagneticFieldStrength))
        self.tntComment=self.tntfile.TNMRComment.split(';')
        self.tntStartTime=self.tntfile.start_time.isoformat()
        self.tntFinishTime=self.tntfile.finish_time.isoformat()
        self.StudyDate=self.tntFinishTime
        self.DataType='?' 
        self.Manufacturer='NIST Tecmag' 
        self.SeriesDescription='SE' 
        self.InstitutionName='NIST'          
        self.PatientName='Phantom'    
        self.bValue=0
        self.PixelBandwidth=40000      
        self.PixelSpacing=(0,1)    #column spacing, distance between neighboring columns 
        self.ImageType="MRI" 
        self.ReceiveCoilName='Agilent3T' 
                        #self.RowDirection.append(np.asfarray(ImageFile.ImageOrientationPatient[:3])) if hasattr(ImageFile,"ImageOrientationPatient") else self.RowDirection.append(np.array([1.,0.,0.]))
        self.RepetitionTime=5 
        self.EchoTime=.01 
        self.FlipAngle=90 
        #InPlanePhaseEncodingDirection) 
        self.InversionTime=0.1 
        self.SliceThickness=3 
        self.SliceLocation=0 
        self.tntGradOrientation=self.tntfile.gradOrientation.replace("'", '').replace('b','')   
        self.FoVX=160 
        self.FoVY=160
        self.HorizontalLabel='Readout {}'.format(self.tntGradOrientation[0])
        self.VerticalLabel='Phase {}'.format(self.tntGradOrientation[1])
        self.addPlotData(self.tntData,self.plotOrientation)
        self.tntTables=self.tntfile.DELAY
        try:
            self.sliceFrequencies=self.tntTables.get('sliceFrequencies')
            self.slicePositions=self.sliceFrequencies/(self.Gammaf*self.Gslice)
            self.message('Slice frequencies found, sliceFrequencies(Hz)={}, Slice positions(mm)={}'.format(np.array2string(self.sliceFrequencies,max_line_width=None, precision=2),np.array2string(self.slicePositions*1000,max_line_width=None, precision=2)))
        except:
            self.message('Slice frequencies not found')
        if self.ProtocolName=="SEMS":
          self.teDelay=self.tntTables.get('teDelay')
          self.TEArray=2*(self.teDelay+self.TE0)
          self.parameterArray=self.TEArray
          self.message('Set image protocol to T2SE, TEArray(ms)={}.'.format(np.array2string(self.TEArray*1000,max_line_width=None, precision=2)))
        if self.ProtocolName=="GEMS" :
          self.teDelay=self.tntTables.get('teDelay')
          self.TEArray=(self.teDelay+self.GEMSTE0)      #only one teDelay in GEMS
          self.parameterArray=self.TEArray
          self.message('Set image protocol to GEMS T2*, TEArray(ms)={}.'.format(np.array2string(self.TEArray*1000,max_line_width=None, precision=2)))
        if self.ProtocolName=="SEMS_IR":
          self.tiDelay=self.tntTables.get('tiDelay')
          self.TIArray=self.tiDelay+self.TI0
          self.parameterArray=self.TIArray
          self.message('Set image protocol to SEMS_IR, TIArray(ms)={}.'.format(np.array2string(self.TIArray*1000,max_line_width=None, precision=2)))
        if self.ProtocolName=="SCOUT":
            self.processScout(nRF=-1)
            self.message('Scout image')
        if self.ProtocolName=="PGSE_Dif":
          for com in self.tntComment:
              if com.find('bValues') != -1:
                  bvalues=com.split('=')[1]
                  self.bValueArray=np.fromstring(bvalues, sep=',')
                  self.parameterArray=self.bValueArray
                  self.message('Image protocol PGSE_Dif, b-ValueArray(s/mm^2)={}.'.format(np.array2string(self.bValueArray,max_line_width=None, precision=2)))
        #self.tntfile.TMAG[name]
        
      def flipROandPhase(self):
        self.tntData=np.flip(self.tntData, axis=1)
        self.tntData=np.flip(self.tntData, axis=2)
        self.addPlotData(self.tntData,self.plotOrientation)  
                
      def setImageOrientationPar_RO_Ph_Sl(self):
        self.plotOrientation=(3,1,2,0)
        self.addPlotData(self.tntData,self.plotOrientation)
      
      def setImageOrientation(self):
        orient=list(self.imageOrientation)
        slider, ok = QInputDialog.getItem(self, "Image orientation", "select slider axis", orient, 0, False)
        orient.remove(slider)
        xaxis, ok = QInputDialog.getItem(self, "Image orientation", "select x axis", orient, 0, False)
        orient.remove(xaxis)
        yaxis, ok = QInputDialog.getItem(self, "Image orientation", "select y axis", orient, 0, False)
        orient.remove(yaxis)
        hidden=orient[0]
        hdim=self.Or2Index[hidden]
        hmax=self.tntData.shape[hdim]
        hiddenIndex, ok = QInputDialog.getInt(self, "Image orientation", "input hidden index: "+hidden, 0, 0, hmax)
        self.hiddenIndex='{}={}'.format(hidden,hiddenIndex)
        #self.sliderIndex='{}={}'.format(slider,hiddenIndex)
        self.plotOrientation=(self.Or2Index[slider],self.Or2Index[xaxis],self.Or2Index[yaxis],self.Or2Index[hidden])

        self.addPlotData(self.tntData,self.plotOrientation)
          
      def setFoV(self):
        '''Sets field of view and e-display ReconMag image with new FoVs'''
        fov, OK= QInputDialog.getDouble(self,"Input field of view in X-direction", "FoVX(mm)", value=self.FoVX, decimals=1)
        if OK:
            self.FoVX=fov
        fov, OK= QInputDialog.getDouble(self,"Input field of view in Y-direction", "FoVY(mm)", value=self.FoVY, decimals=1)
        if OK:
            self.FoVY=fov
        self.xscale=self.FoVX/self.fftData.shape[1]
        self.yscale=self.FoVY/self.fftData.shape[2]
        self.HorizontalUnits='mm'
        self.VerticalUnits='mm'
        self.plotFFTMag()
                            

      def dimValue(self, d, ii):
          '''Returns value of array descriptor for dimension=d, index=ii'''
          ad=0.0
          units=''
          try:
              if d=='Slice':
                  ad=self.slicePositions[ii]*1000
                  units='mm'
              if d=='Parameter':
                  ad=self.parameterArray[ii]
          except:
            pass
          return ad, units
               
      def addVerticalMarkers(self): 
          self.imv.addItem(self.infv1)
          self.imv.addItem(self.infv2)
          
      def addHorizontalMarkers(self): 
          self.imv.addItem(self.infh1)
          self.imv.addItem(self.infh2)
          
      def removeMarkers(self): 
          self.imv.removeItem(self.infh1)
          self.imv.removeItem(self.infh2)
          self.imv.removeItem(self.infv1)
          self.imv.removeItem(self.infv2)
                    
      def updateVMarkerLabel(self):
        self.markerdV=self.infv2.value()-self.infv1.value()
      def updateHMarkerLabel(self):
        self.markerdH=self.infh2.value()-self.infh1.value()               
 
      def showTNTComment(self): 
          self.commentTextBox=self.InfoWindow()
          self.commentTextBox.setWindowTitle('tnt file Comment')
          self.commentTextBox.show()
          comment= "\n".join(str(element) for element in self.tntComment)
          self.commentTextBox.info.setText(comment)
                           
      def exportImages(self):
        self.imv.export("geeks")
        
      def toggleScaledVoxelView(self):
          if self.xscale==1 and self.yscale==1:
              self.xscale=self.FoVX/self.fftData.shape[1]
              self.yscale=self.FoVY/self.fftData.shape[2]
              self.imv.getView().setLabel('bottom',"X","mm")
              self.imv.getView().setLabel('left',"Y","mm")
          else:
              self.xscale=1
              self.yscale=1
              self.imv.getView().setLabel('bottom',"H","voxel index")
              self.imv.getView().setLabel('left',"V","voxel index")
          self.plotFFTMag()
        
      def plotMag(self, sliceindex=0, autolevel=True):
        self.setWindowTitle(self.programID+self.fileName+', RawImage Magnitude, shape={}'.format(self.rawData.shape))
        self.dataType='RawMag'        
        self.imv.getView().setLabel('bottom',self.HorizontalLabel,self.HorizontalUnits)
        self.imv.getView().setLabel('left',self.VerticalLabel,self.VerticalUnits)
        self.imv.setImage(np.absolute(self.rawData[:,:,:,self.paramIndex]),axes=self.dataAxes, autoLevels=autolevel)
        if sliceindex>=0 :
            self.imv.setCurrentIndex(sliceindex)
            self.imv.updateImage()
      def plotPhase(self, sliceindex=0, autolevel=True):
        self.setWindowTitle(self.programID+self.fileName+', RawImage Phase, shape={}'.format(self.rawData.shape))
        self.dataType='RawPhase'       
        self.imv.getView().setLabel('bottom',self.HorizontalLabel,self.HorizontalUnits)
        self.imv.getView().setLabel('left',self.VerticalLabel,self.VerticalUnits)
        self.imv.setImage(np.angle(self.rawData[:,:,:,self.paramIndex]),axes=self.dataAxes,autoLevels=autolevel)
        if sliceindex >=0:
            self.imv.setCurrentIndex(sliceindex)
            self.imv.updateImage()
      def plotFFTMag(self, sliceindex=0, autolevel=True):
        self.setWindowTitle(self.programID+self.fileName+', ReconstructedImage Magnitude, shape={}'.format(self.fftData.shape))
        self.dataType='ReconMag'
        self.imv.getView().setLabel('bottom',self.HorizontalLabel,self.HorizontalUnits)
        self.imv.getView().setLabel('left',self.VerticalLabel,self.VerticalUnits)
        self.imv.setImage(np.absolute(self.fftData[:,:,:,self.paramIndex]),axes=self.dataAxes,scale = (self.xscale,self.yscale),autoLevels=autolevel)
        if sliceindex >=0:
            self.imv.setCurrentIndex(sliceindex)
            self.imv.updateImage()
      def plotFFTPhase(self, sliceindex=0, autolevel=True):
        self.setWindowTitle(self.programID+self.fileName+', ReconstructedImage Phase, shape={}'.format(self.fftData.shape))
        self.dataType='ReconPhase'
        self.imv.getView().setLabel('bottom',self.HorizontalLabel,self.HorizontalUnits)
        self.imv.getView().setLabel('left',self.VerticalLabel,self.VerticalUnits)
        self.imv.setImage(np.angle(self.fftData[:,:,:,self.paramIndex]),axes=self.dataAxes,scale = (self.xscale,self.yscale),autoLevels=autolevel)
        if sliceindex >=0:
            self.imv.setCurrentIndex(sliceindex)
            self.imv.updateImage()

      def replotData(self):
        try:   
            if self.dataType=='RawMag':
                  self.plotMag(sliceindex=self.sliceIndex,autolevel=False)
            if self.dataType=='RawPhase':
                  self.plotPhase(sliceindex=self.sliceIndex,autolevel=False)
            if self.dataType=='ReconMag':
                  self.plotFFTMag(sliceindex=self.sliceIndex,autolevel=False)
            if self.dataType=='ReconPhase':
                  self.plotFFTPhase(sliceindex=self.sliceIndex,autolevel=False) 
        except:
            pass
                                    
      def addPlotData(self, data, imageorient=(0,1,2,3), nhidden=0, shape=0):
            '''inputs self.tntdata (4d complex data) being sent in from TNMR or .tnt file; self.rawData is 4d data sent to the plot widget with the first dimension
             on the bottom slider, second and third dimension on plot horizontal and vertical axes, 4th dimension on the left slider, 
            an FFT is applied to generate self.fftData with optional zero padding to obtain a square image of dim=shape,
            the FFT magnitude is then plotted'''

            self.rawData=np.transpose(data,imageorient)[:,:,:,:]     #4d data for the plot widget, may be transposed but usually it is [slice,RO,Phase,Param]
            self.paramSlider.setMaximum(data.shape[3]-1)  #left slider is for the 4th dimension          
            dat=self.rawData
            if shape!=0:        #if desired pad zeros equally on all sides to make square reconstructed image with dimension=shape
                n1=int(shape -dat.shape[1])
                n2=int(shape -dat.shape[2])
                if n1>0:
                    dat=np.pad(dat,((0,0),(int(n1/2),n1-int(n1/2)),(0,0)))
                if n2>0:
                    dat=np.pad(dat,((0,0),(0,0),(int(n2/2),n2-int(n2/2)),(0,0)))
            dat=np.fft.fftshift(dat,axes=(self.roIndex,self.pIndex))   #shift to make k=0 at upper left
            self.fftData=np.fft.fftshift(np.fft.fft2(dat,axes=(self.roIndex,self.pIndex)),axes=(self.roIndex,self.pIndex))     #FFT then shift o center image
            self.plotFFTMag()
 
      def fft2dRawData(self,data, index1=1, index2=2):
            '''FFTs raw data along indices given, assumes k=0 is in the center'''
            dat=np.fft.fftshift(data,axes=(index1,index2))   #shift to make k=0 at upper left
            self.fftData=np.fft.fftshift(np.fft.fft2(dat,axes=(index1,index2)),axes=(index1,index2))   
                     
      def gaussianFilter(self, data):
            self.rawData=data
            self.fftData=np.fft.fftshift(np.fft.fft2(data),axes=(1,2))
            self.setWindowTitle(self.programID+'Reconstructed Image Magnitude')
            self.imv.setImage(np.absolute(self.fftData),axes=self.dataAxes)
            
      def set3DColor(self):  
            self.view3DColor = QColorDialog.getColor()
            self.view3d()
 
      def set3DTransparency(self):  
        tr,ok = QInputDialog.getDouble(self, 'Input 3d transparency','1=transparent 10=not transparent', value=1.0, min=0.1, max=10,decimals=2)
        if not ok:
          return
        self.view3DTransparency=tr
        self.view3d() 
                                             
      def view3d(self, invert=False, transparency=2):
        '''creates 3d rendering of data currently in the image view'''  
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
        data=np.copy(self.imv.image)
        #scale image to levels set on the histogram
        self.min_level, self.max_level = self.imv.ui.histogram.getLevels()
        data-=self.min_level
        data[data< 0]=0.0
        data[data>(self.max_level-self.min_level)]=0.0
        data=data/np.amax(data)
        if invert:
            data=(1.0-data)
            data[data>1]=0
        d2 = np.empty(data.shape + (4,), dtype=np.ubyte)
        d2[..., 0] = data * self.view3DColor.red()
        d2[..., 1] = data * self.view3DColor.green()
        d2[..., 2] = data * self.view3DColor.blue()
        d2[..., 3] = (data)**self.view3DTransparency * 255.   #sets transparency  
        d2[:, 0:3, 0:3] = [255,0,0,20]   #draw axes at corner of box 
        d2[0:3, :, 0:3] = [0,255,0,20]
        d2[0:3, 0:3, :] = [0,0,255,20]    
        self.image3DVol=gl.GLVolumeItem(d2)
        self.image3DVol.translate(-data.shape[0]/2, -data.shape[1]/2,-data.shape[2]/2)
        self.view3Dwin.addItem(self.image3DVol)
        
      def interpolate(self):
        '''scales and interpolates 2d images using PIL.image.resize, changes image array size
        3d image are scaled using scipy.ndimage.zoom'''
        interp= ['bicubic', 'nearest', 'bilinear', 'lanczos']
        scalex,ok = QInputDialog.getDouble(self, "Scale X",'the number of rows will be increased (or decreased) by', value=2.0, min=0.1, max=10.0,decimals=2)
        if not ok:
          return
        scaley,ok = QInputDialog.getDouble(self, "Scale Y",'the number of columns will be increased (or decreased) by', 2.0, 0.1, 10.0,1)
        if not ok:
          return
        scalez,ok = QInputDialog.getDouble(self, "Scale Z",'the number of images in stack will be increased (or decreased) by', 2.0, 0.1, 10.0,1)
        if not ok:
            return
        self.scaledImage=zoom(np.absolute(self.fftData), (scalez,scalex,scaley,1),order=4)
        self.scaledPhase=zoom(np.angle(self.fftData), (scalez,scalex,scaley,1),order=4)
#        self.pw.message('Scaled image stack using scipy.ndimage.zoom, scalex={:.2f}, scaley={:.2f}, scalez={:.2f} ,interpolation=cubic spline image shape= {},{},{} \n'.format(scalex,scaley,scalez,s[0],s[1],s[2]))
        self.setWindowTitle(self.programID+'Reconstructed Image Magnitude, scaled by sx={}, sy={}, sz={}'.format(scalex, scaley,scalez))
        self.imv.setImage(np.absolute(self.scaledImage[:,:,:,0]),axes=self.dataAxes)
        self.scaledData=True
 
      def normalizeImages(self):
          max=np.amax(self.rawData, axis=(1,2))

      def extractCh1(self):
          self.addPlotData(self.tntData[:,0::4,:,:])
      def extractCh2(self):
          self.addPlotData(self.tntData[:,1::4,:,:])
      def extractCh3(self):
          self.addPlotData(self.tntData[:,2::4,:,:])
      def extractCh4(self):
          self.addPlotData(self.tntData[:,3::4,:,:])
                    
      def separateImages(self,n=3):
          '''Separates 3 consecutive images in echo, used for tri-scouts'''
          if n==0:
              n=3
          dshape=self.rawData.shape
          columns=int(dshape[1]/n)
          rows=dshape[2]
          self.rawData.shape=(3, columns, rows)
          self.addPlotData(self.rawData)
          
      def zeroPad(self,): 
          '''adds zeros to equalize image dimensions'''
          dshape=self.tntData.shape
          columns=dshape[1]  #
          rows=dshape[2]
          if rows< columns:
              npad=int((columns-rows)/2)
              self.tntData=np.pad(self.tntData, ((0,0), (0,0), (npad, npad),(0,0)), 'constant', constant_values=((0, 0),(0, 0),(0, 0),(0, 0)))
              self.message('Zero padded phase encode by {}, array dimension={}'.format(2*npad, str(self.tntData.shape)))
          self.addPlotData(self.tntData)
          
      def processScout(self, nRF=0): 
          '''Process Magnetica scout which has 3 views concatinated'''
          if nRF !=-1:
            self.tntData=self.tntData[:,nRF::4,:,:]     #extract channel 1 from  CH1,2,3,4 data for Magnetica data
          dshape=self.tntData.shape
          columns=int(dshape[1]/3)  #Saggital, coronal, axial need to be separated 
          rows=dshape[2]
          self.tntData=np.reshape(self.tntData,(3, columns, rows,dshape[3]))
          self.message('<b>msMRI Scout processing:</b> Extracted CH1, separated Coronal,Saggital,Axial slices')
          if rows< columns:
              npad=int((columns-rows)/2)
              self.tntData=np.pad(self.tntData, ((0,0), (0,0), (npad, npad),(0,0)), 'constant', constant_values=((0, 0),(0, 0),(0, 0),(0, 0)))
              self.message('Magnetica data: Zero padded phase encode by {}, array dimension={}'.format(2*npad, str(self.tntData.shape)))
          self.addPlotData(self.tntData)
          
      def processGEFlash(self, nRF=0):
          '''Process Magnetica style data which has 4 receive channels''' 
          if self.ProtocolName=="GE_FLASH" :
               nRF=3  
          self.tntData=self.tntData[:,nRF::4,:,:]     #extract channel nRF from  CH1,2,3,4 data 
          dshape=self.tntData.shape
          columns=dshape[1]
          rows=dshape[2]
          #self.parWin.message('<b>msMRI GE processing:</b> Extracted CH1')
          if rows< columns:
              npad=int((columns-rows)/2)
              self.tntData=np.pad(self.tntData, ((0,0), (0,0), (npad, npad),(0,0)), 'constant', constant_values=((0, 0),(0, 0),(0, 0),(0, 0)))
              self.message('Magnetica data: Zero padded phase encode by {}, array dimension={}'.format(2*npad, str(self.tntData.shape)))
          self.addPlotData(self.tntData)
          
      def zeroNoise(self):
          '''Converts all points that are outside the histogram levels to NaN
          These points will not be included in 3D plots or mapping''' 
          smin, smax = self.imv.ui.histogram.getLevels()
          self.fftData[np.absolute(self.fftData)< smin]=np.nan 
          self.fftData[np.absolute(self.fftData)> smax]=np.nan 
          self.plotFFTMag()
          
      def zeroKspacePoint(self):
          self.zeroKspacePointonMouseClick=True      
      
      def zeroDataBorder(self, d1=1,d2=2,n=10):
          n, ok = QInputDialog.getInt(self, "Zero k-space border", "input number of border voxels to zero: ", 20, 0, 100)
          if ok:
            self.tntData[:,0:n,:,:]=0
            self.tntData[:,:,0:n,:]=0
            self.tntData[:,-n:,:,:]=0
            self.tntData[:,:,-n:,:]=0
            self.addPlotData(self.tntData)
                    
      def printFullTNTHeader(self):
          """Expose members of the TMAG and TMG2 structures as attributes"""
          self.tntHeaderWin=self.InfoWindow(self)
          self.tntHeaderWin.setWindowTitle('TNT Header')
          self.tntHeaderWin.show()
          #self.tntHeaderWin.info.setText(str(self.fileName))
          for name in self.tntfile.TMAG.dtype.names:
               self.tntHeaderWin.info.append(name + '= '+ str(self.tntfile.TMAG[name])) 
          for name in self.tntfile.TMG2.dtype.names:
               self.tntHeaderWin.info.append(name + '= '+ str(self.tntfile.TMG2[name]))
          self.tntHeaderWin.info.append('TNT default tables; ' + str(self.tntfile.DELAY))
          
      def writeAnimatedGIF(self, filename='', duration=0.2):         
        filename = self.fileName.replace('.tnt', ".gif")
        im=np.copy(self.imv.image)
        im=np.swapaxes(im, 1,2)
        im-=im.min()    #set lowest value to 0
        im*=255.0/im.max()  #set highest value to 255
        imageio.mimsave(filename, im.astype(np.uint8), format='GIF', duration=duration)
                          
      def writeDicomFiles(self, filename):
            '''Writes current image stack to a stack of DICOM files'''
            f = QFileDialog.getSaveFileName(self,'Enter DICOM filename', '', "DICOM Files (*.dcm)")
            if f[0]=='':
                return 'cancel'
            if self.dataType=='ReconMag':
                dData=np.absolute(self.fftData[:,:,:,0])
            if self.dataType=='ReconPhase':
                dData=np.angle(self.fftData[:,:,:,0])   
            if self.dataType=='RawMag':
                dData=np.absolute(self.rawData[:,:,:,0])
            if self.dataType=='RawPhase':
                dData=np.angle(self.rawData[:,:,:,0])   
            if self.scaledData==True:
                 dData=self.scaledImage
            self.DICOMfileName=str(f[0])      #make sure fileName is a string, not a Qstring        
            nimages=dData.shape[0]
            ncolumns=dData.shape[1]
            nrows=dData.shape[2]
            for i in range(nimages):     #write out images in separate dicom files
                pixel_array=dData[i,:,:]
                fileName = self.DICOMfileName.replace('.dcm', str(i) + ".dcm")
                fileName = fileName.replace('/','\\')        
                # Populate required values for file meta information
                file_meta = Dataset()
                file_meta.MediaStorageSOPClassUID = b'Secondary Capture Image Storage'
                file_meta.MediaStorageSOPInstanceUID = b'1.3.6.1.4.1.9590.100.1.1.111165684411017669021768385720736873780'
                file_meta.ImplementationClassUID = b'1.3.6.1.4.1.9590.100.1.0.100.4.0'
                # Create the FileDataset instance (initially no data elements, but file_meta supplied)
                ds = FileDataset(fileName, {}, file_meta=file_meta, preamble=b"\0"*128)
                ds.Modality = b'MR'
                ds.ContentDate = str(datetime.date.today()).replace('-','')
#                 ds.ContentTime = str(time.time()) #milliseconds since the epoch
                ds.StudyInstanceUID =    b'1.3.6.1.4.1.9590.100.1.1.124313977412360175234271287472804872093'
                ds.SeriesInstanceUID = b'1.3.6.1.4.1.9590.100.1.1.369231118011061003403421859172643143649'
                ds.SOPInstanceUID =        b'1.3.6.1.4.1.9590.100.1.1.111165684411017669021768385720736873780'
                ds.SOPClassUID = b'MR Image Storage'
                ds.SecondaryCaptureDeviceManufctur = b'Python 2.7.3'
            ## These are the necessary imaging components of the FileDataset object.
                ds.SamplesPerPixel = 1
                ds.PhotometricInterpretation = b"MONOCHROME2"
                ds.PixelRepresentation = 0
                ds.HighBit = 15
                ds.BitsStored = 16
                ds.BitsAllocated = 16
                ds.SmallestImagePixelValue = b'\\x00\\x00'
                ds.LargestImagePixelValue = b'\\xff\\xff'
                #Add MRI data elements 
#                ds.bValue= self.bValue[i]
                ds.Columns=ncolumns
                ds.Rows= nrows                
                ds.PixelSpacing=[self.FoVX/ncolumns,self.FoVY/nrows]
#                ds.FoVX=self.FoVX
#                ds.FoVY=self.FoVY
#                 ds.ImageOrientationPatient=self.ImagePosition[i].tolist()
                ds.InstitutionName=self.InstitutionName
#                 ds.PixelBandwidth= self.PixelBandwidth[i] 
#                 ds.PixelSpacing =[self.PixelSpacingY[i],self.PixelSpacingX[i]]         

#                 ds.PatientName = self.PatientName[i]
#                 ds.PatientID = "123456"
#                 ds.ProtocolName=self.ProtocolName[i]
#                 ds.RepetitionTime = self.TR[i]
#                 ds.SeriesDescription=self.SeriesDescription[i]
#                 ds.EchoTime = self.TE[i]
#                 ds.FlipAngle= self.FA[i]
#                 ds.InversionTime=self.TI[i] 
#                 ds.ImageOrientationPatient[3:6]=self.ColumnDirection[i].tolist() 
#                 ds.ImageOrientationPatient[:3]=self.RowDirection[i].tolist()
#                 ds.MagneticFieldStrength = self.MagneticFieldStrength[i]
#                 ds.Manufacturer = self.Manufacturer[i]            
#                 pixel_array=np.transpose(self.PA[i])        #numpy    pixel array
                if self.dataType== "RawPhase" or self.dataType== "ReconPhase":
                    pixel_array=(pixel_array+np.pi) *10000    # phase data -pi to pi is converted to (0 to 2pi)*1000
                    pixel_array = pixel_array.astype(np.uint16)
                    #print "Adjusting phase to 16 bit integer"
                if self.dataType== "RawMag" or self.dataType== "ReconMag":
                    scale=2**15/np.amax(pixel_array)
                    pixel_array*=scale       #rescale so that max value is 16 bits                           
                    pixel_array = pixel_array.astype(np.uint16)
                #ds.PixelData = pixel_array.tostring()    # image byte data
                ds.PixelData = pixel_array.tobytes()    # image byte data
                #ds.PixelData = pixel_array    # image byte data
                # Set the transfer syntax
                ds.is_little_endian = True
                ds.is_implicit_VR = True
                ds.save_as(fileName, write_like_original=False)
                #pydicom.filewriter.write_file(fileName, ds, write_like_original=False)
                
      def makeImageStack(self): 
            """Makes a image stack for PhantomViewer from current plot, replicates a stack of DICOM images and headers"""
            self.stackIndex=3   
            sI,ok = QInputDialog.getInt(self, "Stack Index",'stack dimension 2=slice, 3=parameter', value=3, min=0, max=3)
            if ok:
                self.stackIndex=sI
            nslice=self.sliceIndex
            self.iS=ImageList.ImageList()
            print (self.ProtocolName)
            for i in range (self.rawData.shape[self.stackIndex]):
                #System and scan information
                self.iS.FileName.append(self.fileName) 
                self.iS.StudyDate.append(self.StudyDate) 
                self.iS.StudyDate.append(self.StudyDate) 
                self.iS.Manufacturer.append(self.Manufacturer) 
                self.iS.SeriesDescription.append(self.SeriesDescription) 
                self.iS.InstitutionName.append(self.InstitutionName)     
                self.iS.MagneticFieldStrength.append(self.MagneticFieldStrength) 
                self.iS.ImagingFrequency.append(self.ImagingFrequency)         
                self.iS.PatientName.append(self.PatientName)    
                self.iS.DataType.append(self.DataType)
                if self.stackIndex==3:  
                    self.iS.PA.append(np.absolute(self.fftData[nslice,:,:,i])) #add pixel arrays
                if self.stackIndex==0:  
                    self.iS.PA.append(np.absolute(self.fftData[i,:,:,0])) #add pixel arrays
                self.iS.Comment.append(self.tntComment) 

                self.iS.Columns.append(self.rawData.shape[1])
                self.iS.Rows.append(self.rawData.shape[2])
                self.iS.header.append('header') 
                self.iS.PixelBandwidth.append(self.PixelBandwidth)      
                self.iS.PixelSpacingX.append(self.FoVX/self.rawData.shape[1])    #column spacing, distance between neighboring columns
                self.iS.PixelSpacingY.append(self.FoVY/self.rawData.shape[2])   #row spacing, distance between neighboring rows  
                self.iS.ProtocolName.append(self.ProtocolName)    
                self.iS.ImagePosition.append(np.array([0,0,0]))
                self.iS.ImageType.append(self.ImageType) 
                self.iS.ReceiveCoilName.append(self.ReceiveCoilName) 
                self.iS.RowDirection.append(np.array([1.,0.,0.]))
                self.iS.ColumnDirection.append(np.array([0.,1.,0.]))
                self.iS.TR.append(self.RepetitionTime) 
                if (self.ProtocolName=='SEMS' or self.ProtocolName=='GEMS') and self.TEArray.shape[0]>0:
                    self.iS.TE.append(1000*self.TEArray[i])
                else:
                     self.iS.TE.append(self.EchoTime)
                if self.ProtocolName=='SEMS_IR' and self.stackIndex==3:
                    self.iS.TI.append(1000*self.TIArray[i])
                else:
                     self.iS.TI.append(self.InversionTime)
                
                if self.ProtocolName=='PGSE_Dif' and self.stackIndex==3:
                    self.iS.bValue.append(self.bValueArray[i])
                else:
                    self.iS.bValue.append(0.0)
                self.iS.FA.append(self.FlipAngle) 
                self.iS.InPlanePhaseEncodingDirection.append("")
                self.iS.SliceThickness.append(self.SliceThickness) 
                self.iS.SliceLocation.append(self.SliceLocation)
                self.iS.ScaleSlope.append (1.0)  #default value is 1
                self.iS.ScaleIntercept.append (0.0)    #default value is 0 
    
                self.iS.FoVX.append(self.FoVX) 
                self.iS.FoVY.append(self.FoVY)
            return  self.iS 
        
      def paramSliderChange(self):
          self.paramIndex=self.paramSlider.value()
          self.lblParamIndex.setText(str(self.paramIndex))
          self.replotData()
          
      def sliceSliderChange(self):
          self.sliceIndex=int(self.imv.currentIndex)         
       
      def toggleHistogramLogMode(self): 
          pass
          # self.imv.getHistogramWidget().item.axis.setLogMode(True,False )
          # self.imv.getHistogramWidget().item.vb.setLimits(yMin=1, yMax=16000)
             
      def makeDaughter(self):
            self.daughter.show()
          
      def message(self,m):   
            self.statusBar.showMessage(m)
            self.messages+= '\n' + m

      def showMessages(self):
          self.messageTextBox.show()
          self.messageTextBox.info.setText(self.messages)
                                  
      def mouseButtonPress(self, evt):
          '''Mouse click event button
          Used for zeroing kspace points or other processing options'''
          mouseEvent= evt[0]  ## using signal proxy turns original arguments into a tuple
          button=mouseEvent.button()
          pos=mouseEvent.pos()
          scenePos=mouseEvent.scenePos()
          #print(pos, scenePos, self.imv.view.vb.mapSceneToView(pos), self.imv.view.vb.mapToView(pos))
          if button==1 and self.dataType=='RawMag' and self.zeroKspacePointonMouseClick and mouseEvent.double():
                  if self.imv.view.sceneBoundingRect().contains(pos):
                    mousePoint = self.imv.view.vb.mapSceneToView(pos)
                    self.Xindex=int(mousePoint.x()) 
                    self.Yindex=int(mousePoint.y())
                    sliderindex=self.imv.timeIndex(self.imv.timeLine)[0]
                    if self.Xindex>=0 and self.Yindex>=0:
                        self.rawData[sliderindex,self.Xindex,self.Yindex,self.paramIndex]=0
                        self.statusBar.showMessage('zero:slice {}, X {}, Y {},Param {}'.format(sliderindex, self.Xindex,self.Yindex, self.paramIndex))
                        self.fft2dRawData(self.rawData, index1=1, index2=2)     #FFT then shift o center image
                        self.plotMag(sliceindex=sliderindex, autolevel=False)       #replot kspace data and re FFT
                        
          if button==1 and mouseEvent.double()and not self.zeroKspacePointonMouseClick:
                  if self.imv.view.sceneBoundingRect().contains(pos):
                    mousePoint = self.imv.view.vb.mapToView(pos)
                    self.Xindex=mousePoint.x() 
                    self.Yindex=mousePoint.y()
                    sliderindex=self.imv.timeIndex(self.imv.timeLine)[0]
                    self.imv.view.removeItem(self.measureLine)
                    if self.Xindex>=0 and self.Yindex>=0:
                        self.measureLine = pg.LineROI(self.lineIntPoint, (self.Xindex,self.Yindex), .1, pen=self.p4) 
                        self.imv.view.addItem(self.measureLine)
                        self.measureLineLength=np.sqrt((self.Xindex-self.lineIntPoint[0])**2+(self.Yindex-self.lineIntPoint[1])**2)
                        self.lineIntPoint=(self.Xindex,  self.Yindex)
                        
                        
          
      def mouseMove(self,evt): 
            '''mouse move event to display location and values'''
            pos = evt[0]  ## using signal proxy turns original arguments into a tuple
            im=self.imv.getImageItem().image
            if isinstance(im,np.ndarray):
                sliderindex=self.imv.timeIndex(self.imv.timeLine)[0]
                sliderDim=self.Index2Or[self.plotOrientation[0]]
                sliderArrayValue=self.dimValue(sliderDim, sliderindex)[0]
                sliderArrayDim=self.dimValue(sliderDim, sliderindex)[1]
                xDim=self.Index2Or[self.plotOrientation[1]]
                yDim=self.Index2Or[self.plotOrientation[2]]
                try:
                    if self.imv.view.sceneBoundingRect().contains(pos):
                        mousePoint = self.imv.view.vb.mapToView(pos)
                        self.Xindex=mousePoint.x() 
                        self.Yindex=mousePoint.y()
                        if self.Xindex>=0 and self.Xindex<im.shape[0] and self.Yindex>=0 and self.Yindex<im.shape[1]:
                            self.statusBar.showMessage('Slider_{}{}={:.2f}{}; X={}={:.2f}; Y={}={:.2f}; Value={:.2f}, Parameter={}, dH(mm)={:.2f}, dV(mm)={:.2f}, line length={:.2f}' \
                                .format(sliderDim, sliderindex, sliderArrayValue,sliderArrayDim, xDim, self.Xindex, yDim, self.Yindex, im[int(self.Xindex),int(self.Yindex)], self.paramIndex, self.markerdH, self.markerdV, np.sqrt((self.measureLine.size()[0]**2+self.measureLine.size()[1]**2))),5000)#
                except:
                    raise
                    pass
                    
              
      class InfoWindow(QMainWindow):
           def __init__(self, parent = None):
              '''Define info window/textbox, rv is the parent ROIView window'''    
              super(QMainWindow, self).__init__()
              self.resize(800,300)
              self.info = QTextEdit()
              self.setCentralWidget(self.info)
              self.setWindowTitle('Info')
              #self.pw=pw     
              
              
if __name__ == '__main__':
    app = QApplication(sys.argv)
    #clipboard=app.clipboard()
    #app.setStyleSheet("QWidget{font-size: 8pt;}") 
    tnmrViewerd = TNMRviewer()     #daughter windo in case need to view data in 2 windows
    tnmrViewer = TNMRviewer(daughter=tnmrViewerd)
    tnmrViewer.show()
    sys.exit(app.exec_())
                