'''
Created on April 29, 2017
Classes and routines to analyze system phantom fiducial array and compute geometric distortion
   
@author: stephen russek
updated 3-15-2020
'''

# -*- coding: utf-8 -*-
import sys
import numpy as np
import copy
from pyqt import *         #imports required PyQt modules
import pydicom    #import pydicom to read DICOM data sets and DICOMDIR
import pyqtgraph as pg    #uses pyqtgraph PlotWidget for graphs and ImageView windows for images, http://www.pyqtgraph.org/documentation/introduction.html
import pyqtgraph.functions as fn
import pyqtgraph.opengl as gl
import ImageList
import scipy.ndimage
import scipy.optimize
import skimage.measure
import matplotlib.pyplot as plt
import scipy.ndimage
import scipy.optimize
from scipy import signal
from scipy.spatial import procrustes
from Procrustes import mlprocrustes
import skimage.measure
import matplotlib.pyplot as plt
import GaussianFit
import SystemPhantom
import lmfit
                                 

class fiducialWindow():
  def __init__(self, pw=None,parent = None):
    '''Defines  window to view fiducial array and calculate geometric distortion'''    
    self.win = QMainWindow()
    #self.win.setWindowFlags(Qt.WindowStaysOnTopHint)
    self.win.resize(800,600)
    self.imv = pg.ImageView( view = pg.PlotItem())
    #self.imv.getView().invertY(False)
    #self.imv.ui.histogram.plot.setLogMode(None,True)    #set the histogram y axis to a log scale    
    self.imv.ui.roiBtn.setText("Line scan/ROI")
    self.win.setCentralWidget(self.imv)
    self.win.setWindowTitle('Fiducial Array')
    self.parentWindow=pw    #requires access to parent window to update ROI positions
    #point = self.rect().topRight()
    #self.win.move(point + QPoint(self.width()/2, 0))
    #cross hair
    self.vLine = pg.InfiniteLine(angle=90, movable=False)
    self.hLine = pg.InfiniteLine(angle=0, movable=False)
    self.imv.addItem(self.vLine, ignoreBounds=True)
    self.imv.addItem(self.hLine, ignoreBounds=True)
    self.label = pg.LabelItem(justify='top')
    self.label.setPos(0,0)
    self.imv.addItem(self.label)
    self.proxy = pg.SignalProxy(self.imv.view.scene().sigMouseMoved, rateLimit=60, slot=self.mouseMoved)
#    self.proxy2 = pg.SignalProxy(self.imv.view.scene().sigMouseClicked, rateLimit=60, slot=self.mouseClicked) 
    self.menu = self.win.menuBar()
    self.penr = fn.mkPen(255, 0, 0)
    self.peng = fn.mkPen(0, 255, 0)
    self.peny = fn.mkPen(255, 255, 0)
    self.penm = fn.mkPen(255, 0, 255)
    self.bblabelStyle = {'color':'w', 'font-size': '20px'}  #label and title styles for black and white background plots
    self.bbtitleStyle = {'color':'w', 'font-size': '20px'}
    self.wblabelStyle = {'color':'k', 'font-size': '20px'}
    self.wbtitleStyle = {'color':'k', 'font-size': '20px'}
    self.imv.getView().setLabel('bottom',"H","mm")
    self.imv.getView().setLabel('left',"V","mm")
    self.imv.sigTimeChanged.connect(self.indexChanged)
    self.messages=messageWindow()
    self.view3DColor = QColor(255, 255 ,255 , alpha=10)
    self.view3DBackground = QColor(155, 155 ,255 , alpha=10)
    self.view3DTransparency = 3   #set transparency scaling for 3dview, data=normlized_image**transparency
    self.view3Dinvert = False    #flag to invert contrast in 3d image
    self.scale3D=1.0    #sets overall brightness of 3d images
    self.viewPlane='Coronal'    #viewplane can be coronal sagittal or axial; analysis array is always x,y,z
    self.bViewSphereCenters=False
    self.pgROIs=[]    #list of pyqtgraph ROIs to be included on image
    self.bPlotROIs=True
    self.finalCropRadius=6.0    #final radius to crop spheres, determined voxels uses for center of mass
    self.sysPhantom=SystemPhantom.SystemPhantom()   #copy of system phantom for a reference
    #standard 3d image has 0,0,0 in upper left, as index changes x increases, y decreases, z decreases
    self.dntodx=np.array([1,-1,-1])   #required to convert index change to distance change

    self.menu = self.win.menuBar()
    
    self.imageMenu = self.menu.addMenu('&Images')
    self.actionViewImage = QAction('View phantom image', self.win)
    self.actionViewImage.triggered.connect(self.viewImage)
    self.imageMenu.addAction(self.actionViewImage) 
     
    self.actionViewCroppedImage = QAction('View cropped image', self.win)
    self.actionViewCroppedImage.triggered.connect(self.viewCroppedImage)
    self.imageMenu.addAction(self.actionViewCroppedImage)
    
    self.actionViewConvolvedImage = QAction('View convolved image', self.win)
    self.actionViewConvolvedImage.triggered.connect(self.viewConvolvedImage)
    self.imageMenu.addAction(self.actionViewConvolvedImage) 
 

    self.actionView3DImage = QAction('View 3d image', self.win)
    self.actionView3DImage.triggered.connect(self.view3DImage)
    self.imageMenu.addAction(self.actionView3DImage)
     
    self.actionView3DCroppedImage = QAction('View 3d cropped image', self.win)
    self.actionView3DCroppedImage.triggered.connect(self.view3DCroppedImage)
    self.imageMenu.addAction(self.actionView3DCroppedImage)
              
    self.actionviewSphereCenters = QAction('View sphere centers', self.win)
    self.actionviewSphereCenters.triggered.connect(self.viewSphereCenters)
    self.imageMenu.addAction(self.actionviewSphereCenters)
    
    self.actionviewMaskMag = QAction('View spherical mask', self.win)
    self.actionviewMaskMag.triggered.connect(self.viewMaskMag)
    self.imageMenu.addAction(self.actionviewMaskMag)
    
    self.actionviewkMaskMag = QAction('View spherical k-space mask', self.win)
    self.actionviewkMaskMag.triggered.connect(self.viewkMaskMag)
    self.imageMenu.addAction(self.actionviewkMaskMag)
    
    self.actionClear_All_Images = QAction('Clear All Images', self.win)
    self.actionClear_All_Images.setShortcut('Ctrl+C')
    self.actionClear_All_Images.setStatusTip('Clear All Images')
    #self.actionClear_All_Images.triggered.connect(clearImages)
    self.imageMenu.addAction(self.actionClear_All_Images)
    
    self.imageMenu = self.menu.addMenu('&View')
    self.actionViewAxial = QAction('Axial View', self.win)
    self.actionViewAxial.triggered.connect(self.viewAxial)
    self.imageMenu.addAction(self.actionViewAxial)
    self.actionViewCoronal = QAction('Coronal View', self.win)
    self.actionViewCoronal.triggered.connect(self.viewCoronal)
    self.imageMenu.addAction(self.actionViewCoronal)
    self.actionViewSagittal = QAction('Sagittal View', self.win)
    self.actionViewSagittal.triggered.connect(self.viewSagittal)
    self.imageMenu.addAction(self.actionViewSagittal)    
    self.actionViewMessages = QAction('View message/report window', self.win)
    self.actionViewMessages.triggered.connect(self.viewMessages)
    self.imageMenu.addAction(self.actionViewMessages)
    
    self.imageMenu = self.menu.addMenu('&Analysis')
    self.actionMakeMask = QAction('Make mask,convolve,find centers', self.win)
    self.actionMakeMask.triggered.connect(self.makeMaskandConvolve)
    self.imageMenu.addAction(self.actionMakeMask)
    self.actionRefineCenters = QAction('Refine centers', self.win)
    self.actionRefineCenters.triggered.connect(self.sphereRefinement)
    self.imageMenu.addAction(self.actionRefineCenters)
    self.actionCoordTransform = QAction('Transform Coordinates', self.win)
    self.actionCoordTransform.triggered.connect(self.coordTransform)
    self.imageMenu.addAction(self.actionCoordTransform)
    
    
    self.imageMenu = self.menu.addMenu('&ROIs')
    self.actionupdateROICenters = QAction('Update ROI centers from main window', self.win)
    self.actionupdateROICenters.triggered.connect(self.updateROICenters)
    self.imageMenu.addAction(self.actionupdateROICenters)
    self.actionPlotROIs = QAction('Show ROIs', self.win)
    self.actionPlotROIs.triggered.connect(self.plotROIs)
    self.imageMenu.addAction(self.actionPlotROIs)
    self.actionClearROIs = QAction('Clear ROIs', self.win)
    self.actionClearROIs.triggered.connect(self.clearROIs)
    self.imageMenu.addAction(self.actionClearROIs)
    self.actionSetROId = QAction('Set ROI diameter', self.win)
    self.actionSetROId.triggered.connect(self.setROId)
    self.imageMenu.addAction(self.actionSetROId)
    self.actionPrintCData = QAction('Save convolution data', self.win)
    self.actionPrintCData.triggered.connect(self.printCData)
    self.imageMenu.addAction(self.actionPrintCData)
    self.actionview3DColor = QAction('view3DColor', self.win)
    self.actionview3DColor.triggered.connect(self.v3DColor)
    self.imageMenu.addAction(self.actionview3DColor)
    self.action3DTransparency = QAction('3DTransparency', self.win)
    self.action3DTransparency.triggered.connect(self.v3DTransparency)
    self.imageMenu.addAction(self.action3DTransparency)
    self.action3DIntensity = QAction('3DIntensity', self.win)
    self.action3DIntensity.triggered.connect(self.v3DIntensity)
    self.imageMenu.addAction(self.action3DIntensity)      

  def viewMessages(self):
    '''show message text box which contains output results'''
    self.messages.win.show()

  def viewSphereCenters(self):
    '''toggle view of current sphere centers in cropped image view'''
    self.bViewSphereCenters= not self.bViewSphereCenters
    self.viewCroppedImage()
    
  def viewAxial(self):
    self.viewPlane='Axial'
    self.fovh=self.fovx
    self.fovv=self.fovy
    self.viewImage()
    self.plotROIs()
  def viewCoronal(self):
    self.viewPlane='Coronal'
    self.fovh=self.fovx
    self.fovv=self.fovz
    self.viewImage()
    self.plotROIs()
  def viewSagittal(self):
    self.viewPlane='Sagittal'
    self.fovh=self.fovz
    self.fovv=self.fovy
    self.viewImage()
    self.plotROIs()
                     
  def getDouble(self, text1='QInputDialog', text2='input floating point value', default=0.0, min=0, max=100, decimals=2):
      double,ok = QInputDialog.getDouble(self.win, text1,text2, default, min, max, decimals)
      if ok:
        return double
      else:
        return None

  def displayRealImage(self,image,title):
    '''display 3d image, stack/3d dimension needs to be first'''
    self.windowTitle=title
    self.win.setWindowTitle(title)
    if self.bViewSphereCenters:#add bright spots at array origin, array enter, and at sphere centers
      icenter=(np.asarray(image.shape)/2).astype('int')
      imax=np.amax(image)
      im=image+self.sphereCenterImage
      im[0,0,0]=2*imax
      im[icenter[0], icenter[1],icenter[2]]=2*imax
    else:
      im=image
    if self.viewPlane=='Coronal':#transform xyz to yxz (first index is depth, 2nd is horizontal, third is vertical)
      self.imv.setImage(np.transpose(im, (1,0,2)),xvals=self.zarray, pos=(-self.fovx/2,-self.fovz/2), scale = (self.dx,self.dz))
      self.imv.getView().setLabel('bottom','X(mm)')
      self.imv.getView().setLabel('left','-Z(mm)')
    if self.viewPlane=='Axial':#transform xyz to zxy
      #self.imv.setImage(np.transpose(im, (2,0,1)),xvals=self.zarray, pos=(-self.fovx/2,-self.fovy/2), scale = (self.dx,self.dy))
      self.imv.setImage(np.flip(np.transpose(im, (2,0,1)),2),xvals=self.zarray, pos=(-self.fovx/2,-self.fovy/2), scale = (self.dx,self.dy))
      self.imv.getView().setLabel('bottom','X',"mm")
      self.imv.getView().setLabel('left','Y',"mm")
    if self.viewPlane=='Sagittal':  #transform xyz to xzy
      self.imv.setImage(np.flip(np.flip(np.transpose(im, (0,2,1)),2),1),xvals=self.zarray, pos=(-self.fovz/2,-self.fovy/2), scale = (self.dz,self.dy))
      self.imv.getView().setLabel('bottom','Z',"mm")
      self.imv.getView().setLabel('left','Y',"mm")
    self.imv.setCurrentIndex(int(im.shape[0]/2))
    
  def viewImage(self): 
    self.displayRealImage(self.image,'Fiducial Array: Phantom Image'+self.header)
    self.currentArray=self.image
    
  def viewCroppedImage(self): 
    self.displayRealImage(self.croppedImage,'Fiducial Array: Cropped Image'+self.header)
    self.currentArray=self.croppedImage
    
  def viewConvolvedImage(self):   
    self.displayRealImage(self.convolveArray,'Convolution'+self.header)
    self.currentArray=self.convolveArray

  def view3DImage(self):
    self.view3d(self.image)         

  def viewMaskMag(self):
    self.displayRealImage(self.maskMag,'Fiducial Array: Mask Image'+self.header) 
    self.currentArray=self.maskMag    
  
  def viewkMaskMag(self):
    image=np.absolute(self.kImage)
    self.imv.setImage(image,xvals=self.zarray, pos=(-np.pi/self.dx,-np.pi/self.dx), scale = (np.pi/self.fovx,np.pi/self.fovx))
    self.imv.getView().setLabel('bottom','kx',"1/mm")
    self.imv.getView().setLabel('left','kz',"1/mm")
    self.imv.setCurrentIndex(int(image.shape[0]/2))
      
  def view3DCroppedImage(self):
    self.view3d(self.croppedImage)         

  def addImage(self,ims,cROIs, ny=0):
    '''passes in an image stack along with all metadata''' 
    self.ims=ims    #image stack
    self.filename=str(ims.FileName[1])
    self.fileDir=self.filename[0:self.filename.rfind("\\")]   
    self.fileDir=self.fileDir[self.fileDir.rfind("\\"):]  #use folder Name for titles and ID
    self.currentROIs=cROIs
    self.image=np.swapaxes(self.ims.np3dArray(),0,1)    #input stack has stack dimension first, in standard coronal this should be y axis 2    
    self.currentArray=self.image
    self.fovx=self.ims.FoVX[1]   #isotropic field of view hopefully
    self.fovz=self.ims.FoVY[1]   
    self.fovy=float(self.ims.SliceThickness[1]*self.image.shape[1])  #assume coronal view with y being the stacking direction
    self.fovh=self.fovx
    self.fovv=self.fovz
    self.dx=self.fovx/self.image.shape[0]          #x voxel size
    self.dy=self.fovy/self.image.shape[1]          #y voxel size
    self.dz=self.fovz/self.image.shape[2]          #z voxel size
    self.dr=np.array([self.dx,self.dy,self.dz])
    self.sphereCenterImage=np.zeros(self.image.shape)   #array with voxels set bright to show approx sphere centers
    self.nsize=self.image.shape[0]    #array size
    self.yZero=self.fovy/2.0-self.dy*ny #voxel corresponding 
    self.voxOrigin=np.array([-self.fovx/2,self.fovy/2,self.fovz/2])   #standard coronal image 0,0 corresponds to -x (R), +y(A), +z(S) corner of the phantom
    self.zarray=np.arange(self.nsize)*self.ims.SliceThickness[1]
    self.header=self.filename +'\n' + 'Image Shape: {}, FoVx,y,z(mm)={:.2f},{:.2f},{:.2f}, dx,dy,dz(mm)={:.3f},{:.3f},{:.3f}, \
        Voxel 0,0,0 corner pos(mm)={:.3f},{:.3f},{:.3f}\n'.format(
        str(self.image.shape),self.fovx,self.fovy,self.fovz,self.dx,self.dy,self.dz,self.voxOrigin[0],self.voxOrigin[1],self.voxOrigin[2])
    self.messages.message(self.header)    
    self.sphereCenters=np.zeros((len(self.currentROIs.ROIs),3))  #sphere centers, initialize to nominal ROI positions
    self.sphereRCenters=np.zeros(len(self.currentROIs.ROIs))  #distance of sphere centers from center of image, convenient parameter for graphs
    self.sphereVolumes=np.zeros(len(self.currentROIs.ROIs))  #sphere volumes relative to central sphere
    self.initialSphereCenters=np.zeros((len(self.currentROIs.ROIs),3))  #sphere centers prescribed by the phantom
    self.sphereIntensity=np.zeros(len(self.currentROIs.ROIs))  #Integrated intensity of each sphere
    self.dRc=np.zeros((len(self.currentROIs.ROIs),3))  #array of displacement vector between prescribed and actual sphere centers
    self.sphereRadius=np.zeros((len(self.currentROIs.ROIs),3))  #array of sphere/ellipsoid radii
    self.fiducialColors=[pg.mkColor(roi.color) for roi in self.currentROIs.ROIs]    #make a list of colors so ROIs can be distinguished
    self.fiducialSymbols=[roi.symbol for roi in self.currentROIs.ROIs]    #make a list of symbols
    for n, roi in enumerate(self.currentROIs.ROIs): #start by assuming center of masses are at nominal positions
      self.sphereCenters[n,0]=roi.Xcenter
      self.sphereCenters[n,1]=roi.Ycenter+self.yZero
      self.sphereCenters[n,2]=roi.Zcenter
      self.initialSphereCenters[n,0]=self.currentROIs.initialROIs[n].Xcenter
      self.initialSphereCenters[n,1]=self.currentROIs.initialROIs[n].Ycenter
      self.initialSphereCenters[n,2]=self.currentROIs.initialROIs[n].Zcenter
    self.croppedImage= self.extractFidSpheres(self.image,r=10)
    self.viewImage()
    self.plotROIs()


  def cropFiducials(self):
    radius=self.getDouble('Enter radius', 'mm', default=8, min=0, max=50, decimals=2)
    self.croppedImage=self.extractFidSpheres(self.image, r=radius)
    self.viewCroppedImage()
    return
   
  def extractFidSpheres(self, image, r=10):
    '''Extract fiducial spheres from 3d image, crop to ~radius r in mm'''
    nr=int(r/self.dx) #radius in voxels for extraction region around each sphere
    self.maxS=np.amax(image)
    ncenter=int(image.shape[0]/2)
    arr=np.zeros(self.image.shape)
    sum=np.zeros(len(self.currentROIs.ROIs))    #total signal for fiducial sphere volume
    com=np.zeros((len(self.currentROIs.ROIs),3))  #center of Mass
    #makes array arr by cropping image in a radius nr about each roi center
    for n, roi in enumerate(self.currentROIs.ROIs):
      nroi=self.vIndex(self.sphereCenters[n])
      for i in range(-nr,nr+1):
        for j in range (-nr,nr+1):
          for k in range (-nr,nr+1):
            if i**2+j**2+k**2<=r**2:
              nx=nroi[0]+j
              ny=nroi[1]+i
              nz=nroi[2]+k
              arr[nx,ny,nz]=self.image[nx,ny,nz]
              sum[n]+=arr[nx,ny,nz]
              com[n]+=arr[nx,ny,nz]*self.vCenter(np.array([nx,ny,nz]))
      #Make sphere center indicator image
      nin=self.vIndex(self.sphereCenters[n])
      self.sphereCenterImage[nin[0],nin[1],nin[2]]=2*self.maxS    
    self.sphereCoMs=com/sum.reshape(len(self.currentROIs.ROIs),1)  #center of mass is weighted sum over sphere
    return arr

  def vCenter(self, nr):
    '''returns voxel center in mm of voxel at nx, ny, nz'''
    r=self.voxOrigin+self.dr*(nr+0.5)*self.dntodx
    return r
  
  def vIndex(self, r):
    '''returns voxel index which contains point r in mm '''
    nr=np.absolute(r-self.voxOrigin)/self.dr
    return nr.astype(int)
  
  def rtoIndex(self, r):
    '''returns float voxel coord corresponding to r(mm) or integer part gives
    index which contains point r in mm '''
    nr=np.absolute(r-self.voxOrigin)/self.dr
    return nr
        

  def makeMaskandConvolve(self):  
    self.makeMask()
    self.convolveWithMask()
      
  def makeMask(self):
    '''Creates a spherical mask'''
    self.win.setWindowTitle('Fiducial Analysis: Making fiducial sphere mask ....')
    self.makeFTSphereImage(a=5)  #creates kspace and synthetic images 'self.synImage' 'self.kImage'


  def FTsphere(self,kx,ky,kz,a,x0,y0,z0, sigma=0.2, filterOn=True): 
    '''FT of sphere of radius a at x0,y0,z0; kx, ky, kz are 3d matrices'''
    k2=kx**2+ky**2+kz**2
    ka=np.sqrt(k2)*a
    ka[ka==0]=1E-6
    J32ka=(np.sin(ka)-ka*np.cos(ka))/ka**3 
    ftsphere=(4*np.pi*a**3)*np.exp(1j*x0*kx)*np.exp(1j*y0*ky)*np.exp(1j*z0*kz)*J32ka
#     Jka=scipy.special.jn(1.5,ka)
#     ftsphere=(2*np.pi*a**2/ka)**1.5*np.exp(1j*x0*kx)*np.exp(1j*y0*ky)*np.exp(1j*z0*kz)*Jka
    if filterOn:
      filter=np.exp(-sigma**2*k2)
      return ftsphere*filter
    else:
      return ftsphere

  def makeFTSphereImage(self, a=1, x0=0.0, y0=0.0, z0=0.0):
    '''Calculates discrete FT of sphere of radius a with center at x0,y0,z0, returns kspace image'''
    npoints=self.image.shape[0]  #assume all array dimension the same for now
    n=int(npoints/2)
    x=range(-n,n)
    y=range(-n,n)
    z=range(-n,n)
    kx, ky, kz = np.meshgrid(x,y,z)
    self.kx=2*np.pi*kx/self.fovx    #kx*dx goes from -pi to pi
    self.ky=2*np.pi*ky/self.fovy
    self.kz=2*np.pi*kz/self.fovz
    dx=self.dx/2    #translate 1/2 voxel since images are even and center is between 2 voxels
    dy=self.dy/2
    dz=self.dz/2
    self.kImage=self.FTsphere(self.kx,self.ky,self.kz,a,dx,dy,dz)
    self.maskComplex=np.fft.fftshift(np.fft.ifftn(self.kImage))    #3d fft then shift so zero frequency is in center
    self.maskMag=np.absolute(self.maskComplex)
    immax=np.amax(self.maskMag)
    self.maskMag[self.maskMag<immax/50]=0.0 #set small values to zero
    
  def convolveWithMask(self):
    self.win.setWindowTitle('Fiducial Analysis: Convolving image with mask ....')
    self.sphereCenterImage=np.zeros(self.image.shape)   #image showing voxel closest to sphere center
    a=int(self.maskMag.shape[0]/2-10)
    b=int(self.maskMag.shape[0]/2+10)
    self.sphereMask=self.maskMag[a:b,a:b,a:b]
    
    #pass 1 find approximate centers
    self.convolveArray = signal.fftconvolve(self.croppedImage, self.sphereMask, mode='same')
    #find local maximal voxel for each ROI
    d=8   #2d = size of neighborhood to around fiducial sphere to analyze
    for n, roi in enumerate(self.currentROIs.ROIs):   
      r=self.sphereCenters[n]
      nr=self.vIndex(r)
      neighborhood=self.convolveArray[nr[0]-d:nr[0]+d,nr[1]-d:nr[1]+d,nr[2]-d:nr[2]+d]
      indexofMax=np.asarray(np.unravel_index(np.argmax(neighborhood), neighborhood.shape)) #index of maximum in nb  
      nin=indexofMax+[nr[0]-d,nr[1]-d,nr[2]-d]
      self.sphereCenters[n]=self.vCenter(nin) #approximate sphere center as voxel ith largest correlation
    #pass 2 refine  centers      
    self.sphereCenterRefinement(finalradius=7, quiet=True, refinement=1)
    #pass 3 refine  centers      
    self.sphereCenterRefinement(finalradius=7, quiet=True,refinement=2)
    #pass 4 refine  centers and output     
    self.sphereCenterRefinement(finalradius=7, quiet=False,refinement=3)

  def sphereRefinement(self):
    self.sphereCenterRefinement(finalradius=7, quiet=False)
    
  def sphereCenterRefinement(self,finalradius=7, quiet=True, refinement=1):
    #refinement pass to find  centers, fits convolution to gaussian functions
    self.win.setWindowTitle('Fiducial Analysis: Convolving image with mask: refinement ....')
    d=10   #2d = size of neighborhood to around fiducial sphere to analyze
    self.croppedImage=self.extractFidSpheres(self.image, r=finalradius)  
    self.convolveArray = signal.fftconvolve(self.croppedImage, self.sphereMask, mode='same')
    self.fitArray=np.zeros((len(self.currentROIs.ROIs),2*d,6))    #saves raw data that is fit to gaussians
#    self.displayRealImage(self.convolveArray,'Convolution')
    for n, roi in enumerate(self.currentROIs.ROIs):   
      r=self.sphereCenters[n]
      nr=self.vIndex(r)
      neighborhood=self.convolveArray[nr[0]-d:nr[0]+d,nr[1]-d:nr[1]+d,nr[2]-d:nr[2]+d]
      roi.array=self.croppedImage[nr[0]-d:nr[0]+d,nr[1]-d:nr[1]+d,nr[2]-d:nr[2]+d]  #this array is used to determine local signal intesity
      indexofMax=np.asarray(np.unravel_index(np.argmax(neighborhood), neighborhood.shape)) #index of maximum in nb  
      nin=indexofMax+[nr[0]-d,nr[1]-d,nr[2]-d]
      self.sphereCenters[n]=self.vCenter(nin) #approximate sphere center as voxel with largest correlation
      self.sphereCenterImage[nin[0],nin[1],nin[2]]=2*self.maxS
      xsum=np.sum(neighborhood, axis=(1,2))
      ysum=np.sum(neighborhood, axis=(0,2))
      zsum=np.sum(neighborhood, axis=(0,1))      
      xi=self.voxOrigin[0]+self.dx*np.arange(nr[0]-d,nr[0]+d)  #+self.dx/2.0
      xc,sigmaX=self.findSphereCenters(xi,xsum)
      yi=self.voxOrigin[1]-self.dy*np.arange(nr[1]-d,nr[1]+d)  #-self.dy/2.0
      yc,sigmaY=self.findSphereCenters(yi,ysum)
#      print n , y, np.sum(neighborhood, axis=(0,2))/1000000, yc
      zi=self.voxOrigin[2]-self.dz*np.arange(nr[2]-d,nr[2]+d)  #-self.dz/2.0
      zc,sigmaZ=self.findSphereCenters(zi,zsum)
      self.sphereCenters[n]=[xc,yc,zc]
      self.sphereRadius[n]=[sigmaX,sigmaY,sigmaZ]
      self.dRc[n]=[xc-roi.Xcenter,yc-roi.Ycenter,zc-roi.Zcenter]
      self.fitArray[n,:,0]=xi
      self.fitArray[n,:,1]=xsum
      self.fitArray[n,:,2]=yi
      self.fitArray[n,:,3]=ysum
      self.fitArray[n,:,4]=zi
      self.fitArray[n,:,5]=zsum
    self.coordTransform(refinement=refinement)
    if quiet != True:
      self.showSpherePositions()
    self.clearROIs()
    self.plotROIs()
    self.croppedImage=self.extractFidSpheres(self.image, r=self.finalCropRadius)  
    self.win.setWindowTitle('Fiducial Analysis: Updated sphere centers, integrated signal, and ellipticities')

  def coordTransform(self, refinement=1):
    '''performs procrustean transformation to translate, rotate, stretch  ideal sphere centers self.initialSphereCenters
     to match measured sphere centers self.sphereCenters'''
    self.pTransformD, self.pTransformSphereCenterRef, self.pTransform = mlprocrustes(self.sphereCenters, self.initialSphereCenters, scaling=True, reflection=False)
    self.messages.message('Procustes transform refinement{}: Disparity={:.3e}, scale={:.4f}, translation='.format(refinement,self.pTransformD,self.pTransform['scale'])
      +np.array2string(self.pTransform['translation'], precision=4) + ', rotation= \n' + np.array2string(self.pTransform['rotation'], precision=6) + '\n')
    self.geoDistortion= self.sphereCenters-self.pTransformSphereCenterRef
  
  def printCData(self):
    fileName = QFileDialog.getSaveFileName(parent=None, caption="Report File Name", directory = '', selectedFilter = ".txt")
    if not fileName:  #if cancel is pressed return
      return None
    n=20
    fout=self.fitArray[n,:,:]
    np.savetxt(str(fileName), fout, delimiter= ',', fmt='%10.3f', header="ROI"+str(n))
      
  def showSpherePositions(self):
    '''print sphere positions in mm, voxel coordinates, center of mass, and deviation from prescribed positions
     after Procrustes transformation''' 
    self.messages.message('ROI:\t X(mm)\t Y(mm)\t Z(mm)\t nx\t ny\t nz\t sx(mm)\t sy(mm)\t sz(mm)\t CoMx\t CoMy\t CoMz\t Sum\t R0/R0n\t R(mm)\t dX(mm)\t dY(mm)\t dZ(mm) \n')
    for n, roi in enumerate(self.currentROIs.ROIs):  
      x=self.sphereCenters[n,0]
      y=self.sphereCenters[n,1]
      z=self.sphereCenters[n,2]
      Rs=(x**2+y**2+z**2)**0.5
      ijk=self.rtoIndex(self.sphereCenters[n])  #sphere center in index coordinates
      x0=self.sphereCenters[28,0]  #apparent position of the central sphere 28
      y0=self.sphereCenters[28,1]
      z0=self.sphereCenters[28,2]
#       dx=self.dRc[n,0]-x0
#       dy=self.dRc[n,1]-y0
#       dz=self.dRc[n,2]-z0
      dx=self.geoDistortion[n,0]
      dy=self.geoDistortion[n,1]
      dz=self.geoDistortion[n,2]
      sx=self.sphereRadius[n,0]
      sy=self.sphereRadius[n,1]
      sz=self.sphereRadius[n,2]
      self.CoMminusCovCenter=self.sphereCoMs-self.sphereCenters   #compare center of mass and convolution sphere centers
      self.sphereVolumes[n]=sx*sy*sz/(self.sphereRadius[28,0]*self.sphereRadius[28,1]*self.sphereRadius[28,2])
      asum=np.sum(roi.array)/1E6
      self.sphereIntensity[n]=asum
      Rn=roi.Rcenter()-self.currentROIs.ROIs[28].Rcenter()  #prescribed distance between centers
      self.sphereRCenters[n]=((x-x0)**2+(y-y0)**2+(z-z0)**2)**0.5 #calculate apparent sphere distance from center sphere
      RnLength=np.linalg.norm(Rn)
      if RnLength>0:  #calculate apparent distance from center sphere, normalize by prescribed distance
        dR0=self.sphereRCenters[n]/RnLength
      else:
        dR0=1.0
      self.messages.message('{}\t {:.3f}\t {:.3f}\t {:.3f}\t{:.3f}\t {:.3f}\t {:.3f}\t{:.3f}\t {:.3f}\t {:.3f}\t{:.3f}\t {:.3f}\t {:.3f}\t {:.4f}\t {:.4f}\t {:.3f}\t {:.3f}\t {:.3f}\t{:.3f}\n'
          .format(str(n+1),x,y,z,ijk[0],ijk[1],ijk[2],sx,sy,sz,self.sphereCoMs[n,0], self.sphereCoMs[n,1], self.sphereCoMs[n,2],asum, dR0, Rs, dx,dy,dz))    
    #make summary plots
    self.distortionPlot=plotWindow(self)
    self.distortionPlot.win.show()
    #plot distortion map
    self.distortionPlot.plotdx()


    
           
  def updateROICenters(self):
    '''updates ROIs positions from main window assuming user has center ROIs about phantom origin'''
    for n, roi in enumerate(self.currentROIs.ROIs): #start by assuming center of masses are at nominal positions
      self.sphereCenters[n,0]=roi.Xcenter
      self.sphereCenters[n,1]=roi.Ycenter
      self.sphereCenters[n,2]=roi.Zcenter
    self.croppedImage= self.extractFidSpheres(self.image,r=10)
    self.viewImage()
    self.plotROIs()

  def findSphereCenters(self,xi,S):
    params,paramlist= GaussianFit.initialize(xi,S)
    #params.update_constraints()
    result=GaussianFit.fit(params,xi,S)
    x0=result.params['x0'].value
    sigma=result.params['sigma'].value
    #print lmfit.fit_report(result)
    return x0, sigma

        
  def translateKImage(self,x0,y0,z0):
    n=int(self.ofovX/self.pixdx/2)  #subarray number of rows
    nref=int(self.refImage.shape[0]/2)
    tx=np.exp(1j*x0*self.kx)   #translations
    ty=np.exp(1j*y0*self.ky)
    tz=np.exp(1j*z0*self.kz)
    kimage=self.kImage*tx*ty*tz
    synImageComplex=np.flip(np.transpose(np.fft.fftshift(np.fft.ifft2(kimage))))    #2d fft then shift so zero frequency is in center
    im=np.absolute(synImageComplex[n-nref:n+nref,n-nref:n+nref])
    im=self.nImage(im)
    return im, kimage
    
  def transFTAlignment (self, stepsize=0.1):
      '''Translates over array with stepsize displacements and maximizes correlation with image,
      All ROIs are shifted to give maximum correlation, self.image is the new synthetic image'''
      sr=10
      image, kimage=self.translateKImage(0,0)
      self.difference= self.difIm(image, self.refImage)  
      for i in range(-sr,sr+1):
        for j in range (-sr,sr+1):
          x0=i*stepsize #xoffset
          z0=j*stepsize #zoffset
          im, kim=self.translateKImage(x0,z0)
          dif = self.difIm(im, self.refImage)
          if dif < self.difference:
            self.difference=dif
            self.xshift=x0
            self.zshift=z0
            image=im
            kimage=kim
      for roi in self.currentROIs.ROIs:
        roi.Xcenter+=self.xshift
        roi.Zcenter+=self.zshift
      self.kImage=kimage
      return image
    
  def makeRealImage(self, kimage, filter=1.0):
      n=int(self.ofovX/self.pixdx/2)  #subarray number of rows
      nref=int(self.refImage.shape[0]/2)
      synImageComplex=np.flip(np.transpose(np.fft.fftshift(np.fft.ifft2(kimage*filter))))    #2d fft then shift so zero frequency is in center
      im=np.absolute(synImageComplex[n-nref:n+nref,n-nref:n+nref])
      im=self.nImage(im)
      return im 
    
  def plotROIs(self):
    self.bPlotROIs=True
    if self.viewPlane=='Coronal':
      yi=self.vCenter(np.array([0,self.imv.currentIndex,0]))[1]
      for n, roi in enumerate(self.currentROIs.ROIs):
        d1=roi.d1  
        x=self.sphereCenters[n,0]
        y=self.sphereCenters[n,1]
        z=-self.sphereCenters[n,2]    #display has z inverted
        roipen= fn.mkPen(roi.color)
        roipen.setWidth(2)
        if np.absolute(y-yi)<20:
          r = fCircleROI((x-d1/2,z-d1/2),d1,roi.Label, labelcolor=pg.mkColor(roi.color),pen=roipen, movable=False)
          self.imv.getView().addItem(r)
          self.pgROIs.append(r)
          self.imv.getView().addItem(r.label)
    if self.viewPlane=='Axial':
      zi=self.vCenter(np.array([0,0,self.imv.currentIndex]))[2]
      for n, roi in enumerate(self.currentROIs.ROIs):  
        d1=roi.d1  
        x=self.sphereCenters[n,0]
        y=self.sphereCenters[n,1]
        z=self.sphereCenters[n,2]    #display has z inverted
        roipen= fn.mkPen(roi.color)
        roipen.setWidth(2)
        if np.absolute(z-zi)<20:
          r = fCircleROI((x-d1/2,y-d1/2),d1,roi.Label, labelcolor=pg.mkColor(roi.color),pen=roipen, movable=False)
          self.imv.getView().addItem(r)
          self.pgROIs.append(r)
          self.imv.getView().addItem(r.label)
    if self.viewPlane=='Sagittal':
      xi=self.vCenter(np.array([self.imv.currentIndex,0,0]))[0]
      for n, roi in enumerate(self.currentROIs.ROIs):  
        d1=roi.d1
        x=self.sphereCenters[n,0]
        y=self.sphereCenters[n,1]
        z=self.sphereCenters[n,2]    #display has z inverted
        roipen= fn.mkPen(roi.color)
        roipen.setWidth(2)
        if np.absolute(x-xi)<20:
          r = fCircleROI((z-d1/2,y-d1/2),d1,roi.Label, labelcolor=pg.mkColor(roi.color),pen=roipen, movable=False)
          self.imv.getView().addItem(r)
          self.pgROIs.append(r)
          self.imv.getView().addItem(r.label)

  def clearROIs(self): #remove ROIs
      self.bPlotROIs=False
      for roi in self.pgROIs:
          self.imv.getView().removeItem(roi)
          self.imv.getView().removeItem(roi.label)
      self.pgROIs=[]

  def setROId(self):
       d1=self.getDouble(text1='ROI diameter', text2='Input ROI diameter(mm)', default=10, min=0)
       if d1!= None:
        self.clearROIs()
        for n, roi in enumerate(self.currentROIs.ROIs):  
           roi.d1=d1
        self.plotROIs()
          
  def indexChanged(self):
    '''when stack index changes, replot ROIs'''
    self.clearROIs()
    self.plotROIs()

  def view3d(self, image3D):
      '''creates 3d rendering of current image stack, scale increase brightness'''
      self.Image3D=image3D
      if self.bViewSphereCenters:
        self.sphereCenterImage[0,0,0]=2*np.amax(image3D)
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
      data=self.scale3D*image3D.astype(float) /float(image3D.max())  #normalize data to 1
      if self.bViewSphereCenters:
        mask=2*self.scale3D*self.sphereCenterImage.astype(float) /float(self.sphereCenterImage.max())
      if self.view3Dinvert :
        data=1-data
      if self.bViewSphereCenters:
        d2 = np.empty(data.shape + (4,), dtype=np.ubyte)
        d2[..., 0] = data * self.view3DColor.red()
        d2[..., 1] = data * self.view3DColor.green()
        d2[..., 2] = mask * self.view3DColor.blue()
        d2[..., 3] = (data)**self.view3DTransparency * 255.   #sets transparency  
        d2[:, 0:3, 0:3] = [255,0,0,20]   #draw axes at corner of box 
        d2[0:3, :, 0:3] = [0,255,0,20]
        d2[0:3, 0:3, :] = [0,0,255,20] 
      else:
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
    
  def v3DColor(self):  
    self.view3DColor = QColorDialog.getColor()
    self.view3d(self.Image3D)
    
  def v3DTransparency(self):  
    t  =  self.getDouble(text1='3D transparency', text2='Enter value 0 (solid) to 10 transparent', default=2, min=0, max=10, decimals=2)
    if t!=None:
      self.view3DTransparency = t
      self.view3d(self.Image3D)

  def v3DIntensity(self):  
    t =  self.getDouble(text1='3D intensity', text2='Enter value 0 (dark) to 3 bright', default=2, min=0, max=10, decimals=2)
    if t!=None:
      self.scale3D = t
      self.view3d(self.Image3D)
                   
  def mouseMoved(self,evt): 
    '''mouse move event to move crosshairs and display location and values'''
    pos = evt[0]  ## using signal proxy turns original arguments into a tuple
    if self.imv.view.sceneBoundingRect().contains(pos):
      mousePoint = self.imv.view.vb.mapSceneToView(pos)
      if self.viewPlane=='Coronal':
        xp=mousePoint.x()
        zp=-mousePoint.y()
        yp=self.vCenter(np.array([0,self.imv.currentIndex,0]))[1]   #convert index to position
      if self.viewPlane=='Axial':
        xp=mousePoint.x()
        yp=-mousePoint.y()
        zp=self.vCenter(np.array([0,0,self.imv.currentIndex]))[2]   #convert index to position
      if self.viewPlane=='Sagittal':
        zp=mousePoint.x()
        yp=-mousePoint.y()
        xp=self.vCenter(np.array([self.imv.currentIndex,0,0,]))[0]   #convert index to position
      ijk=self.rtoIndex([xp,yp,zp])
      try:
        if abs(mousePoint.x()) < self.fovh/2 and abs(mousePoint.y()) < self.fovh/2:
          value=self.currentArray[int(ijk[0]),int(ijk[1]),int(ijk[2])] 
          #self.label.setText("x={:.2f}, y={:.2f}, value={:.2f}".format(mousePoint.x(),mousePoint.y(), value ))
          self.win.setWindowTitle("i,j,k={:.2f},{:.2f},{:.2f}; x={:.2f}, y={:.2f}, z={:.2f}, value={:.2f}".format(ijk[0],ijk[1],ijk[2],xp,yp,zp, value ))
        else:
          self.win.setWindowTitle(self.windowTitle)
      except:
        pass
      self.vLine.setPos(mousePoint.x())
      self.hLine.setPos(mousePoint.y())
    
      

      
class fCircleROI(pg.EllipseROI):
    """Defines a circular ROI using pyqtgraph's EllipseROI"""
    def __init__(self, pos, size,label,labelcolor=(0,0,0), **args):   #passing in calling form, could be a problem
        pg.ROI.__init__(self, pos, size, **args)
        lblColor=(255,0,255)
        self.aspectLocked = True
        self.label = pg.TextItem(label, labelcolor, anchor = (0,0))
        self.label.setPos(pos[0],pos[1])
        self.path=None    #added when updating to pyqtgraph 0.11 Do not know why it is not there
        
class messageWindow():
  def __init__(self, image=None, parent = None):
    '''Defines message window, contains output data'''    
    self.win = QMainWindow()
    #self.win.setWindowFlags(Qt.WindowStaysOnTopHint)        #having the window always on top can be a problem
    self.win.resize(1500,1000)
    self.text=QTextEdit()
    self.win.setCentralWidget(self.text)
    self.win.setWindowTitle('Fiducial analysis')
    self.menu = self.win.menuBar()
    
    self.fileMenu = self.menu.addMenu('&file')
    self.actionSaveFile = QAction('Save to file', self.win)
    self.actionSaveFile.triggered.connect(self.saveFile)
    self.fileMenu.addAction(self.actionSaveFile)
       
  def saveFile(self):
    f = QFileDialog.getSaveFileName(None, "Report File Name", '',".txt")
    if not f:  #if cancel is pressed return
      return None     
    if type(f)==tuple:    #passes  string with PyQt4 and a tuple with PyQt5
      fileName=f[0]
    else:
      fileName=f
    file= open(fileName, 'w')
    file.write(self.text.toPlainText())
    file.close()
      
  
  def message(self, s):
    self.text.insertPlainText(s)     
 
class plotWindow(QMainWindow):
      def __init__(self, pw, image=None, parent=None):
            '''Defines window for plotting data, pw is the parent window'''    
            super(plotWindow,self).__init__()
            self.dplot = pg.PlotWidget()
            self.win=self
            self.win.setCentralWidget(self.dplot)
            self.win.resize(1000,500)
            self.win.setWindowTitle('Fiducial Data')
            self.pw=pw
            self.menu = self.menuBar()
            self.BG='black'  #flag to set background color
            self.bgLabelColor='#FFF'
            self.tfont=QFont('Times', 20) #Times, Helvetica
            self.ticktextSize=20    #tick text size
            self.dplot.getAxis("bottom").tickFont = self.tfont
            self.dplot.getAxis("bottom").setStyle(tickTextOffset=15, tickTextHeight=100, tickTextWidth=100)
            self.dplot.getAxis("left").setStyle(tickTextOffset=5)
            self.dplot.getAxis("left").setWidth(90)
            self.dplot.getAxis("left").tickFont = self.tfont
            self.dplot.getAxis("right").tickFont = self.tfont
            self.yLabelSize='24pt'
            self.xLabelSize='24pt'
            self.titleSize='24pt'
            self.dplot.showGrid(x=True, y=True)
            self.dplot.setBackground(background='k')
            self.symbolSize=16
            self.symb=['o', 't', 's', 'd',  '+']  #symbol list for multiple plots circle, triangle, start, diamond
            self.symbolColor=['k', 'b', 'r','g', 'c', 'm', 'k']
            self.nplot=0
            self.penm = fn.mkPen(255, 0, 255)
            cc='rgb'+str(SystemPhantom.fidColor['center'])
            cl='rgb'+str(SystemPhantom.fidColor['left'])
            cr='rgb'+str(SystemPhantom.fidColor['right'])
            cp='rgb'+str(SystemPhantom.fidColor['posterior'])
            ca='rgb'+str(SystemPhantom.fidColor['anterior'])
            cs='rgb'+str(SystemPhantom.fidColor['superior'])
            ci='rgb'+str(SystemPhantom.fidColor['inferior'])
            self.dpTitle='<span style="color:'+cc+';">Center, <\span><span style="color:'+cs+';">Superior,<\span> <span style="color:'+ci+';">Inferior,\
              <\span><span style="color:'+cr+';"> Right,<\span><\span><span style="color:'+cl+';"> Left,<\span>\
              <\span><span style="color:'+ca+';"> Anterior,<\span><\span><span style="color:'+cp+';"> Posterior<\span>'
            self.imageMenu = self.menu.addMenu('&File')    
            self.actionSaveData = QAction('Save data', self.win)
            self.imageMenu.addAction(self.actionSaveData)
            self.actionSaveData.triggered.connect(self.saveData)
            self.imageMenu = self.menu.addMenu('&Data')    
            self.actionPlotdx = QAction('Plot dx', self.win)
            self.imageMenu.addAction(self.actionPlotdx)
            self.actionPlotdx.triggered.connect(self.plotdx)
            self.actionPlotdy = QAction('Plot dy', self.win)
            self.imageMenu.addAction(self.actionPlotdy)
            self.actionPlotdy.triggered.connect(self.plotdy)
            self.actionPlotdz = QAction('Plot dz', self.win)
            self.imageMenu.addAction(self.actionPlotdz)
            self.actionPlotdz.triggered.connect(self.plotdz)
            self.actionPlotIntensity = QAction('Plot sphere intensity', self.win)
            self.imageMenu.addAction(self.actionPlotIntensity)
            self.actionPlotIntensity.triggered.connect(self.plotIntensity)
            self.actionPlotVolume = QAction('Plot volume disortion', self.win)
            self.imageMenu.addAction(self.actionPlotVolume)
            self.actionPlotVolume.triggered.connect(self.plotVolume)
            self.actionPlotCoMminusConv = QAction('Plot CoM-convolution centers', self.win)
            self.imageMenu.addAction(self.actionPlotCoMminusConv)
            self.actionPlotCoMminusConv.triggered.connect(self.plotCoMminusConv)
            self.actionShowSpreadsheet = QAction('Show spreadsheet', self.win)
            self.imageMenu.addAction(self.actionShowSpreadsheet)
            self.actionShowSpreadsheet.triggered.connect(self.showSpreadsheet)
            
            self.imageMenu = self.menu.addMenu('&Plot')    
            self.actionToggleBG = QAction('Toggle background', self.win)
            self.imageMenu.addAction(self.actionToggleBG)
            self.actionToggleBG.triggered.connect(self.toggleBackground)
            self.actionSymbolSize = QAction('Symbol Size', self.win)
            self.imageMenu.addAction(self.actionSymbolSize)
            self.actionSymbolSize.triggered.connect(self.setSymbolSize)
            self.actionLabelSize = QAction('Label Size', self.win)
            self.imageMenu.addAction(self.actionLabelSize)
            self.actionLabelSize.triggered.connect(self.setLabelSize)
            
      def plotdx(self):
        xdata=self.pw.sphereRCenters
        ydata=self.pw.geoDistortion[:,0]
        ylabel="L/R distortion dx (mm)"
        sd=np.std(ydata)
        wintitle='Fiducial Geometric Disortion:X-direction, standard deviation(mm)= ' + '{:.4f}, '.format(sd)
        self.plotData(xdata, ydata, ylabel=ylabel,wintitle=wintitle,clearPlot=True,symbolBrush=self.pw.fiducialColors)
      def plotdy(self):
        xdata=self.pw.sphereRCenters
        ydata=self.pw.geoDistortion[:,1]
        ylabel="A/P distortion dy (mm)"
        sd=np.std(ydata)
        wintitle='Fiducial Geometric Disortion:Y-direction, standard deviation(mm)= ' + '{:.3f}, '.format(sd)
        self.plotData(xdata, ydata, ylabel=ylabel,wintitle=wintitle, clearPlot=True,symbolBrush=self.pw.fiducialColors)  
      def plotdz(self):
        xdata=self.pw.sphereRCenters
        ydata=self.pw.geoDistortion[:,2]
        ylabel="S/I distortion dz (mm)"
        sd=np.std(ydata)
        wintitle='Fiducial Geometric Disortion:Z-direction, standard deviation= ' + '{:.3f}'.format(sd)
        self.plotData(xdata, ydata, ylabel=ylabel,wintitle=wintitle, clearPlot=True,symbolBrush=self.pw.fiducialColors)      
      def plotIntensity(self):
        xdata=self.pw.sphereRCenters
        ydata=self.pw.sphereIntensity/np.amax(self.pw.sphereIntensity)
 
        ylabel="Integrated Sphere Intensity"
        m=np.mean(ydata)
        std=np.std(ydata)
        wintitle='Mean={:.3f} ,std={:.3f}:'.format(m,std) 
        self.plotData(xdata, ydata, ylabel=ylabel, wintitle=wintitle, clearPlot=True, symbolBrush=self.pw.fiducialColors)
      def plotVolume(self):
        xdata=self.pw.sphereRCenters
        ydata=self.pw.sphereVolumes
        ylabel="Relative volume"
        wintitle='Fiducial relative volume:'
        self.plotData(xdata, ydata, ylabel=ylabel, wintitle=wintitle, clearPlot=True, symbolBrush=self.pw.fiducialColors)
      def plotCoMminusConv(self):
        xdata=self.pw.sphereRCenters
        ylabel="CoM-Convolution center (&mu;m)"
        wintitle='Fiducial relative volume:'
        self.plotData(xdata, self.pw.CoMminusCovCenter[:,0]*1000, ylabel=ylabel, wintitle=wintitle, clearPlot=True, symbolBrush=self.pw.fiducialColors)
        self.plotData(xdata, self.pw.CoMminusCovCenter[:,1]*1000, ylabel=ylabel, wintitle=wintitle, clearPlot=False, symbolBrush=self.pw.fiducialColors)
        self.plotData(xdata, self.pw.CoMminusCovCenter[:,2]*1000, ylabel=ylabel, wintitle=wintitle, clearPlot=False, symbolBrush=self.pw.fiducialColors)
        
      def showSpreadsheet(self):
        self.pw.viewMessages()

      def plotData(self,x,y,ylabel='ylabel',wintitle='', clearPlot=False,symbolBrush=None,xlabelsize='18pt'):
        if clearPlot:
          self.dplot.clear()
        self.dplot.plot(x, y, pen=None, symbolBrush=symbolBrush, symbol=self.pw.fiducialSymbols,symbolSize=self.symbolSize)
        self.setWindowTitle(wintitle +ylabel +', ' + self.pw.filename)
        self.dplot.setTitle(self.dpTitle,size=self.titleSize)
        self.tfont.setPixelSize(self.ticktextSize)
        self.dplot.getAxis("bottom").setTickFont(self.tfont)
        self.dplot.getAxis("bottom").setStyle(tickTextOffset = int(self.ticktextSize/2))
        self.dplot.getAxis("left").setTickFont(self.tfont)
        self.dplot.getAxis("left").setStyle(tickTextOffset = int(self.ticktextSize/2))
        self.dplot.setLabel('bottom',"|R-R<sub>0</sub>|(mm)",color=self.bgLabelColor, font=self.xLabelSize)
        self.dplot.setLabel('left',ylabel,color=self.bgLabelColor, font=self.yLabelSize )
        self.dplot.getAxis('right').setStyle(showValues=False)
        self.dplot.showAxis('right')
        self.dplot.getAxis('top').setStyle(showValues=False)
        self.dplot.showAxis('top')

          
      def saveData(self):        
           pass
         
      def toggleBackground(self):
        if self.BG=='black':  #switch to white background
          self.dplot.setBackground('w')
          self.dplot.getAxis("bottom").setPen('k')
          self.dplot.getAxis("bottom").setTextPen('k')
          self.dplot.getAxis("left").setPen('k')
          self.dplot.getAxis("left").setTextPen('k')
          self.dplot.getAxis("right").setPen('k')
          self.BG='white'
          self.bgLabelColor='k'
          self.textPen=self.penm
        else: #switch to black background
          self.dplot.setBackground('k')
          self.dplot.getAxis("bottom").setPen('w')
          self.dplot.getAxis("bottom").setTextPen('w')
          self.dplot.getAxis("left").setPen('w')
          self.dplot.getAxis("left").setTextPen('w')
          self.dplot.getAxis("right").setPen('w')
          self.bgLabelColor='w'
          self.BG='black'
        
      def setSymbolSize(self):
        fsize,ok = QInputDialog.getInt(self, 'Input symbol', 'size', value=self.symbolSize, min=0, max=100, step=1)
        if ok:
          self.symbolSize=fsize
        
      def setLabelSize(self):
        lsize,ok = QInputDialog.getInt(self, 'Input X-label', 'size', value=18, min=0, max=100, step=1)
        if ok:
          self.xLabelSize=str(lsize)+'pt'
        lsize,ok = QInputDialog.getInt(self, 'Input Y-label', 'size', value=18, min=0, max=100, step=1)
        if ok:
          self.yLabelSize=str(lsize)+'pt'
                             
if __name__ == '__main__':
    app = QApplication(sys.argv)
    main = fiducialWindow()
    main.win.show() 
    sys.exit(app.exec_())    


