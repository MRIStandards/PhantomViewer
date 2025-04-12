'''
Created on April 29, 2017
Classes and routines to analyze ACR-like resolution inset and compute resolution
   
@author: stephen russek
'''

# -*- coding: utf-8 -*-
import sys
import numpy as np
import copy
from PyQt5 import QtGui, QtCore
import pydicom    #import pydicom to read DICOM data sets and DICOMDIR
import pyqtgraph as pg    #uses pyqtgraph PlotWidget for graphs and ImageView windows for images, http://www.pyqtgraph.org/documentation/introduction.html
import pyqtgraph.functions as fn
import ImageList
import scipy.ndimage
import scipy.optimize
import skimage.measure
import matplotlib.pyplot as plt
import scipy.ndimage
import scipy.optimize
import skimage.measure
import matplotlib.pyplot as plt

#parameters defining ACR-type resolution inset
RI1UR = np.array([[0.0,0.0],[10.58,.0],[19.76,0],[27.56,0], [33.96,0]] )    # upper right coordinate in mm
RI1d=np.array([0.8,0.7,0.6,0.5,0.4])  #hole diameters in mm
RI1angle=10.0   #angle of angled array in degrees
RInholes=155    #number of holes in ACR resolution inset
RIfrw=45.0      #resolution inset frame width
RIfrh=14        #resolution inset frame height
RIframeLL=np.array([-RIfrw/2,-RIfrh/2])
RIsize=np.array([RIfrw,RIfrh])
RIxoffset=16.5


def createResArray():
  holearray=[]  #array of holes (x,y,diameter) all in mm
  for n,d in enumerate(RI1d):
    for i, j in np.ndindex((4,4)):
      holearray.append(np.array([-i*2*d+RI1UR[n,0]-RIxoffset,j*2*d+RI1UR[n,1],d]))
    for i, j in np.ndindex((4,4)):
      x,y=rotateArray(i*2*d,j*2*d,RI1angle+90.0)
      if not(j==0 and i==0):
        holearray.append(np.array([x+RI1UR[n,0]-RIxoffset,y+RI1UR[n,1],d]))
  return holearray

def rotateArray(x,y,angle): 
  theta=angle*np.pi/180.0
  xr=np.cos(theta)*x+np.sin(theta)*y
  yr=-np.sin(theta)*x+np.cos(theta)*y
  return xr,yr


                                          
class fCircleROI(pg.EllipseROI):
    """Defines a circular ROI using pyqtgraph's EllipseROI"""
    def __init__(self, pw, pos, size, **args):   #passing in calling form, could be a problem
        pg.ROI.__init__(self, pos, size, **args)
        self.aspectLocked = True
        self.pw=pw    #parent window

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
            if ev.button() == QtCore.Qt.LeftButton:
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
        if self.translatable and self.isMoving and ev.buttons() == QtCore.Qt.LeftButton:
            snap = True if (ev.modifiers() & QtCore.Qt.ControlModifier) else None
            newPos = self.mapToParent(ev.pos()) + self.cursorOffset
            self.translate(newPos - self.pos(), snap=snap, finish=False)
            #self.ROIViewInstance.translateROIs(newPos - self.pos(),snap, self.Index)

    def setMouseHover(self, hover):
        '''Inform the ROI that the mouse is(not) hovering over it'''
        if self.mouseHovering == hover:
            return
        self.mouseHovering = hover
        if hover:
            self.currentPen = fn.mkPen(255, 255, 0)

        else:
            self.currentPen = fn.mkPen(255, 0, 0)
        self.update() 
        
class fRectROI(pg.RectROI):
    """Defines a rectangular ROI using pyqtgraph's RectROI"""
    def __init__(self,pw, pos, size,   **args):   
        pg.ROI.__init__(self, pos, size, **args)
        self.aspectLocked = False
        self.pw=pw    #parent window

    def mouseDragEvent(self, ev):       #Dragging ROI event: translates ROIs
        if ev.isStart():  
            if ev.button() == QtCore.Qt.LeftButton:
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
        if self.translatable and self.isMoving and ev.buttons() == QtCore.Qt.LeftButton:
            snap = True if (ev.modifiers() & QtCore.Qt.ControlModifier) else None
            newPos = self.mapToParent(ev.pos()) + self.cursorOffset
            #self.translate(newPos - self.pos(), snap=snap, finish=False)
            self.pw.translateROIs(newPos - self.pos())


    def setMouseHover(self, hover):
        '''Inform the ROI that the mouse is(not) hovering over it'''
        if self.mouseHovering == hover:
            return
        self.mouseHovering = hover
        if hover:
            self.currentPen = fn.mkPen(255, 255, 0)

        else:
            self.currentPen = fn.mkPen(255, 0, 0)
            #self.ROIViewInstance.currentROI=-1
        self.update() 

class messageWindow():
  def __init__(self, image=None, parent = None):
    '''Define message window'''    
    self.win = QtGui.QMainWindow()
    #self.win.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)        #having the window always on top can be a problem
    self.win.resize(500,200)
    self.text=QtGui.QTextEdit()
    self.win.setCentralWidget(self.text)
    self.win.setWindowTitle('Resolution analysis') 
    
  def message(self, s):
    self.text.insertPlainText(s)     
 
class resolutionWindow():
  def __init__(self, image=None, parent = None):
    '''Define  window to view and compare resolution images: original, reference, synthetic, and difference images'''    
    self.win = QtGui.QMainWindow()
    #self.win.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)
    self.win.resize(800,600)
    self.imv = pg.ImageView( view = pg.PlotItem())
    #self.imv.getView().invertY(False)
    #self.imv.ui.histogram.plot.setLogMode(None,True)    #set the histogram y axis to a log scale    
    self.imv.ui.roiBtn.setText("Line scan/ROI")
    self.win.setCentralWidget(self.imv)
    self.win.setWindowTitle('Resolution Inset')
    #point = self.rect().topRight()
    #self.win.move(point + QtCore.QPoint(self.width()/2, 0)) 
    self.menu = self.win.menuBar()
    self.penr = fn.mkPen(255, 0, 0)
    self.peng = fn.mkPen(0, 255, 0)
    self.peny = fn.mkPen(255, 255, 0)
    self.imv.getView().setLabel('bottom',"H","mm")
    self.imv.getView().setLabel('left',"V","mm")
    self.messages=messageWindow()
    #self.messages.win.show()

    self.menu = self.win.menuBar()
    
    self.imageMenu = self.menu.addMenu('&Images')
    self.actionSelectImages = QtGui.QAction('Select/Add Images', self.win)
    self.actionSelectImages.setShortcut('Ctrl+S')
    self.actionSelectImages.setStatusTip('Select/Add Images')
    self.actionSelectImages.triggered.connect(self.openFile)
    self.imageMenu.addAction(self.actionSelectImages) 
    
    self.actionClear_All_Images = QtGui.QAction('Clear All Images', self.win)
    self.actionClear_All_Images.setShortcut('Ctrl+C')
    self.actionClear_All_Images.setStatusTip('Clear All Images')
    #self.actionClear_All_Images.triggered.connect(clearImages)
    self.imageMenu.addAction(self.actionClear_All_Images)
    
    self.actionViewImage = QtGui.QAction('View Original Image', self.win)
    self.actionViewImage.triggered.connect(self.viewImage)
    self.imageMenu.addAction(self.actionViewImage)
    
    self.actionViewSynImage = QtGui.QAction('View Synthetic Image', self.win)
    self.actionViewSynImage.triggered.connect(self.viewSynImage)
    self.imageMenu.addAction(self.actionViewSynImage)

    self.actionViewRefImage = QtGui.QAction('View Reference Image', self.win)
    self.actionViewRefImage.triggered.connect(self.viewRefImage)
    self.imageMenu.addAction(self.actionViewRefImage)
    
    self.actionViewDifImage = QtGui.QAction('View Difference Image', self.win)
    self.actionViewDifImage.triggered.connect(self.viewDifImage)
    self.imageMenu.addAction(self.actionViewDifImage)
    
    self.actionViewkImage = QtGui.QAction('View k-space Image', self.win)
    self.actionViewkImage.triggered.connect(self.viewkImage)
    self.imageMenu.addAction(self.actionViewkImage)
        
    self.actionSaveRIFile = QtGui.QAction('Save RI file', self.win)
    self.actionSaveRIFile.triggered.connect(self.saveRIFile)
    self.imageMenu.addAction(self.actionSaveRIFile)

    self.imageMenu = self.menu.addMenu('&Analysis')
    self.actionOptimizeSynImage = QtGui.QAction('Align Synthetic Image', self.win)
    self.actionOptimizeSynImage.triggered.connect(self.alignSyntheticImage)
    self.imageMenu.addAction(self.actionOptimizeSynImage)
             
    self.actioncreatePCImage = QtGui.QAction('Create Point Cloud Synthetic Image', self.win)
    self.actioncreatePCImage.triggered.connect(self.createPCImage)
    self.imageMenu.addAction(self.actioncreatePCImage)
    
    self.actionCreateKImage = QtGui.QAction('Create FT(circle) Synthetic Image', self.win)
    self.actionCreateKImage.triggered.connect(self.createKspaceImage)
    self.imageMenu.addAction(self.actionCreateKImage)
        
    self.actionFindRes = QtGui.QAction('Find Resolution', self.win)
    self.actionFindRes.triggered.connect(self.blurFTImages)
    self.imageMenu.addAction(self.actionFindRes)
    
    self.imageMenu = self.menu.addMenu('&ROIs')
    self.actionPlotROIs = QtGui.QAction('Plot/Hide ROIs', self.win)
    self.actionPlotROIs.triggered.connect(self.plotROIs)
    self.imageMenu.addAction(self.actionPlotROIs)
            
    self.pgROIs=[]    #list of ROIs
    self.riview=46.0   #dimension in mm of reference and synthetic RI images
    self.xshift=0 #measures shift in ROIs required to match synthetic and reference images
    self.zshift=0
    self.resolution=0   #calculated resolution in mm
     
    if image ==None:
      self.image= np.zeros((100,100))
    else:
      self.image=image
    self.imv.setImage(self.image,pos = (-1,-1), scale = (2,2),)
    
  def getDouble(self, text1='QInputDialog', text2='input floating point value', default=0.0, min=0, max=100, decimals=2):
      double,ok = QtGui.QInputDialog.getDouble(self.win, text1,text2, default, min, max, decimals)
      if ok:
        return double
      else:
        return None
      
  def plotOverlay(self,holearray):
    self.holearray=holearray
    for hole in holearray:
      r = MyCircleOverlay(pw=self,pos=(hole[0],hole[1]), size=hole[2], pen=self.peny, movable=True)
      self.imv.getView().addItem(r)
    r = pg.ROI(pos = (RIframeLL[0], RIframeLL[1]), size=(RIsize[0], RIsize[1]), pen=self.penr, movable=True)
    self.imv.getView().addItem(r)

  def clearROIs(self):
      for roi in self.pgROIs:
        self.imv.getView().removeItem(roi)
      self.pgROIs=[]
               
  def plotROIs(self):
    if len(self.pgROIs)==0: #if no ROIs plot them
      pen=self.penr
      for roi in self.currentROIs.ROIs:
        if roi.Type=='Sphere':
          if roi.Index==75:
            pen=self.peny
          else:
            pen=self.penr
          r = fCircleROI(pw=self, pos=(roi.Xcenter-roi.d1/2,roi.Zcenter-roi.d1/2), size=roi.d1, pen=pen, movable=True)
          self.imv.getView().addItem(r)
          self.pgROIs.append(r)
        if roi.Type=='Rectangle':
          r = fRectROI(self,[roi.Xcenter-roi.dx/2,roi.Zcenter-roi.dy/2], [roi.dx,roi.dy],angle=roi.theta, pen=self.peny, movable=True)
          self.imv.getView().addItem(r)
          self.pgROIs.append(r)
    else: #remove ROIs
      for roi in self.pgROIs:
          self.imv.getView().removeItem(roi)
      self.pgROIs=[]

  def translateROIs(self, pos):
    for roi in self.currentROIs.ROIs:
       roi.Xcenter+=pos[0]
       roi.Zcenter+=pos[1]
    for roi in self.pgROIs:
      roi.translate(pos)
  
                
  def viewImage(self): 
    self.win.setWindowTitle('Resolution Inset: Cropped Image')
    self.imv.setImage(self.image,pos = (-self.fovX/2,-self.fovY/2), scale = (self.fovX/self.image.shape[0],self.fovY/self.image.shape[1]),)
    self.imv.getView().setLabel('bottom','X',"mm")
    self.imv.getView().setLabel('left','Y',"mm")
    
  def viewSynImage(self, res=0): 
    dif=self.difIm(self.synImage, self.refImage)
    title='Resolution Inset: Synthetic Image, x,z shift(mm)=' 
    title+="{:.2f}".format(self.xshift) +  ',' + "{:.2f}".format(self.zshift) +', |Dif|=' + "{:.2f}".format(dif) 
    title+=', resolution='+"{:.2f}".format(self.resolution)
    self.messages.message(title+'\n')
    self.win.setWindowTitle(title)
    self.imv.setImage(self.synImage,pos = (-self.riview/2,-self.riview/2), scale = (self.riview/self.synImage.shape[0],self.riview/self.synImage.shape[1]),)
    self.imv.getView().setLabel('bottom','X',"mm")
    self.imv.getView().setLabel('left','Y',"mm")
        
  def viewRefImage(self): 
    self.win.setWindowTitle('Resolution Inset: Reference Image:')
    self.imv.setImage(self.refImage,pos = (-self.riview/2,-self.riview/2), scale = (self.riview/self.synImage.shape[0],self.riview/self.synImage.shape[1]),)
    self.imv.getView().setLabel('bottom','X',"mm")
    self.imv.getView().setLabel('left','Y',"mm")

  def viewDifImage(self): 
    dif=self.difIm(self.synImage, self.refImage)
    self.win.setWindowTitle('Resolution Inset: Difference Image, |Diff|=' + "{:.2f}".format(dif))
    self.imv.setImage(self.difImage,pos = (-self.riview/2,-self.riview/2), scale = (self.riview/self.synImage.shape[0],self.riview/self.synImage.shape[1]),)
    self.imv.getView().setLabel('bottom','X',"mm")
    self.imv.getView().setLabel('left','Y',"mm")

  def viewkImage(self): 
    title='Resolution Inset: k-Space Image, x,z shift(mm)=' 
    self.messages.message(title+'\n')
    self.win.setWindowTitle(title)
    self.imv.setImage(np.absolute(self.kImage))
    self.imv.getView().setLabel('bottom','kx',"1/mm")
    self.imv.getView().setLabel('left','ky',"1/mm")
        
  def addImage(self,image,cROIs,fovx=1.0,fovy=1.0, ofovx=200, ofovy=200):
    '''passes in a cropped image of the resolution inset, with its field of view and the original image field of view''' 
    self.image=self.nImage(image)   #normalize input image
    self.currentROIs=copy.deepcopy(cROIs)
    for roi in self.currentROIs.ROIs:   #flipping Zcenter to have upper left at -,- with rotated holes on top
      roi.Zcenter=-roi.Zcenter
    self.fovX=fovx      #field of view of he cropped image
    self.fovY=fovy 
    self.ofovX=ofovx    #field of view of original image
    self.ofovY=ofovy
    self.pixdx=self.fovX/self.image.shape[0]
    self.pixdy=self.fovX/self.image.shape[1]       
    self.win.setWindowTitle('Resolution Inset: Calculating reference and synthetic images ...')
    self.refImage=self.makeRefImage()   #make a reference image by cropping input image to self.riview
    self.createKspaceImage()
    #self.calculatePointClouds(sigma=0) #initial point clouds used to calculate synthetic image
    #self.synImage=self.makeSynImage(0,0)    #calculates synthetic image using 0,0 offset of the resolution inset
    self.difImage=self.synImage-self.refImage
    #display original image   
    self.imv.setImage(self.image,pos = (-self.fovX/2,-self.fovY/2), scale = (self.fovX/self.image.shape[0],self.fovY/self.image.shape[1])) #pos = (-self.fovX/2,self.fovY/2), scale = (self.fovX/self.image.shape[0],self.fovY/self.image.shape[1]),
    self.imv.getView().setLabel('bottom','X',"mm")
    self.imv.getView().setLabel('left','Y',"mm")
    self.win.setWindowTitle('Resolution Inset: Cropped RI Image')

  def nImage(self,image):
    '''normalizes an image to give L1norm=1000'''
    pav=np.sum(np.absolute(image))/1000
    if pav==0.0:
      print ("warning image=0")
      return image
    else:
      return image/pav

  def alignSyntheticImage(self):
      ss= self.getDouble('Enter translation step','translation step(mm)', 0.1, 0, 1, 2)
      #self.calculatePointClouds()
      self.win.setWindowTitle('Aligning Synthetic Image ......')
      self.messages.message('Aligning Synthetic Image, step size =' + str(ss))
      #im=self.makeSynImage(0,0)
      im=self.synImage
      self.initdifference = self.difIm(im, self.refImage) #measures difference using L1 norm
      self.difference=self.initdifference 
      self.messages.message( ', initial difference='+"{:.2f}".format(self.difference) +'\n')    
      self.synImage=self.transFTAlignment(stepsize=ss)
      self.difference = self.difIm(self.synImage, self.refImage)
      self.viewSynImage()
      self.clearROIs()
      self.plotROIs()
      return
     
  def transAlignment (self, stepsize=0.1):
      '''Translates over array with stepsize displacements and maximizes correlation with image,
      All ROIs are shifted to give maximum correlation, self.image is the new synthetic image'''
      sr=6  
      for i in range(-sr,sr+1):
        for j in range (-sr,sr+1):
          x0=i*stepsize #xoffset
          z0=j*stepsize #zoffset
          im=self.makeSynImage(x0,z0)
          dif = self.difIm(im, self.refImage)
          if dif < self.difference:
            self.difference=dif
            self.xshift=x0
            self.zshift=z0
            image=im
      for roi in self.currentROIs.ROIs:
        roi.Xcenter-=self.xshift
        roi.Zcenter-=self.zshift
      return image
      
  def createPCImage (self):
      '''Creates synthetic image from pointclouds, Blurs image with gaussian point spread function'''  
      sigma=self.getDouble('Enter blur', 'mm', default=0.1, min=0, max=10, decimals=2)
      self.calculatePointClouds(sigma=sigma)
      self.synImage=self.makePCSynImage(0,0)
      dif = self.difIm(self.synImage, self.refImage)
      self.messages.message( 'New synthetic image: blur=' + str(sigma) + ', dif=' + str(dif) + '\n')
      self.viewSynImage(res=self.resolution)
  
  def blurImages (self):
      '''Blurs image with gaussian point spread function'''  
      nps=40
      self.resArray=np.zeros((nps,2))
      blur=0.025  #blur in increments of 0.1 mm
      m='Calculating image resolution ......'
      self.win.setWindowTitle(m)
      self.messages.message(m)
      for i in range(0,nps):
        self.calculatePointClouds(sigma=i*blur)
        im=self.makeSynImage(0,0)
        dif = self.difIm(im, self.refImage)
        if i==0:
          dif0=dif
        if dif<dif0:
          self.synImage=im
          self.resolution=i*blur
          dif0=dif
        self.resArray[i,0]=i*blur
        self.resArray[i,1]=dif
      self.calculatePointClouds(sigma=self.resolution)
      self.resPlot=pg.plot(title='Resolution')
      self.resPlot.plot(self.resArray[:,0],self.resArray[:,1])        
      self.viewSynImage(res=self.resolution)
          
  def makeSynImage(self,x0, z0,imType='ft'):
    if imType=='ft':
      self.makeFTcircleArrayImage(x0=x0, y0=z0)
    if imType=='pc':
      self.makePCSynImage(x0, z0)
  
  def makePCSynImage(self,x0, z0):
      '''makes a synthetic image (2d numpy array) from ROIs with a shift of x0,z0
      Creates synthetic image from array of point clouds by histogramming points with given binsize and image range'''
      nbins=self.refImage.shape[0]
      a=self.riview/2
      imageRange=np.array([[-a,a],[-a,a]])
      im=np.zeros((self.refImage.shape[0], self.refImage.shape[1])) #start with a zero array 
      for roi in self.currentROIs.ROIs:
        if roi.Type=="Sphere":
          pc=self.selectpc(roi.d1)  #select pointcloud based on hole diameter
          roiCenter=np.array([roi.Xcenter-x0, roi.Zcenter-z0])
          tpc = pc + roiCenter    #shifting point cloud
          arr=np.histogram2d(tpc[:,0],tpc[:,1], bins=nbins, range=imageRange)[0]
          im += arr
      im=self.nImage(im)    #normalize image
      return im

  def makeRefImage(self):
    '''make sub array with resolution inset, with all non-close voxels set to zero'''
    xr=3  #include xr points on each side of the ROI
    yr=3
    dx=self.fovX/self.image.shape[0]  #pixel spacing in mm in input image
    dy=self.fovY/self.image.shape[1]
    xp=int(self.riview/dx/2)  #subarray number of rows
    yp=int(self.riview/dy/2)
    self.riview=xp*dx*2
    nr=int(self.image.shape[0]/2)
    nc=int(self.image.shape[1]/2)
    arr=np.zeros((2*xp,2*yp))
    for roi in self.currentROIs.ROIs:
      x=int(roi.Xcenter/dx)
      y=int(roi.Zcenter/dy)
      for i in range(-xr,xr+1):
        for j in range (-yr,yr+1):
          try:
            arr[xp+x+i,yp+y+j]=self.image[nc+x+i,nr+y+j]
          except:
            pass
    return self.nImage(arr)
              
  def calculatePointClouds(self, sigma=0):
      self.pcp8= self.pointCloud(d=0.8, sigma=sigma)
      self.pcp7= self.pointCloud(d=0.7, sigma=sigma)
      self.pcp6= self.pointCloud(d=0.6, sigma=sigma)
      self.pcp5= self.pointCloud(d=0.5, sigma=sigma)
      self.pcp4= self.pointCloud(d=0.4, sigma=sigma)

  def selectpc(self, d):
      '''selects a point cloud'''
      if d ==0.8:
        pc=self.pcp8
      if d ==0.7:
        pc=self.pcp7
      if d ==0.6:
        pc=self.pcp6
      if d ==0.5:
        pc=self.pcp5
      if d ==0.4:
        pc=self.pcp4
      return pc
          
  def correlateIm(self,im1,im2):
    return np.sum(im1*im2)
  
  def difIm(self,im1,im2):
    try:
      return (np.sum(np.square(im1-im2)))**0.5
    except:
      return 0 
    #return np.sum(np.absolute(im1-im2))    
            
  def pointCloud(self,pdensity=20000.0, d=1.0, dx=1.0, dy=1.0, sigma=0):
    '''generates a point cloud, pc, with n circular elements of diameter d in a space +-dx by +-dy,
     with gaussian noise sigma, all dimensions assumed to be in mm'''
    n=int(pdensity*np.pi*d**2/4)
    npc=0
    pc=np.zeros((n,2))
    if sigma==0:    #random.normal will bomb if sigma=0
      sigma=0.0001
    while npc<n:
      psf=np.random.normal(scale=sigma,size=2)
      p=2*dx*(np.random.rand(2)-0.5)
      r=p-psf
      if p[0]**2+p[1]**2<d**2/4:
        pc[npc]=r
        npc +=1
    return pc
  
  def createKspaceImage(self):
    '''Creates a kspace image by FT an array of holes'''
    self.makeFTcircleArrayImage()  #creates kspace and synthetic images 'self.synImage' 'self.kImage'
    self.viewSynImage()

    
  def FTcircle(self,kx,ky,a,x0,y0): 
    '''FT of circle of radius a at x0,y0'''
    ka=np.sqrt((kx**2+ky**2))*a
    jka=scipy.special.jn(1,ka)/ka
    ftc=2*np.pi*a*a*np.exp(1j*x0*kx)*np.exp(1j*y0*ky)*jka
    return ftc

  def makeFTcircleArrayImage(self, x0=0.30, y0=0):
    '''Calculates discrete FT of hole array, returns kspace image'''
    npoints=int(self.ofovX/self.pixdx)  #subarray number of rows
    n=int(npoints/2)
    nref=int(self.refImage.shape[0]/2)
    x=range(-n,n)
    y=range(-n,n)
    kx, ky = np.meshgrid(x, y)
    self.kx=2*np.pi*kx/self.ofovX
    self.ky=2*np.pi*ky/self.ofovX
    self.kx[n,n]=1E-9    #cannot have k=0 or divide by zero, make it small 
    self.ky[n,n]=1E-9
    self.kImage=np.zeros((npoints,npoints), dtype=np.complex128)
    for roi in self.currentROIs.ROIs:
      if roi.Type=='Sphere':
        self.kImage+=self.FTcircle(self.kx,self.ky,roi.d1/2,roi.Xcenter,roi.Zcenter) 
    tx=np.exp(1j*x0*self.kx)   #translations
    ty=np.exp(1j*y0*self.ky)
    self.kImage=self.kImage*tx*ty
    self.synImageComplex=np.flip(np.transpose(np.fft.fftshift(np.fft.ifft2(self.kImage))))    #2d fft then shift so zero frequency is in center
    im=np.absolute(self.synImageComplex[n-nref:n+nref,n-nref:n+nref])
    self.synImage=self.nImage(im)

  def translateKImage(self,x0,y0):
    n=int(self.ofovX/self.pixdx/2)  #subarray number of rows
    nref=int(self.refImage.shape[0]/2)
    tx=np.exp(1j*x0*self.kx)   #translations
    ty=np.exp(1j*y0*self.ky)
    kimage=self.kImage*tx*ty
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
        
  def blurFTImages (self):
      '''Blurs image with gaussian point spread function'''  
      nps=100
      self.resArray=np.zeros((nps,2))
      blur=0.005  #blur in increments of 0.1 mm
      m='Calculating image resolution ......'
      self.win.setWindowTitle(m)
      self.messages.message(m)
      k2=2*np.pi**2*(self.kx**2+self.ky**2)
      for i in range(0,nps):
        sigma=i*blur
        filter=np.exp(-sigma**2*k2)
        im=self.makeRealImage(self.kImage, filter)
        dif = self.difIm(im, self.refImage)
        if i==0:
          dif0=dif
        if dif<dif0:
          self.synImage=im
          self.resolution=sigma
          dif0=dif
        self.resArray[i,0]=sigma
        self.resArray[i,1]=dif
      self.resPlot=pg.plot(title='Resolution')
      self.resPlot.plot(self.resArray[:,0],self.resArray[:,1])        
      self.viewSynImage(res=self.resolution)
                                   
  def openFile (self):
    fileName = QtGui.QFileDialog.getOpenFileName(None,"Open Image Files  or DICOMDIR") 
    if not fileName:  #if cancel is pressed return
      return None   
    print('filename= ' + str(fileName))
    fstatus=self.ds.addFile(fileName)
    if fstatus[0]:
        print('success')
        self.plotImage(1)
    else:
        print(fstatus[1])   #Could not open or read file
        
 
  def saveRIFile (self):
    f = QtGui.QFileDialog.getSaveFileName(None,"Save RI file") 
    if not f:  #if cancel is pressed return
      return None   
    f= open(f, 'w')
    f.write('Standard Resolution inset: x,y,d in mm\n')
    for hole in self.holearray:
      s = "{:.4f}".format(hole[0])+ ',' + "{:.4f}".format(hole[1])+ ','+ "{:.4f}".format(hole[2])+ '\n'
      f.write(s)
    f.close()

  def imageMask(self,holes):
    n = 1400
    m=5000
    scale=100   #pixels/mm
    img = np.ones((n,m))
    mask=img==0
    a, b = n/2,m/2
    y,x = np.ogrid[-a:n-a, -b:m-b] 

    #mask = (x-midx)**2 + (y-midy)**2 <= r^2
    for point in holes:
        point=point*scale
        point_mask = (x-point[0])**2 + (y-point[1])**2 <= point[2]*point[2]/4
        mask = np.logical_or(mask, point_mask)
    
    img[mask] = 0
    #plt.imshow(img)
    #plt.title('Ideal Mask')
    #plt.show()       
    

      
if __name__ == '__main__':
    holearray=createResArray()
    app = QtGui.QApplication(sys.argv)
    main = resolutionWindow()
    main.plotOverlay(holearray)
    main.win.show() 
    main.imageMask(holearray)
    sys.exit(app.exec_())    


