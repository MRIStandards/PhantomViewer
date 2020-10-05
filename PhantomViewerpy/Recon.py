# -*- coding: utf-8 -*-
"""
Created on Tue Jun 03 10:56:26 2014
Class to reconstruct and manipulate phantom data
@author: Hannah Erdevig

  execute   "designer\pyuic4 designer\ReconGui.ui -o PhantomViewer\ReconGui.py" from system shell to regenerate ReconGuiGui.py from ReconGuiGui.ui
"""
import sys
import os    #operating system file/directory names
from PyQt4 import QtGui, QtCore
from ReconGui import Ui_ReconGUI  # GUI module
from ImageList import ImageList     # file import and export helper module
import numpy as np
import scipy
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import pyqtgraph.functions as fn

class Recon(QtGui.QMainWindow):
    def __init__(self , parent = None):
        super(Recon, self).__init__()
        pg.setConfigOption('background', 0.2)   #Background on plots 0 = black, 1 = white
        pg.setConfigOption('foreground', 'w')
        self.ui = Ui_ReconGUI()
        self.ui.setupUi(self)
        self.dataSetIsNew = False
    #window 1
        self.imv1 = self.ui.widget_k1
        #self.imv1.getView().setLabel('bottom',"H","mm")    # labels that keep the window from sizing properly
        #self.imv1.getView().setLabel('left',"V","mm")
#        self.imv1.ui.normBtn.hide()
        self.imv1.ui.roiBtn.setText("Line scan")
        self.imv1.vLine = pg.InfiniteLine(pos=None, angle=90, pen=None, movable=False, bounds=None)   #cross hairs
        self.imv1.hLine = pg.InfiniteLine(pos=None, angle=0, pen=None, movable=False, bounds=None)
        self.imv1.addItem(self.imv1.vLine, ignoreBounds=True)
        self.imv1.addItem(self.imv1.hLine, ignoreBounds=True)
        if self.dataSetIsNew == False:
            self.proxy = pg.SignalProxy(self.imv1.view.scene().sigMouseMoved, rateLimit=60, slot=self.mouseMoved)
    #window 2
        self.imv2 = self.ui.widget_k2
        #self.imv2.getView().setLabel('bottom',"H","mm")
        #self.imv2.getView().setLabel('left',"V","mm")
#        self.imv2.ui.normBtn.hide()
        self.imv2.ui.roiBtn.setText("Line scan")
        self.imv2.vLine = pg.InfiniteLine(pos=None, angle=90, pen=None, movable=False, bounds=None)
        self.imv2.hLine = pg.InfiniteLine(pos=None, angle=0, pen=None, movable=False, bounds=None)
        self.imv2.addItem(self.imv2.vLine, ignoreBounds=True)
        self.imv2.addItem(self.imv2.hLine, ignoreBounds=True)
        if self.dataSetIsNew == False:
            self.proxy2 = pg.SignalProxy(self.imv2.view.scene().sigMouseMoved, rateLimit=60, slot=self.mouseMoved2)
    #window 3
        self.imv3 = self.ui.widget_imag
        #self.imv3.getView().setLabel('bottom',"H","mm")
        #self.imv3.getView().setLabel('left',"V","mm")
#        self.imv3.ui.normBtn.hide()
        self.imv3.ui.roiBtn.setText("Line scan")
        self.imv3.vLine = pg.InfiniteLine(pos=None, angle=90, pen=None, movable=False, bounds=None)
        self.imv3.hLine = pg.InfiniteLine(pos=None, angle=0, pen=None, movable=False, bounds=None)
        self.imv3.addItem(self.imv3.vLine, ignoreBounds=True)
        self.imv3.addItem(self.imv3.hLine, ignoreBounds=True)
        if self.dataSetIsNew == False:
            self.proxy3 = pg.SignalProxy(self.imv3.view.scene().sigMouseMoved, rateLimit=60, slot=self.mouseMoved3)
    #window 4
        self.imv4=self.ui.widget_iphase
        #self.imv4.getView().setLabel('bottom',"H","mm")
        #self.imv4.getView().setLabel('left',"V","mm")
#        self.imv4.ui.normBtn.hide()
        self.imv4.ui.roiBtn.setText("Line scan")
        self.imv4.vLine = pg.InfiniteLine(pos=None, angle=90, pen=None, movable=False, bounds=None)
        self.imv4.hLine = pg.InfiniteLine(pos=None, angle=0, pen=None, movable=False, bounds=None)
        self.imv4.addItem(self.imv4.vLine, ignoreBounds=True)
        self.imv4.addItem(self.imv4.hLine, ignoreBounds=True)
        if self.dataSetIsNew == False:
            self.proxy4 = pg.SignalProxy(self.imv4.view.scene().sigMouseMoved, rateLimit=60, slot=self.mouseMoved4)
        
        self.nImages = 0
        self.nCurrentImage = 0
        self.dicomHeader = "DICOM Header"
        self.ui.lineEdit_nimages.setText((str(self.nImages)))
        self.ui.label.setText("none")
        self.dsRe = ImageList()    # Use ImageList.py to create list of image data sets
        self.dsIm = ImageList()
        self.dsMg = ImageList()
        self.dsPh = ImageList()
        self.dsOriginalComplex = ImageList()
        self.dsComplex = ImageList()
        self.dsComplexImage = ImageList()
        self.dsImageMag = ImageList()
        self.dsImagePhase = ImageList()
        self.seriesFileNames = []
        self.windows = [0,0,0,0]
        self.dataSet = [0,0,0,0]
    #signals and slots
#         self.ui.actionNew.triggered.connect(self.NewFile)
        self.ui.actionOpenRI.triggered.connect(self.OpenFileReIm)
        self.ui.actionOpenMP.triggered.connect(self.OpenFileMgPh)
        self.ui.actionSave12.triggered.connect(self.writeDicomFiles12)
        self.ui.actionSave34.triggered.connect(self.writeDicomFiles34)
        self.ui.actionSave_to_Numpy_array.triggered.connect(self.writeNumpyFiles12)
        self.ui.actionClear.triggered.connect(self.ClearImages)
        self.ui.actionDeleteCurrent.triggered.connect(self.deleteCurrentImage)
        self.ui.verticalSlider_slice.valueChanged.connect(self.ImageSlider)
        self.ui.radioButton_ri.clicked.connect(self.SwitchDisplaytoRI)#(self.dsRe, self.dsIm))
        self.ui.radioButton_mp.clicked.connect(self.SwitchDisplaytoMP)#(self.dsMg, self.dsPh))
        self.ui.pushButton_reconstructI.clicked.connect(self.ReconstructImageData)
        self.ui.pushButton_reconstructK.clicked.connect(self.ReconstructRawData)
        self.ui.pushButton_reset.clicked.connect(self.ResetData)
        self.ui.pushButton_apply.clicked.connect(self.EditData)
        
#     def NewFile (self):
#         self.dataSetIsNew = True
#         self.dsRe = ImageList()    # Use ImageList.py to create list of image data sets
#         self.dsIm = ImageList()
#         self.dsMg = ImageList()
#         self.dsPh = ImageList()
#         self.dsOriginalComplex = ImageList()
#         self.dsComplex = ImageList()
#         self.dsComplexImage = ImageList()
#         self.dsImageMag = ImageList()
#         self.dsImagePhase = ImageList()
#         self.dsMg.PA.append(np.zeros([256,256]))
#         self.dsPh.PA.append(np.zeros([256,256]))
#         self.dsRe.PA.append(np.zeros([256,256]))
#         self.dsIm.PA.append(np.zeros([256,256]))
#         self.dsImageMag.PA.append(np.zeros([256,256]))
#         self.dsImagePhase.PA.append(np.zeros([256,256]))
#         self.dataSet = self.dsMg, self.dsPh, self.dsImageMag, self.dsImagePhase
#         self.windows = [1,1,1,1]
#         self.nImages = 1
#         self.ui.lineEdit_nimages.setText(str(self.nImages))
#         self.ui.verticalSlider_slice.setMinimum(1)       #set slider to go from 1 to the number of images
#         self.ui.verticalSlider_slice.setMaximum(self.nImages)
#         self.nCurrentImage = 1
#         self.ui.verticalSlider_slice.setValue(self.nCurrentImage)
#         self.DisplayCurrentImage(self.imv1, self.dsMg)
#         self.DisplayCurrentImage(self.imv2, self.dsPh)
#         self.DisplayCurrentImage(self.imv3, self.dsImageMag)
#         self.DisplayCurrentImage(self.imv4, self.dsImagePhase)
#         self.ui.label_k1.setText("K-Space [Magnitude]")
#         self.ui.label_k2.setText("K-Space [Phase]")
         
    def OpenFileReIm (self):
        self.dataSetIsNew = False
        self.dsRe = ImageList()    # Use ImageList.py to create list of image data sets
        self.dsIm = ImageList()
        self.dsMg = ImageList()
        self.dsPh = ImageList()
        self.dsOriginalComplex = ImageList()
        self.dsComplex = ImageList()
        self.dsComplexImage = ImageList()
        self.dsImageMag = ImageList()
        self.dsImagePhase = ImageList()
        
        #REAL FILES
        self.fileNamesRe = QtGui.QFileDialog.getOpenFileNames(self,"Open Real Image Files", "/home/file", "Image Files (*.dcm *.DCM *.bmp *.tif *.fdf)")
        if not self.fileNamesRe:  #if cancel is pressed return
            return None        
        self.seriesFileNames.extend(self.fileNamesRe)     #concatenate new file list with previous file list
        for i in range(len(self.fileNamesRe)):
            fileName = self.fileNamesRe[i]
            self.dsRe.addFile(fileName)
            self.dsOriginalComplex.addFile(fileName)
            self.dsComplex.addFile(fileName)
            self.dsComplexImage.addFile(fileName)
            self.dsImageMag.addFile(fileName)
            self.dsImagePhase.addFile(fileName)
            self.dsMg.addFile(fileName)
            self.dsPh.addFile(fileName)
        self.windows[0] = 1
        self.dataSet[0] = self.dsRe
        self.nImages=self.nImages+len(self.fileNamesRe)
        self.ui.lineEdit_nimages.setText(str(self.nImages))
        self.ui.verticalSlider_slice.setMinimum(1)       #set slider to go from 1 to the number of images
        self.ui.verticalSlider_slice.setMaximum(self.nImages)
        self.nCurrentImage=1
        self.ui.verticalSlider_slice.setValue(self.nCurrentImage)
        self.DisplayCurrentImage(self.imv1, self.dsRe)
        self.ui.label_k1.setText("K-Space [Real]")
        
        #IMAGINARY FILES
        self.fileNamesIm = QtGui.QFileDialog.getOpenFileNames(self,"Open Imaginary Image Files", "/home/file", "Image Files (*.dcm *.DCM *.bmp *.tif *.fdf)")
        if not self.fileNamesIm:  #if cancel is pressed return
            return None        
        self.seriesFileNames.extend(self.fileNamesIm)     #concatenate new file list with previous file list
        for i in range(len(self.fileNamesIm)):
            fileName = self.fileNamesIm[i]
            self.dsIm.addFile(fileName)
        self.windows[1] = 1
        self.dataSet[1] = self.dsIm
        self.ui.lineEdit_nimages.setText(str(self.nImages))
        self.ui.verticalSlider_slice.setMinimum(1)       #set slider to go from 1 to the number of images
        self.ui.verticalSlider_slice.setMaximum(self.nImages)
        self.nCurrentImage=1
        self.ui.verticalSlider_slice.setValue(self.nCurrentImage)
        self.DisplayCurrentImage(self.imv2,self.dsIm)
        self.ui.label_k2.setText("K-Space [Imaginary]")
        
        #CREATE COMPLEX ARRAY
#         self.dsMg = self.dataSet[0]
#         self.dsPh = self.dataSet[0]
#         self.dsImageMag = self.dataSet[0]
#         self.dsImagePhase = self.dataSet[0]
        for i in range(1, len(self.dsMg.PA)):
            self.dsOriginalComplex.PA[i] = self.dsRe.PA[i] + 1j*self.dsIm.PA[i]
            self.dsComplex.PA[i] = self.dsRe.PA[i] + 1j*self.dsIm.PA[i]
            self.dsMg.PA[i] = np.absolute(self.dsComplex.PA[i])
            self.dsPh.PA[i] = np.angle(self.dsComplex.PA[i])
            self.dsComplexImage.PA[i] = np.fft.fft2(self.dsComplex.PA[i])
            self.dsImageMag.PA[i] = np.absolute(self.dsComplexImage.PA[i])
            self.dsImagePhase.PA[i] = np.angle(self.dsComplexImage.PA[i])
        
    def OpenFileMgPh (self):
        self.dataSetIsNew = False
        self.dsRe = ImageList()    # Use ImageList.py to create list of image data sets
        self.dsIm = ImageList()
        self.dsMg = ImageList()
        self.dsPh = ImageList()
        self.dsOriginalComplex = ImageList()
        self.dsComplex = ImageList()
        self.dsComplexImage = ImageList()
        self.dsImageMag = ImageList()
        self.dsImagePhase = ImageList()
        
        #MAGNITUDE FILES
        self.fileNamesMg = QtGui.QFileDialog.getOpenFileNames(self,"Open Magnitude Image Files", "/home/file", "Image Files (*.dcm *.DCM *.bmp *.tif *.fdf)")
        if not self.fileNamesMg:  #if cancel is pressed return
            return None        
        self.seriesFileNames.extend(self.fileNamesMg)     #concatenate new file list with previous file list
        d1Mg= []  #d1 is 3d data stack for 3d images
        for i in range(len(self.fileNamesMg)):
            fileName = self.fileNamesMg[i]
            self.dsMg.addFile(fileName)
            self.dsOriginalComplex.addFile(fileName)
            self.dsComplex.addFile(fileName)
            self.dsComplexImage.addFile(fileName)
            self.dsImageMag.addFile(fileName)
            self.dsImagePhase.addFile(fileName)
            self.dsRe.addFile(fileName)
            self.dsIm.addFile(fileName)
            d1Mg.append(self.dsMg.PA[i+1])
        self.windows[0] = 1
        self.dataSet[0] = self.dsMg
        self.nImages=self.nImages+len(self.fileNamesMg)
        self.ui.lineEdit_nimages.setText(str(self.nImages))
        self.ui.verticalSlider_slice.setMinimum(1)       #set slider to go from 1 to the number of images
        self.ui.verticalSlider_slice.setMaximum(self.nImages)
        self.nCurrentImage=1
        self.ui.verticalSlider_slice.setValue(self.nCurrentImage)
        self.DisplayCurrentImage(self.imv1,self.dsMg)
        self.ui.label_k1.setText("K-Space [Magnitude]")
        self.image3D = np.dstack(d1Mg)
#         self.msgPrint("image size" + str(d1Mg[0].shape) + "; image3D size" + str(self.image3D.shape)+ os.linesep)
        
        #PHASE FILES
        self.fileNamesPh = QtGui.QFileDialog.getOpenFileNames(self,"Open Phase Image Files", "/home/file", "Image Files (*.dcm *.DCM *.bmp *.tif *.fdf)")
        extension = self.fileNamesPh[0].split(".")[-1]       
        if not self.fileNamesPh:  #if cancel is pressed return
            return None
        self.seriesFileNames.extend(self.fileNamesPh)     #concatenate new file list with previous file list
        for i in range(len(self.fileNamesPh)):
            fileName = self.fileNamesPh[i]
            self.dsPh.addFile(fileName)
            if extension.toLower() == "dcm" or extension==fileName or extension.toLower() == "ima": # if .dcm , .ima or  no extension try dicom
                self.dsPh.PA[i+1]= self.dsPh.PA[i+1].astype(float)/10000.0 - np.pi       #DICOM stores phase in integer*10000, convert to float
        self.windows[1] = 1
        self.dataSet[1] = self.dsPh
        self.ui.lineEdit_nimages.setText(str(self.nImages))
        self.ui.verticalSlider_slice.setMinimum(1)       #set slider to go from 1 to the number of images
        self.ui.verticalSlider_slice.setMaximum(self.nImages)
        self.nCurrentImage=1    # 0th item is dummy image
        self.ui.verticalSlider_slice.setValue(self.nCurrentImage)
        self.DisplayCurrentImage(self.imv2, self.dsPh)
        self.ui.label_k2.setText("K-Space [Phase]")
        
        #CREATE COMPLEX ARRAY
        for i in range(1, len(self.dsMg.PA)):
            self.dsOriginalComplex.PA[i] = np.multiply(self.dsMg.PA[i],np.exp(1j*self.dsPh.PA[i]))
            self.dsComplex.PA[i] = np.multiply(self.dsMg.PA[i],np.exp(1j*self.dsPh.PA[i]))
            self.dsRe.PA[i] = self.dsComplex.PA[i].real
            self.dsIm.PA[i] = self.dsComplex.PA[i].imag
            self.dsComplexImage.PA[i] = np.fft.fft2(self.dsComplex.PA[i])
            self.dsComplexImage.PA[i] = np.fft.fftshift(self.dsComplexImage.PA[i])
            self.dsImageMag.PA[i] = np.absolute(self.dsComplexImage.PA[i])
            self.dsImagePhase.PA[i] = np.angle(self.dsComplexImage.PA[i])     

    def writeNumpyFiles12 (self):
        fileName = QtGui.QFileDialog.getSaveFileName(parent=None, caption="Numpy File Name")
        if not fileName:  #if cancel is pressed return
            return None
        np.save(str(fileName), self.dsComplex.PA[1])   #write current image list in DICOM format to filename+ imagenumber + .npy
        file = open(str(fileName)+'.txt','w') 
        file.write(self.dsComplex.header[1]) 
        file.close() 
        

    def writeDicomFiles12 (self):
        fileName = QtGui.QFileDialog.getSaveFileName(parent=None, caption="Dicom File Name")
        if not fileName:  #if cancel is pressed return
            return None
        self.dataSet[0].writeDicomFiles(str(fileName) + "kmag")   #write current image list in DICOM format to filename+ imagenumber + .dcm
        self.dataSet[1].writeDicomFiles(str(fileName) + "kphase")   #write current image list in DICOM format to filename+ imagenumber + .dcm
        
    def writeDicomFiles34 (self):
        fileName = QtGui.QFileDialog.getSaveFileName(parent=None, caption="Dicom File Name")
        if not fileName:  #if cancel is pressed return
            return None
        self.dataSet[2].writeDicomFiles(str(fileName) + "rmag")   #write current image list in DICOM format to filename+ imagenumber + .dcm
        self.dataSet[3].writeDicomFiles(str(fileName) + "rphase")   #write current image list in DICOM format to filename+ imagenumber + .dcm
            
    def EditData (self):
        if self.ui.radioButton_setValue.isChecked():
            val = float(self.ui.lineEdit_setValue.text())
            if self.ui.groupBox_voxelSelection.isChecked():             #if voxel selection is selected
                rstart, rend = int(self.ui.lineEdit_voxelRowStart.text()), int(self.ui.lineEdit_voxelRowEnd.text())
                cstart, cend = int(self.ui.lineEdit_voxelColStart.text()), int(self.ui.lineEdit_voxelColEnd.text())
                if self.ui.checkBox_w1.isChecked():                     # only apply editing to window if corresponding checkBox is checked
                    for i in range(1, len(self.dsMg.PA)):               # for each image in the stack (same for all types)
                        for c in range(cstart-1, cend):                 # user inputs pixels numbered 1:imagesize
                            for r in range(rstart-1, rend):
                                self.dataSet[0].PA[i][c, r] = val
                    self.DisplayCurrentImage(self.imv1,self.dataSet[0])
                if self.ui.checkBox_w2.isChecked():
                    for i in range(1, len(self.dsMg.PA)):               # for each image in the stack (same for all types)
                        for c in range(cstart-1, cend):                 # user inputs pixels numbered 1:imagesize
                            for r in range(rstart-1, rend):
                                self.dataSet[1].PA[i][c, r] = val
                    self.DisplayCurrentImage(self.imv2,self.dataSet[1])
                    
            if self.ui.groupBox_regionSelection.isChecked():            #if region selection is selected
                xCenter, yCenter, radius = int(self.ui.lineEdit_Xpos.text()), int(self.ui.lineEdit_Ypos.text()), int(self.ui.lineEdit_radius.text())
                # selects a circular region in coordinate system where x and y range from [-ImageSize/2, ImageSize/2]
                if self.ui.checkBox_w1.isChecked():                     # only apply editing to window if corresponding checkBox is checked
                    for i in range(1, len(self.dataSet[0].PA)):               # iterate through arrays                                
                        if (self.ui.radioButton_interior.isChecked()):     # voxel is inside circle
                            lx, ly = self.dataSet[0].PA[i].shape
                            X, Y = np.ogrid[0:lx, 0:ly]
                            mask = (X - lx/2-xCenter)**2 + (Y - ly/2+yCenter)**2 <= radius**2
                            self.dataSet[0].PA[i][mask] = val
                        elif (self.ui.radioButton_exterior.isChecked()):    # voxel is outside circle
                            lx, ly = self.dataSet[0].PA[i].shape
                            X, Y = np.ogrid[0:lx, 0:ly]
                            mask = (X - lx/2-xCenter)**2 + (Y - ly/2+yCenter)**2 > radius**2
                            self.dataSet[0].PA[i][mask] = val
                    self.DisplayCurrentImage(self.imv1,self.dataSet[0])
                    
                if self.ui.checkBox_w2.isChecked():                     # only apply editing to window if corresponding checkBox is checked
                    for i in range(1, len(self.dataSet[1].PA)):               # iterate through arrays                                
                        if (self.ui.radioButton_interior.isChecked()):     # voxel is inside circle
                            lx, ly = self.dataSet[1].PA[i].shape
                            X, Y = np.ogrid[0:lx, 0:ly]
                            mask = (X - lx/2-xCenter)**2 + (Y - ly/2+yCenter)**2 <= radius**2
                            self.dataSet[1].PA[i][mask] = val
                        elif (self.ui.radioButton_exterior.isChecked()):    # voxel is outside circle
                            lx, ly = self.dataSet[1].PA[i].shape
                            X, Y = np.ogrid[0:lx, 0:ly]
                            mask = (X - lx/2-xCenter)**2 + (Y - ly/2+yCenter)**2 > radius**2
                            self.dataSet[1].PA[i][mask] = val
                    self.DisplayCurrentImage(self.imv2,self.dataSet[1])
                
        if self.ui.radioButton_addValue.isChecked():
            val = float(self.ui.lineEdit_addValue.text())
            if self.ui.groupBox_voxelSelection.isChecked():             #if voxel selection is selected
                rstart, rend = int(self.ui.lineEdit_voxelRowStart.text()), int(self.ui.lineEdit_voxelRowEnd.text())
                cstart, cend = int(self.ui.lineEdit_voxelColStart.text()), int(self.ui.lineEdit_voxelColEnd.text())
                if self.ui.checkBox_w1.isChecked():                     # only apply editing to window if corresponding checkBox is checked
                    for i in range(1, len(self.dsMg.PA)):               # for each image in the stack (same for all types)
                        for c in range(cstart-1, cend):                 # user inputs pixels numbered 1:imagesize
                            for r in range(rstart-1, rend):
                                self.dataSet[0].PA[i][c, r] += val
                    self.DisplayCurrentImage(self.imv1,self.dataSet[0])
                if self.ui.checkBox_w2.isChecked():
                    for i in range(1, len(self.dsMg.PA)):               # for each image in the stack (same for all types)
                        for c in range(cstart-1, cend):                 # user inputs pixels numbered 1:imagesize
                            for r in range(rstart-1, rend):
                                self.dataSet[1].PA[i][c, r] += val
                    self.DisplayCurrentImage(self.imv2,self.dataSet[1])
                    
            if self.ui.groupBox_regionSelection.isChecked():            #if region selection is selected
                xCenter, yCenter, radius = int(self.ui.lineEdit_Xpos.text()), int(self.ui.lineEdit_Ypos.text()), int(self.ui.lineEdit_radius.text())
                # selects a circular region in coordinate system where x and y range from [-ImageSize/2, ImageSize/2]
                if self.ui.checkBox_w1.isChecked():                     # only apply editing to window if corresponding checkBox is checked
                    for i in range(1, len(self.dataSet[0].PA)):               # iterate through arrays                                
                        if (self.ui.radioButton_interior.isChecked()):     # voxel is inside circle
                            lx, ly = self.dataSet[0].PA[i].shape
                            X, Y = np.ogrid[0:lx, 0:ly]
                            mask = (X - lx/2-xCenter)**2 + (Y - ly/2+yCenter)**2 <= radius**2
                            self.dataSet[0].PA[i][mask] += val
                        elif (self.ui.radioButton_exterior.isChecked()):    # voxel is outside circle
                            lx, ly = self.dataSet[0].PA[i].shape
                            X, Y = np.ogrid[0:lx, 0:ly]
                            mask = (X - lx/2-xCenter)**2 + (Y - ly/2+yCenter)**2 > radius**2
                            self.dataSet[0].PA[i][mask] += val
                    self.DisplayCurrentImage(self.imv1,self.dataSet[0])
                    
                if self.ui.checkBox_w2.isChecked():                     # only apply editing to window if corresponding checkBox is checked
                    for i in range(1, len(self.dataSet[1].PA)):               # iterate through arrays                                
                        if (self.ui.radioButton_interior.isChecked()):     # voxel is inside circle
                            lx, ly = self.dataSet[1].PA[i].shape
                            X, Y = np.ogrid[0:lx, 0:ly]
                            mask = (X - lx/2-xCenter)**2 + (Y - ly/2+yCenter)**2 <= radius**2
                            self.dataSet[1].PA[i][mask] += val
                        elif (self.ui.radioButton_exterior.isChecked()):    # voxel is outside circle
                            lx, ly = self.dataSet[1].PA[i].shape
                            X, Y = np.ogrid[0:lx, 0:ly]
                            mask = (X - lx/2-xCenter)**2 + (Y - ly/2+yCenter)**2 > radius**2
                            self.dataSet[1].PA[i][mask] += val
                    self.DisplayCurrentImage(self.imv2,self.dataSet[1])
            
        #RECREATE COMPLEX ARRAY
        if self.dataSet[0] == self.dsMg:
            for i in range(1, len(self.dsMg.PA)):
                self.dsComplex.PA[i] = np.multiply(self.dsMg.PA[i],np.exp(1j*self.dsPh.PA[i]))
                self.dsRe.PA[i] = self.dsComplex.PA[i].real
                self.dsIm.PA[i] = self.dsComplex.PA[i].imag
                self.dsComplexImage.PA[i] = np.fft.fft2(self.dsComplex.PA[i])
                self.dsComplexImage.PA[i] = np.fft.fftshift(self.dsComplexImage.PA[i])
                self.dsImageMag.PA[i] = np.absolute(self.dsComplexImage.PA[i])
                self.dsImagePhase.PA[i] = np.angle(self.dsComplexImage.PA[i])
                
        if self.dataSet[0] == self.dsRe:
            for i in range(1, len(self.dsMg.PA)):
                self.dsComplex.PA[i] = self.dsRe.PA[i] + 1j*self.dsIm.PA[i]
                self.dsMg.PA[i] = np.absolute(self.dsComplex.PA[i])
                self.dsPh.PA[i] = np.angle(self.dsComplex.PA[i])
                self.dsComplexImage.PA[i] = np.fft.fft2(self.dsComplex.PA[i])
                self.dsComplexImage.PA[i] = np.fft.fftshift(self.dsComplexImage.PA[i])
                self.dsImageMag.PA[i] = np.absolute(self.dsComplexImage.PA[i])
                self.dsImagePhase.PA[i] = np.angle(self.dsComplexImage.PA[i])
                              
    def ResetData (self):
        for i in range(1, len(self.dsMg.PA)):
            self.dsRe.PA[i] = self.dsOriginalComplex.PA[i].real
            self.dsIm.PA[i] = self.dsOriginalComplex.PA[i].imag
            self.dsMg.PA[i] = np.absolute(self.dsOriginalComplex.PA[i])
            self.dsPh.PA[i] = np.angle(self.dsOriginalComplex.PA[i])
            self.dsComplexImage.PA[i] = np.fft.fft2(self.dsOriginalComplex.PA[i])
            self.dsComplexImage.PA[i] = np.fft.fftshift(self.dsComplexImage.PA[i])
            self.dsImageMag.PA[i] = np.absolute(self.dsComplexImage.PA[i])
            self.dsImagePhase.PA[i] = np.angle(self.dsComplexImage.PA[i])
        if self.ui.radioButton_mp.isChecked():
            self.dataSet[0] = self.dsRe
            self.dataSet[1] = self.dsIm
        if self.ui.radioButton_ri.isChecked():
            self.dataSet[0] = self.dsMg
            self.dataSet[1] = self.dsPh
        self.dataSet[2] = self.dsImageMag
        self.dataSet[3] = self.dsImagePhase
        self.DisplayCurrentImage(self.imv1, self.dataSet[0])
        self.DisplayCurrentImage(self.imv2, self.dataSet[1])
        if self.windows[2] == 1:
            self.DisplayCurrentImage(self.imv3, self.dataSet[2])
        if self.windows[3] == 1:
            self.DisplayCurrentImage(self.imv4, self.dataSet[3])

    def ReconstructImageData (self):
        self.windows[2] = 1
        self.windows[3] = 1
        self.dataSet[2] = self.dsImageMag
        self.dataSet[3] = self.dsImagePhase
        self.DisplayCurrentImage(self.imv3, self.dataSet[2])
        self.DisplayCurrentImage(self.imv4, self.dataSet[3])
        
    def ReconstructRawData (self):
        for i in range(1, len(type.PA)):
            self.dsComplexImage.PA[i] = np.multiply(self.dsImageMag.PA[i],np.exp(1j*self.dsImagePhase.PA[i]))
            self.dsComplexImage.PA[i] = np.fft.ifft2(self.dsComplex.PA[i])
            self.dsComplexImage.PA[i] = np.fft.fftshift(self.dsComplexImage.PA[i])
            self.dsRe.PA[i] = self.dsComplex.PA[i].real
            self.dsIm.PA[i] = self.dsComplex.PA[i].imag
        self.windows[2] = 1
        self.windows[3] = 1
        self.dataSet[2] = self.dsImageMag
        self.dataSet[3] = self.dsImagePhase
        self.DisplayCurrentImage(self.imv3, self.dataSet[2])
        self.DisplayCurrentImage(self.imv4, self.dataSet[3])
        
    def SwitchDisplaytoRI (self):#,data1,data2):   # updates windows 1&2 to display current k-space pair
        if (self.windows[0] == 1):
            self.dataSet[0] = self.dsRe
            self.DisplayCurrentImage(self.imv1, self.dsRe)#data1)
            self.ui.label_k1.setText("K-Space [Real]")
        if (self.windows[1] == 1):
            self.dataSet[1] = self.dsIm
            self.DisplayCurrentImage(self.imv2, self.dsIm)#data2
            self.ui.label_k2.setText("K-Space [Imaginary]")
            
    def SwitchDisplaytoMP (self):#,data1,data2):   # updates windows 1&2 to display current k-space pair
        if (self.windows[0] == 1):
            self.dataSet[0] = self.dsMg
            self.DisplayCurrentImage(self.imv1, self.dsMg)#data1)
            self.ui.label_k1.setText("K-Space [Magnitude]")
        if (self.windows[1] == 1):
            self.dataSet[1] = self.dsPh
            self.DisplayCurrentImage(self.imv2, self.dsPh)#data2
            self.ui.label_k2.setText("K-Space [Phase]")
        
    def ImageSlider (self):
        self.nCurrentImage = self.ui.verticalSlider_slice.value()
        if (self.windows[0] == 1):
            self.DisplayCurrentImage(self.imv1, self.dataSet[0])
        if (self.windows[1] == 1):
            self.DisplayCurrentImage(self.imv2, self.dataSet[1])
        if (self.windows[2] == 1):
            self.DisplayCurrentImage(self.imv3, self.dataSet[2])
        if (self.windows[3] == 1):
            self.DisplayCurrentImage(self.imv4, self.dataSet[3])
        
    def DisplayCurrentImage (self,win,dstype):  # display the current image from the data "dstype" in window "win"
        i = self.nCurrentImage
        self.ui.label.setText(str(self.nCurrentImage))
        if self.dataSetIsNew == False:
            self.ui.lineEdit_date.setText(format(dstype.StudyDate[i]))
            self.ui.textEdit_filename.setText(self.seriesFileNames[i-1]) if i > 0 else self.ui.textEdit_filename.setText("")
            self.ui.lineEdit_manufacturer.setText(dstype.Manufacturer[i]) 
            self.ui.lineEdit_series.setText(dstype.SeriesDescription[i]) 
            self.ui.lineEdit_institution.setText(dstype.InstitutionName[i]) 
            self.ui.lineEdit_fieldT.setText(str(dstype.MagneticFieldStrength[i]))
            self.ui.lineEdit_receivecoil.setText(str(dstype.ReceiveCoilName[i]))    
            #self.ui.lblPatient.setText(dstype.PatientName[i]) 
            self.ui.lineEdit_protocol.setText(str(dstype.ProtocolName[i])) 
            self.ui.lineEdit_bandwidth.setText(str(dstype.PixelBandwidth[i])) 
            self.ui.lineEdit_TE.setStyleSheet("background-color: white") if (self.checkEqual(dstype.TE)) else self.ui.lineEdit_TE.setStyleSheet("background-color: yellow")  
            self.ui.lineEdit_TE.setText(str(dstype.TE[i]))
            self.ui.lineEdit_TR.setStyleSheet("background-color: white") if (self.checkEqual(dstype.TR)) else self.ui.lineEdit_TR.setStyleSheet("background-color: yellow")
            self.ui.lineEdit_TR.setText(str(dstype.TR[i])) 
            self.ui.lineEdit_imagesize_col.setText(str(dstype.Columns[i]))  
            self.ui.lineEdit_imagesize_row.setText(str(dstype.Rows[i]))
            self.ui.lineEdit_TI.setStyleSheet("background-color: white") if (self.checkEqual(dstype.TI)) else self.ui.lineEdit_TI.setStyleSheet("background-color: yellow")    
            self.ui.lineEdit_TI.setText(str(dstype.TI[i])) 
            self.ui.lineEdit_slicethick.setText(str(dstype.SliceThickness[i]))
            self.ui.lineEdit_sliceloc.setStyleSheet("background-color: white") if (self.checkEqual(dstype.SliceLocation)) else self.ui.lineEdit_sliceloc.setStyleSheet("background-color: yellow")
            self.ui.lineEdit_sliceloc.setText(str(dstype.SliceLocation[i])) 
            self.ui.lineEdit_pixsize_row.setText(str(dstype.PixelSpacingX[i])) 
            self.ui.lineEdit_pixsize_col.setText(str(dstype.PixelSpacingY[i])) 
    #         self.ui.lblFA.setStyleSheet("background-color: white") if (self.checkEqual(dstype.FA)) else self.ui.lblFA.setStyleSheet("background-color: yellow")    
    #         self.ui.lblFA.setText(str(dstype.FA[i]))
            self.ui.lineEdit_phasedir.setText(str(dstype.InPlanePhaseEncodingDirection[i])) 
            self.ui.lineEdit_FoVX.setText(str(dstype.FoVX[i]))
            self.ui.lineEdit_FoVY.setText(str(dstype.FoVY[i]))    
            self.ui.lineEdit_b.setText(str(dstype.bValue[i]))
            self.ui.textEdit_header.setText(dstype.header[i])
        data = dstype.PA[i]  #not sure why we need to transpose]
#         xscale = dstype.PixelSpacingX[i] if (dstype.PixelSpacingX[i] > 0.) else 1
#         yscale = dstype.PixelSpacingY[i] if (dstype.PixelSpacingY[i] > 0.) else 1
#         xmin = -dstype.FoVX[i]/2   #set origin to center of image, need to upgrade to set by DICOM tag
#         ymin = -dstype.FoVY[i]/2
#         textEdit_results was lblUpperLeft in the line below
        #self.ui.textEdit_results.setText("UL=" + "{:.1f}".format(dstype.ImagePosition[i][0]) + "," + "{:.1f}".format(dstype.ImagePosition[i][1]) + "," + "{:.1f}".format(dstype.ImagePosition[i][2]))
#         setImage(img, autoRange=True, autoLevels=True, levels=None, axes=None, xvals=None, pos=None, scale=None, transform=None, autoHistogramRange=True)
        win.setImage(data, autoRange=True, autoLevels=True, autoHistogramRange=True)#pos = (xmin,ymin), scale = (xscale,yscale), autoHistogramRange=True)
#         self.imv1.getView().setLabel('bottom',self.DirectionLabel(dstype.RowDirection[i]),"mm")
#         self.imv1.getView().setLabel('left',self.DirectionLabel(dstype.ColumnDirection[i]),"mm")
 
    def checkEqual(self, lst):    #returns True if all elemments (except the 0th element) of the list are equal
        return lst[2:] == lst[1:-1]  
 
    def ClearImages (self):  #Deletes all images except default image at index 1
        self.ds = ImageList()                         #list of data sets, can be dicom, tiff, fdf
        self.dsRe = ImageList()    # Uses separate module to create list of image data sets
        self.dsIm = ImageList()
        self.dsMg = ImageList()
        self.dsPh = ImageList()
        self.dsComplex = ImageList()
        self.dsOriginalComplex = ImageList()
        self.dsComplexImage = ImageList()
        self.dsImageMag = ImageList()
        self.dsImagePhase = ImageList()
        del self.seriesFileNames[:]
        self.windows = [0,0,0,0]
        self.nCurrentImage=0
        self.nImages=0
#         self.image3D.zeros[1,1,1]
        self.DisplayCurrentImage(self.imv1, self.ds)
        self.DisplayCurrentImage(self.imv2, self.ds)
        self.DisplayCurrentImage(self.imv3, self.ds)
        self.DisplayCurrentImage(self.imv4, self.ds)
        self.ui.lineEdit_nimages.setText(str(self.nImages))
        self.ui.verticalSlider_slice.setMaximum(0)
 
    def deleteCurrentImage(self):
        if self.nCurrentImage > 0:
            self.dsRe.deleteImage(self.nCurrentImage)
            self.dsIm.deleteImage(self.nCurrentImage)
            self.dsMg.deleteImage(self.nCurrentImage)
            self.dsPh.deleteImage(self.nCurrentImage)
            self.dsComplex.deleteImage(self.nCurrentImage)
            self.dsOriginalComplex.deleteImage(self.nCurrentImage)
            self.dsComplexImage.deleteImage(self.nCurrentImage)
            self.dsImageMag.deleteImage(self.nCurrentImage)
            self.dsImagePhase.deleteImage(self.nCurrentImage)
            self.nImages -= 1
            self.ui.lineEdit_nimages.setText(str(self.nImages))
            self.ui.verticalSlider_slice.setMinimum(1)       #set slider to go from 1 to the number of images
            self.ui.verticalSlider_slice.setMaximum(self.nImages)
            if self.nImages == 0:
                self.nCurrentImage=0
                self.windows = [0,0,0,0]
                self.ds = ImageList()
                self.DisplayCurrentImage(self.imv1, self.ds)
            else:
                self.nCurrentImage = 1
            self.ui.verticalSlider_slice.setValue(self.nCurrentImage)
            self.DisplayCurrentImage(self.imv1, self.dataSet[0])
            self.DisplayCurrentImage(self.imv2, self.dataSet[1])
       
    def ViewDicomHeader (self):        
        if self.ui.rbViewDicomHeader.isChecked():
            self.ui.textEdit_header.setHidden(False)
            dh = str(self.ds.header[self.nCurrentImage])
            if dh == '':
                dh="DICOM Header"
            self.ui.textEdit_header.setText(dh)
        else:
            self.ui.textEdit_header.setHidden(True)
 
#     def View3d(self):
#         w = gl.GLViewWidget()
#         w.opts['distance'] = 200
#         w.show()
#         w.setWindowTitle('3D View')
#         g = gl.GLGridItem()
#         g.scale(10, 10, 10)
#         w.addItem(g)   
#         data=self.image3D
#         #positive = np.log(np.clip(data, 0, data.max())**2)
#         #negative = np.log(np.clip(-data, 0, -data.min())**2)
#         d2 = np.empty(data.shape + (4,), dtype=np.ubyte)
#         d2[..., 0] = data * (255./data.max())
#         d2[..., 1] = data * (255./data.max())
#         d2[..., 2] = d2[...,1]
#         d2[..., 3] = d2[..., 0]*0.3 + d2[..., 1]*0.3
#         d2[..., 3] = (d2[..., 3].astype(float) / 255.) **2 * 255
#           
#         d2[:, 0, 0] = [255,0,0,100]
#         d2[0, :, 0] = [0,255,0,100]
#         d2[0, 0, :] = [0,0,255,100]
#           
#         v = gl.GLVolumeItem(d2)
#         v.translate(-128,-128,0)
#         w.addItem(v)        
#         ax = gl.GLAxisItem()
#         w.addItem(ax)
 
    def mouseMoved(self,evt): #mouse move event to move crosshairs and display location and values
        if self.dataSetIsNew == False:
            if (self.windows[0] == 1):
                self.ds = self.dataSet[0]
                pos = evt[0]  ## using signal proxy turns original arguments into a tuple
                if self.imv1.view.sceneBoundingRect().contains(pos):
                    mousePoint = self.imv1.view.mapSceneToView(pos)
                    self.ui.lineEdit_h.setText("{:.2f}".format(mousePoint.x()))
                    self.ui.lineEdit_v.setText("{:.2f}".format(mousePoint.y()))
                    if abs(mousePoint.x()) < self.ds.FoVX[self.nCurrentImage]/2 and abs(mousePoint.y()) < self.ds.FoVY[self.nCurrentImage]/2:
                        Xindex = int((mousePoint.x()+self.ds.FoVX[self.nCurrentImage]/2)/self.ds.PixelSpacingX[self.nCurrentImage]) #if self.ds.PixelSpacingX[self.nCurrentImage] > 0. else Xindex = int(mousePoint.x())
                        Yindex = int((mousePoint.y()+self.ds.FoVY[self.nCurrentImage]/2)/self.ds.PixelSpacingY[self.nCurrentImage]) #if self.ds.PixelSpacingY[self.nCurrentImage] > 0. else Yindex = int(mousePoint.y())
                        value=  self.ds.PA[self.nCurrentImage][Xindex,Yindex]      
                        self.ui.lineEdit_value.setText("{:.1f}".format(value))
                        rc= self.ReltoGlobal(mousePoint.x(), mousePoint.y(), self.nCurrentImage, self.ds)
                        self.ui.lineEdit_x.setText("{:.2f}".format(rc[0]))
                        self.ui.lineEdit_y.setText("{:.2f}".format(rc[1]))
                        self.ui.lineEdit_z.setText("{:.2f}".format(rc[2]))
                    self.imv1.vLine.setPos(mousePoint.x())
                    self.imv1.hLine.setPos(mousePoint.y())
            
    def mouseMoved2(self,evt): #mouse move event to move crosshairs and display location and values
        if self.dataSetIsNew == False:
            if (self.windows[1] == 1):
                self.ds = self.dataSet[1]
                pos = evt[0]  ## using signal proxy turns original arguments into a tuple
                if self.imv2.view.sceneBoundingRect().contains(pos):
                    mousePoint = self.imv2.view.mapSceneToView(pos)
                    self.ui.lineEdit_h.setText("{:.2f}".format(mousePoint.x()))
                    self.ui.lineEdit_v.setText("{:.2f}".format(mousePoint.y()))
                    if abs(mousePoint.x()) < self.ds.FoVX[self.nCurrentImage]/2 and abs(mousePoint.y()) < self.ds.FoVY[self.nCurrentImage]/2:
                        Xindex = int((mousePoint.x()+self.ds.FoVX[self.nCurrentImage]/2)/self.ds.PixelSpacingX[self.nCurrentImage]) #if self.ds.PixelSpacingX[self.nCurrentImage] > 0. else Xindex = int(mousePoint.x())
                        Yindex = int((mousePoint.y()+self.ds.FoVY[self.nCurrentImage]/2)/self.ds.PixelSpacingY[self.nCurrentImage]) #if self.ds.PixelSpacingY[self.nCurrentImage] > 0. else Yindex = int(mousePoint.y())
                        value=  self.ds.PA[self.nCurrentImage][Xindex,Yindex]      
                        self.ui.lineEdit_value.setText("{:.1f}".format(value))
                        rc= self.ReltoGlobal(mousePoint.x(), mousePoint.y(), self.nCurrentImage, self.ds)
                        self.ui.lineEdit_x.setText("{:.2f}".format(rc[0]))
                        self.ui.lineEdit_y.setText("{:.2f}".format(rc[1]))
                        self.ui.lineEdit_z.setText("{:.2f}".format(rc[2]))
                    self.imv2.vLine.setPos(mousePoint.x())
                    self.imv2.hLine.setPos(mousePoint.y()) 
            
    def mouseMoved3(self,evt): #mouse move event to move crosshairs and display location and values
        if self.dataSetIsNew == False:
            if (self.windows[2] == 1):
                self.ds = self.dataSet[2]
                pos = evt[0]  ## using signal proxy turns original arguments into a tuple
                if self.imv3.view.sceneBoundingRect().contains(pos):
                    mousePoint = self.imv3.view.mapSceneToView(pos)
                    self.ui.lineEdit_h.setText("{:.2f}".format(mousePoint.x()))
                    self.ui.lineEdit_v.setText("{:.2f}".format(mousePoint.y()))
                    if abs(mousePoint.x()) < self.ds.FoVX[self.nCurrentImage]/2 and abs(mousePoint.y()) < self.ds.FoVY[self.nCurrentImage]/2:
                        Xindex = int((mousePoint.x()+self.ds.FoVX[self.nCurrentImage]/2)/self.ds.PixelSpacingX[self.nCurrentImage]) #if self.ds.PixelSpacingX[self.nCurrentImage] > 0. else Xindex = int(mousePoint.x())
                        Yindex = int((mousePoint.y()+self.ds.FoVY[self.nCurrentImage]/2)/self.ds.PixelSpacingY[self.nCurrentImage]) #if self.ds.PixelSpacingY[self.nCurrentImage] > 0. else Yindex = int(mousePoint.y())
                        value=  self.ds.PA[self.nCurrentImage][Xindex,Yindex]      
                        self.ui.lineEdit_value.setText("{:.1f}".format(value))
                        rc= self.ReltoGlobal(mousePoint.x(), mousePoint.y(), self.nCurrentImage, self.ds)
                        self.ui.lineEdit_x.setText("{:.2f}".format(rc[0]))
                        self.ui.lineEdit_y.setText("{:.2f}".format(rc[1]))
                        self.ui.lineEdit_z.setText("{:.2f}".format(rc[2]))
                    self.imv3.vLine.setPos(mousePoint.x())
                    self.imv3.hLine.setPos(mousePoint.y())
             
    def mouseMoved4(self,evt): #mouse move event to move crosshairs and display location and values
        if self.dataSetIsNew == False:
            if (self.windows[3] == 1):
                self.ds = self.dataSet[3]
                pos = evt[0]  ## using signal proxy turns original arguments into a tuple
                if self.imv4.view.sceneBoundingRect().contains(pos):
                    mousePoint = self.imv4.view.mapSceneToView(pos)
                    self.ui.lineEdit_h.setText("{:.2f}".format(mousePoint.x()))
                    self.ui.lineEdit_v.setText("{:.2f}".format(mousePoint.y()))
                    if abs(mousePoint.x()) < self.ds.FoVX[self.nCurrentImage]/2 and abs(mousePoint.y()) < self.ds.FoVY[self.nCurrentImage]/2:
                        Xindex = int((mousePoint.x()+self.ds.FoVX[self.nCurrentImage]/2)/self.ds.PixelSpacingX[self.nCurrentImage]) #if self.ds.PixelSpacingX[self.nCurrentImage] > 0. else Xindex = int(mousePoint.x())
                        Yindex = int((mousePoint.y()+self.ds.FoVY[self.nCurrentImage]/2)/self.ds.PixelSpacingY[self.nCurrentImage]) #if self.ds.PixelSpacingY[self.nCurrentImage] > 0. else Yindex = int(mousePoint.y())
                        value=  self.ds.PA[self.nCurrentImage][Xindex,Yindex]
                        self.ui.lineEdit_value.setText("{:.1f}".format(value))
                        rc= self.ReltoGlobal(mousePoint.x(), mousePoint.y(), self.nCurrentImage, self.ds)
                        self.ui.lineEdit_x.setText("{:.2f}".format(rc[0]))
                        self.ui.lineEdit_y.setText("{:.2f}".format(rc[1]))
                        self.ui.lineEdit_z.setText("{:.2f}".format(rc[2]))
                    self.imv4.vLine.setPos(mousePoint.x())
                    self.imv4.hLine.setPos(mousePoint.y())
 
    def ReltoGlobal (self, h,v,n, dstype):   #given relative coordinate x,y of image n returns np vector of global coordinates 
        rc= ((h+dstype.FoVX[n]/2) * dstype.RowDirection[n]+(v+dstype.FoVX[n]/2)*dstype.ColumnDirection[n])+dstype.ImagePosition[n]
        return rc

# UNUSED FUNCTION
    def GlobaltoRel(self,r,n, dstype):    #Given r vector in global coordinates returns h,v in image plane of image n
        h=np.dot(r-dstype.ImageCenter[n],dstype.RowDirection[n])  
        v=np.dot(r-dstype.ImageCenter[n],dstype.ColumnDirection[n])
        return [h,v]

if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    main = Recon()
    main.show()
    sys.exit(app.exec_())

    