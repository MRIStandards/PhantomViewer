# -*- coding: utf-8 -*-
"""
Created on Sat Nov 09 11:00:01 2013
Reads in and display ROI properties
Uses ROIPropertiesGui.py created from ROIropertiesGui.ui by QT4
execute   "pyuic4 ROIPropertiesGui.ui -o ROIPropertiesGui.py" from system shell to regenerate ROIPropertiesGui.py from ROIPropertiesGui.ui
@author: stephen russek
"""

import sys

from PyQt5 import QtGui, QtCore
from PyQt5.QtWidgets import QMainWindow
from ROIPropertiesGui import Ui_ROIPropertiesGui
import numpy as np
import VPhantom

class ROIProperties(QMainWindow):
  def __init__(self ,phantom, parent = None):
    super(ROIProperties, self).__init__()
    self.ui = Ui_ROIPropertiesGui()
    self.ui.setupUi(self)
    self.setWindowTitle('ROI Properties')
    self.Phantom= phantom
    for i in range(14):
        self.ui.tblT1.setColumnWidth(i,70)
        self.ui.tblT2.setColumnWidth(i,70)
    self.tblT1={"T1Array-Concentration":0,"T1Array-T1":1,"T1Array-R1":2,"T1Array-Xcenter":3,"T1Array-Ycenter":4,"T1Array-Zcenter":5}
    self.T1parameter={0:"Concentration",1:"T1",2:"R1",3:"Xcenter",4:"Ycenter",5:"Zcenter"}
    self.tblT2={"T2Array-Concentration":0,"T2Array-T2":1,"T2Array-R1":2,"T2Array-Xcenter":3,"T2Array-Ycenter":4,"T2Array-Zcenter":5}
  #signals and slots
    self.ui.actionOpen_ROI_File.triggered.connect(self.openROIFile)
    self.ui.actionSave_ROI_File.triggered.connect(self.saveROIFile)
    self.ui.pbT1ReflectX.clicked.connect(self.T1ReflectX)

    self.SetROIs()
 
  def openROIFile (self, direct =''):
    self.Phantom.T1ROIs.ROIs = []   #reset ROI list
    if direct == False:
        direct = ''
    self.fileName = QtGui.QFileDialog.getOpenFileName(self,"Open ROI File", direct, "ROI File (*.dat)")
    if not self.fileName:  #if cancel is pressed return
      return None
    f = open(str(self.fileName), 'r')
    for line in f:
      parameter = line[:line.find('=')]
      values=line[line.find('=')+1:]
      if parameter == "PhantomName":
        self.ui.txtPhantomName.setText(str(values))
      if parameter == "B0":
        self.ui.txtField.setText(str(values))
      if parameter == "Temperature":
        self.ui.txtTemperature.setText(str(values))
      if parameter == "Comment":
        self.ui.txtComment.setText(str(values))
      if parameter == "NumberofROIsets":
        self.ui.txtnROIsets.setText(str(values))
      if parameter == "ROIName":
        ROIName = parameter
      if parameter == "T1Array-nROIs":
        self.Phantom.T1ROIs.nROIs=int(values)
        for i in range(int(values)):    #create a list of T1 ROIs
            self.Phantom.T1ROIs.ROIs.append(VPhantom.ROI())
            self.Phantom.T1ROIs.ROIs[-1].Type= "T1"
            self.Phantom.T1ROIs.ROIs[-1].Index = i+1
            self.Phantom.T1ROIs.ROIs[-1].Name = "T1-" +str(i+1)
      for te in self.tblT1:
        if parameter == te:     #cycle through T1 table entries
          data = np.fromstring(values, sep=',')
          for i in range(data.size):
            self.ui.tblT1.setItem(self.tblT1[te],i,QtGui.QTableWidgetItem(str(data[i])))
            self.Phantom.T1ROIs.SetROIsParameter( i, self.T1parameter[self.tblT1[te]], data[i])
#            print "Set T1 parameter=" + self.T1parameter[self.tblT1[te]] + " =" + str(data[i]) + "; ROI number=" + str(i)
#            if te == "T1CustomROIT1(ms)":
#                self.ui.tblT1.setItem(self.tblT1[te]+1,i,QtGui.QTableWidgetItem("{:.3f}".format(1000/data[i]))) #Set R1
      for te in self.tblT2:
        if parameter == te:     #cycle through T2 table entries
          data = np.fromstring(values, sep=',')
          for i in range(data.size):
            self.ui.tblT2.setItem(self.tblT2[te],i,QtGui.QTableWidgetItem(str(data[i])))
#            if te == "T2CustomROIT2(ms)":
#                self.ui.tblT2.setItem(self.tblT2[te]+1,i,QtGui.QTableWidgetItem("{:.3f}".format(1000/data[i])))  
  
 
  def saveROIFile (self,phantom):
      fileName = QtGui.QFileDialog.getSaveFileName(parent=None, caption="Report File Name", directory = '', selectedFilter = ".dat")
      if not fileName:  #if cancel is pressed return
        return None
      f= open(fileName, 'w')
      s = self.Phantom.printROIinfo()
      f.write(s)
      f.close()
      print (s)
            
  def SetROIs (self):
    self.ui.txtPhantomName.setText(self.Phantom.phantomName)
    self.ui.txtComment.setText(self.Phantom.Comment)
    self.ui.txtTemperature.setText(str(self.Phantom.Temperature))
    self.ui.txtField.setText(str(self.Phantom.B0))
    for roi in self.Phantom.T1ROIs.ROIs:
            self.ui.tblT1.setItem(0,roi.Index-1,QtGui.QTableWidgetItem("{:.2f}".format(roi.Concentration)))
            self.ui.tblT1.setItem(1,roi.Index-1,QtGui.QTableWidgetItem("{:.2f}".format(roi.T1)))
            self.ui.tblT1.setItem(3,roi.Index-1,QtGui.QTableWidgetItem("{:.2f}".format(roi.Xcenter)))
            self.ui.tblT1.setItem(4,roi.Index-1,QtGui.QTableWidgetItem("{:.2f}".format(roi.Ycenter)))
            self.ui.tblT1.setItem(5,roi.Index-1,QtGui.QTableWidgetItem("{:.2f}".format(roi.Zcenter)))
#            self.ui.tblT1.setItem(self.tblT1[te]+1,i,QtGui.QTableWidgetItem("{:.3f}".format(1000/data[i]))) #Set R1
    for roi in self.Phantom.T2ROIs.ROIs:
            self.ui.tblT2.setItem(0,roi.Index-1,QtGui.QTableWidgetItem("{:.2f}".format(roi.Concentration)))
            self.ui.tblT2.setItem(1,roi.Index-1,QtGui.QTableWidgetItem("{:.2f}".format(roi.T2)))
            self.ui.tblT2.setItem(3,roi.Index-1,QtGui.QTableWidgetItem("{:.2f}".format(roi.Xcenter)))
            self.ui.tblT2.setItem(4,roi.Index-1,QtGui.QTableWidgetItem("{:.2f}".format(roi.Ycenter)))
            self.ui.tblT2.setItem(5,roi.Index-1,QtGui.QTableWidgetItem("{:.2f}".format(roi.Zcenter)))
#            self.ui.tblT1.setItem(self.tblT1[te]+1,i,QtGui.QTableWidgetItem("{:.3f}".format(1000/data[i]))) #Set R1

  def T1ReflectX(self): 
    for roi in self.Phantom.T1ROIs.ROIs:
            roi.Xcenter=-roi.Xcenter
            self.ui.tblT1.setItem(3,roi.Index-1,QtGui.QTableWidgetItem("{:.2f}".format(roi.Xcenter)))
    
                  
if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    test = ROIProperties()
    test.show()
    sys.exit(app.exec_())
