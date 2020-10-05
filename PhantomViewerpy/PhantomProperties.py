# -*- coding: utf-8 -*-
"""
Created on Sat Nov 09 11:00:01 2013
Reads in and display Phantom properties
Uses ROIPropertiesGui.py created from PhantomPropertiesGui.ui by QT4
execute   "designer\pyuic4 designer\PhantomPropertiesGui.ui -o PhantomViewer\PhantomPropertiesGui.py" from system shell to regenerate ROIPropertiesGui.py from ROIPropertiesGui.ui
@author: stephen russek
"""

import sys
try:
  from PyQt4 import QtGui, QtCore
except:
  from PyQt5 import QtGui, QtCore
from PhantomPropertiesGui import Ui_PhantomPropertiesGui
import numpy as np
import VPhantom

class PhantomProperties(QtGui.QMainWindow):
  def __init__(self ,phantom, parent = None):
    super(PhantomProperties, self).__init__()
    self.ui = Ui_PhantomPropertiesGui()
    self.ui.setupUi(self)
    self.setWindowTitle('Phantom Properties')
    self.Phantom= phantom
    
  #signals and slots
    self.ui.actionOpen_Phantom_File.triggered.connect(self.openPhantomFile)
    self.ui.actionSave_Phantom_File.triggered.connect(self.savePhantomFile)

    self.showPhantomInfo(self.Phantom)
 
  def openPhantomFile (self, direct =''):
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
  
 
  def savePhantomFile (self,phantom):
      fileName = QtGui.QFileDialog.getSaveFileName(parent=None, caption="Report File Name", directory = '', selectedFilter = ".dat")
      if not fileName:  #if cancel is pressed return
        return None
      f= open(fileName, 'w')
      s = self.Phantom.printROIinfo()
      f.write(s)
      f.close()
      print (s)
            
  def showPhantomInfo (self, phantom):
    self.Phantom=phantom
    self.ui.txtPhantomName.setText(self.Phantom.phantomName)
    self.ui.txtComment.setText(self.Phantom.Comment)
    self.ui.txtTemperature.setText(str(self.Phantom.Temperature))
    self.ui.txtnROIsets.setText(str(len(self.Phantom.ROIsets)))
    self.ui.txtField.setText(str(self.Phantom.B0))
    self.ui.txtPhantomProperties.setText(self.Phantom.printROIinfo())
    self.ui.lblPhantomImage.setPixmap(QtGui.QPixmap(self.Phantom.phantomImage))
    self.ui.lblPhantomImage.setScaledContents(True)
    self.ui.lblPhantomImage.show()

 
                  
if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    test = PhantomProperties()
    test.show()
    sys.exit(app.exec_())
