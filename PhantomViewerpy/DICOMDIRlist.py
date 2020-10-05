'''
Created on Dec 29, 2014
Uses DICOMDIRGui.py created from DICOMDIRGui.ui by QT4
  execute   "pyuic4 DICOMDIRGui.ui -o DICOMDIRGui.py" 
@author: stephen russek
'''
try:
  from PyQt4 import QtGui, QtCore
except:
  from PyQt5 import QtGui, QtCore
from DICOMDIRGui import Ui_DICOMDIRGui

class DICOMDIRlist(QtGui.QDialog):
    def __init__(self ,parent = None):
        super(DICOMDIRlist, self).__init__()
        self.ui = Ui_DICOMDIRGui()
        self.ui.setupUi(self)
        self.setWindowTitle('DICOMDIR')
        self.list=self.ui.listDICOMDIR
        
        
    def selectedSeries(self):
      '''returns selected series numbers'''
      s=self.list.selectedItems()
      sn=""
      for item in s:
        ss=str(item.text())
        sn =sn + ss[ss.find('=')+1: ss.find(",")+1]
      return sn


 
