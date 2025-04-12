# -*- coding: utf-8 -*-
"""
Created on Tue Jun 03 10:56:26 2014
Main Phantom Viewer class
@author: Hannah Erdevig
Uses PhantomViewerGui created from PhantomVIewerGui.ui by QT4
  execute   "designer\pyuic4 designer\PVSetupGui.ui -o PhantomViewerpy\PVSetupGui.py" 
"""

import sys
from PyQt5 import QtGui, QtCore
from PVSetupGui import Ui_PhantomViewerMainWindow
import PhantomViewer
try:
  import Recon
except:
    pass
from GeneralInstructionGui import Ui_Form_instruction
from PhantomDescriptionGui import Ui_Form_description
import ROIProperties
from AboutGui import Ui_Form

# add in for set pixmap on pop-up windows (so that python doesn't trash then right away)
try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

class MainDialog(QtGui.QMainWindow):
    def __init__(self ,parent = None):
        super(MainDialog, self).__init__()
        self.ui = Ui_PhantomViewerMainWindow()  # make self.ui object an instance of the class defined in PV GUI
        self.ui.setupUi(self)                   
        self.setWindowTitle('Phantom Viewer')
        #signals and slots
        self.ui.pushButton_geometry.clicked.connect(self.Geometry)      # call Geometry() if "Geometric Analysis" button is clicked
        self.ui.pushButton_T1.clicked.connect(self.T1)                  # call T1() if "T1 Analysis" button is clicked
        self.ui.pushButton_T2.clicked.connect(self.T2)                  # call T2() if "T2 Analysis" button is clicked
        self.ui.pushButton_proton.clicked.connect(self.Proton)          # call Proton() if "Proton Density Array/SNR" button is clicked
        self.ui.pushButton_slice.clicked.connect(self.SliceProfile)     # call SliceProfile() if "Slice Profile Inset" button is clicked
        self.ui.pushButton_res.clicked.connect(self.ResolutionInset)    # call ResolutionInset() if "Resolution Inset" button is clicked
        self.ui.pushButton_dicom.clicked.connect(self.Dicom)            # call Dicom() if "Dicom Converter" button is clicked
        self.ui.pushButton_recon.clicked.connect(self.Reconstruction)   # call Reconstruction() if "Reconstruction" button is clicked
        self.ui.action_General.triggered.connect(self.General)                  # call General() if "General Instruction" menu button is pressed
        self.ui.actionPhantom_Description.triggered.connect(self.Description)   # call Description() if "Phantom Description" menu button is pressed
        self.ui.actionROI_properties.triggered.connect(self.ROIProperties)      # call ROIProperties() if "ROI Properties" menu button is pressed
        self.ui.actionAbout_PV.triggered.connect(self.About)                    # call About() if "About Phantom Viewer" menu button is pressed
        
    def Geometry(self): # opens PhantomViewer module with Geometry labels
        self.GeometryView = PhantomViewer.PhantomViewer(self)
        self.GeometryView.setWindowTitle('Geometry Phantom Analysis')
        self.GeometryView.rdPlot.setTitle("Geometry raw data")
        self.GeometryView.fitPlot.setTitle("Geometry Results")
        self.GeometryView.show()
 
    def T1(self):       # opens PhantomViewer module with T1 presets
        self.T1View = PhantomViewer.PhantomViewer( "T1")
        self.T1View.show()
         
    def T2(self):
        self.T2View = PhantomViewer.PhantomViewer('T2' )
        self.T2View.show()
  
    def Proton(self):
        self.ProtonView = PhantomViewer.PhantomViewer("PD-SNR")
        self.ProtonView.show()
   
    def SliceProfile(self):
        self.SliceProf = PhantomViewer.PhantomViewer("SP")
        self.SliceProf.setWindowTitle('Proton Density/SNR Phantom Analysis')
        self.SliceProf.show()
        self.SliceProf.rdPlot.setTitle("Proton Density raw data")
        self.SliceProf.fitPlot.setTitle("Proton Density Results")
        self.SliceProf.DataType = "Proton Density"
        
    def ResolutionInset(self):
        self.ResInset = PhantomViewer.PhantomViewer(self, "Proton Density")
        self.ResInset.setWindowTitle('Proton Density/SNR Phantom Analysis')
        self.ResInset.show()
        self.ResInset.rdPlot.setTitle("Proton Density raw data")
        self.ResInset.fitPlot.setTitle("Proton Density Results")
        self.ResInset.DataType = "Proton Density"
        
    def Dicom(self):
        self.DicomConvert = PhantomViewer.PhantomViewer(self, "Proton Density")
        self.DicomConvert.setWindowTitle('Proton Density/SNR Phantom Analysis')
        self.DicomConvert.show()
        self.DicomConvert.rdPlot.setTitle("Proton Density raw data")
        self.DicomConvert.fitPlot.setTitle("Proton Density Results")
        self.DicomConvert.DataType = "Proton Density"
        
    def Reconstruction(self):
        self.Reconstruct = Recon.Recon(self)
        self.Reconstruct.show()
        
    def General(self):
        Form_instruction = QtGui.QWidget()
        self.instructions = Ui_Form_instruction()
        self.instructions.setupUi(Form_instruction)
        Form_instruction.show()
        Form_instruction.label.setPixmap(QtGui.QPixmap(_fromUtf8("Phantom Viewer accessories/instructions.JPG")))
        
    def Description(self):
        Form_description = QtGui.QWidget()
        self.description = Ui_Form_description()
        self.description.setupUi(Form_description)
        Form_description.show()
        Form_description.label.setPixmap(QtGui.QPixmap(_fromUtf8("Phantom Viewer accessories/description.JPG")))
         
    def ROIProperties(self):
        self.prop = ROIProperties.ROIProperties(self)
        self.prop.show()
        self.prop.setWindowTitle("ROI Properties")
        
    def About(self):
        Form_aboutPV = QtGui.QWidget()
        self.aboutPV = Ui_Form()
        self.aboutPV.setupUi(Form_aboutPV)
        Form_aboutPV.show()
        Form_aboutPV.label.setPixmap(QtGui.QPixmap(_fromUtf8("Phantom Viewer accessories/description.JPG")))
         
        
if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    main = MainDialog()
    main.show()
    sys.exit(app.exec_())

    