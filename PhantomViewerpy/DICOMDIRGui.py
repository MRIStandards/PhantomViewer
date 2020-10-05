# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'designer\DICOMDIRGui.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_DICOMDIRGui(object):
    def setupUi(self, DICOMDIRGui):
        DICOMDIRGui.setObjectName("DICOMDIRGui")
        DICOMDIRGui.setWindowModality(QtCore.Qt.WindowModal)
        DICOMDIRGui.resize(685, 583)
        self.buttonBox = QtWidgets.QDialogButtonBox(DICOMDIRGui)
        self.buttonBox.setGeometry(QtCore.QRect(105, 545, 341, 32))
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.layoutWidget = QtWidgets.QWidget(DICOMDIRGui)
        self.layoutWidget.setGeometry(QtCore.QRect(5, 10, 676, 531))
        self.layoutWidget.setObjectName("layoutWidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.layoutWidget)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.lblDICOMDIR = QtWidgets.QLabel(self.layoutWidget)
        self.lblDICOMDIR.setAutoFillBackground(True)
        self.lblDICOMDIR.setFrameShape(QtWidgets.QFrame.WinPanel)
        self.lblDICOMDIR.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.lblDICOMDIR.setText("")
        self.lblDICOMDIR.setWordWrap(True)
        self.lblDICOMDIR.setObjectName("lblDICOMDIR")
        self.verticalLayout.addWidget(self.lblDICOMDIR)
        self.listDICOMDIR = QtWidgets.QListWidget(self.layoutWidget)
        self.listDICOMDIR.setStyleSheet("selection-background-color: rgb(255, 255, 127);\n"
"selection-color: qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:0, stop:0 rgba(0, 0, 0, 255), stop:1 rgba(255, 255, 255, 255));")
        self.listDICOMDIR.setAlternatingRowColors(True)
        self.listDICOMDIR.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.listDICOMDIR.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectItems)
        self.listDICOMDIR.setResizeMode(QtWidgets.QListView.Adjust)
        self.listDICOMDIR.setObjectName("listDICOMDIR")
        self.verticalLayout.addWidget(self.listDICOMDIR)

        self.retranslateUi(DICOMDIRGui)
        self.buttonBox.accepted.connect(DICOMDIRGui.accept)
        self.buttonBox.rejected.connect(DICOMDIRGui.reject)
        QtCore.QMetaObject.connectSlotsByName(DICOMDIRGui)

    def retranslateUi(self, DICOMDIRGui):
        _translate = QtCore.QCoreApplication.translate
        DICOMDIRGui.setWindowTitle(_translate("DICOMDIRGui", "DICOMDIR"))

