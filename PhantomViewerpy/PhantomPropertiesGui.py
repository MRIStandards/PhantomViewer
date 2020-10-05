# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'designer\PhantomPropertiesGui.ui'
#
# Created by: PyQt4 UI code generator 4.11.4
#
# WARNING! All changes made in this file will be lost!

try:
  from PyQt4 import QtGui, QtCore
except:
  from PyQt5 import QtGui, QtCore

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_PhantomPropertiesGui(object):
    def setupUi(self, PhantomPropertiesGui):
        PhantomPropertiesGui.setObjectName(_fromUtf8("PhantomPropertiesGui"))
        PhantomPropertiesGui.resize(1283, 838)
        self.centralwidget = QtGui.QWidget(PhantomPropertiesGui)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.label_3 = QtGui.QLabel(self.centralwidget)
        self.label_3.setGeometry(QtCore.QRect(25, 5, 87, 16))
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.label_4 = QtGui.QLabel(self.centralwidget)
        self.label_4.setGeometry(QtCore.QRect(545, 15, 49, 16))
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.txtField = QtGui.QTextEdit(self.centralwidget)
        self.txtField.setGeometry(QtCore.QRect(535, 40, 71, 21))
        self.txtField.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txtField.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txtField.setObjectName(_fromUtf8("txtField"))
        self.label_5 = QtGui.QLabel(self.centralwidget)
        self.label_5.setGeometry(QtCore.QRect(262, 5, 55, 16))
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.txtComment = QtGui.QTextEdit(self.centralwidget)
        self.txtComment.setGeometry(QtCore.QRect(185, 30, 336, 36))
        self.txtComment.setObjectName(_fromUtf8("txtComment"))
        self.label_6 = QtGui.QLabel(self.centralwidget)
        self.label_6.setGeometry(QtCore.QRect(630, 15, 95, 16))
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.txtTemperature = QtGui.QTextEdit(self.centralwidget)
        self.txtTemperature.setGeometry(QtCore.QRect(630, 40, 71, 21))
        self.txtTemperature.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txtTemperature.setObjectName(_fromUtf8("txtTemperature"))
        self.label_8 = QtGui.QLabel(self.centralwidget)
        self.label_8.setGeometry(QtCore.QRect(755, 15, 86, 16))
        self.label_8.setObjectName(_fromUtf8("label_8"))
        self.txtnROIsets = QtGui.QTextEdit(self.centralwidget)
        self.txtnROIsets.setGeometry(QtCore.QRect(750, 40, 71, 21))
        self.txtnROIsets.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txtnROIsets.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txtnROIsets.setObjectName(_fromUtf8("txtnROIsets"))
        self.txtPhantomName = QtGui.QTextEdit(self.centralwidget)
        self.txtPhantomName.setGeometry(QtCore.QRect(10, 30, 161, 36))
        self.txtPhantomName.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txtPhantomName.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.txtPhantomName.setObjectName(_fromUtf8("txtPhantomName"))
        self.txtPhantomProperties = QtGui.QTextEdit(self.centralwidget)
        self.txtPhantomProperties.setGeometry(QtCore.QRect(0, 80, 821, 701))
        self.txtPhantomProperties.setLineWrapMode(QtGui.QTextEdit.NoWrap)
        self.txtPhantomProperties.setObjectName(_fromUtf8("txtPhantomProperties"))
        self.lblPhantomImage = QtGui.QLabel(self.centralwidget)
        self.lblPhantomImage.setGeometry(QtCore.QRect(830, 80, 450, 300))
        self.lblPhantomImage.setAutoFillBackground(True)
        self.lblPhantomImage.setFrameShape(QtGui.QFrame.WinPanel)
        self.lblPhantomImage.setText(_fromUtf8(""))
        self.lblPhantomImage.setScaledContents(False)
        self.lblPhantomImage.setObjectName(_fromUtf8("lblPhantomImage"))
        PhantomPropertiesGui.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(PhantomPropertiesGui)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1283, 21))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        self.menuFile = QtGui.QMenu(self.menubar)
        self.menuFile.setObjectName(_fromUtf8("menuFile"))
        PhantomPropertiesGui.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(PhantomPropertiesGui)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        PhantomPropertiesGui.setStatusBar(self.statusbar)
        self.actionOpen_Phantom_File = QtGui.QAction(PhantomPropertiesGui)
        self.actionOpen_Phantom_File.setObjectName(_fromUtf8("actionOpen_Phantom_File"))
        self.actionSave_Phantom_File = QtGui.QAction(PhantomPropertiesGui)
        self.actionSave_Phantom_File.setObjectName(_fromUtf8("actionSave_Phantom_File"))
        self.menuFile.addAction(self.actionOpen_Phantom_File)
        self.menuFile.addAction(self.actionSave_Phantom_File)
        self.menubar.addAction(self.menuFile.menuAction())

        self.retranslateUi(PhantomPropertiesGui)
        QtCore.QMetaObject.connectSlotsByName(PhantomPropertiesGui)

    def retranslateUi(self, PhantomPropertiesGui):
        PhantomPropertiesGui.setWindowTitle(_translate("PhantomPropertiesGui", "Phantom Properties", None))
        self.label_3.setText(_translate("PhantomPropertiesGui", "Phantom Name", None))
        self.label_4.setText(_translate("PhantomPropertiesGui", "Field (T)", None))
        self.label_5.setText(_translate("PhantomPropertiesGui", "Comment", None))
        self.label_6.setText(_translate("PhantomPropertiesGui", "Temperature (C)", None))
        self.label_8.setText(_translate("PhantomPropertiesGui", "# of ROI sets", None))
        self.menuFile.setTitle(_translate("PhantomPropertiesGui", "File", None))
        self.actionOpen_Phantom_File.setText(_translate("PhantomPropertiesGui", "Open Phantom File", None))
        self.actionSave_Phantom_File.setText(_translate("PhantomPropertiesGui", "Save Phantom File", None))

