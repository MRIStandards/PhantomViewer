# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'designer\InfoGui.ui'
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

class Ui_InfoWindow(object):
    def setupUi(self, InfoWindow):
        InfoWindow.setObjectName(_fromUtf8("InfoWindow"))
        InfoWindow.resize(896, 747)
        self.centralwidget = QtGui.QWidget(InfoWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.verticalLayout_2 = QtGui.QVBoxLayout(self.centralwidget)
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.label = QtGui.QLabel(self.centralwidget)
        self.label.setAutoFillBackground(True)
        self.label.setFrameShape(QtGui.QFrame.WinPanel)
        self.label.setFrameShadow(QtGui.QFrame.Sunken)
        self.label.setText(_fromUtf8(""))
        self.label.setObjectName(_fromUtf8("label"))
        self.verticalLayout.addWidget(self.label)
        self.lblInfo = QtGui.QLabel(self.centralwidget)
        self.lblInfo.setAutoFillBackground(True)
        self.lblInfo.setText(_fromUtf8(""))
        self.lblInfo.setScaledContents(True)
        self.lblInfo.setObjectName(_fromUtf8("lblInfo"))
        self.verticalLayout.addWidget(self.lblInfo)
        self.horizontalLayout.addLayout(self.verticalLayout)
        self.verticalScrollBar = QtGui.QScrollBar(self.centralwidget)
        self.verticalScrollBar.setOrientation(QtCore.Qt.Vertical)
        self.verticalScrollBar.setObjectName(_fromUtf8("verticalScrollBar"))
        self.horizontalLayout.addWidget(self.verticalScrollBar)
        self.verticalLayout_2.addLayout(self.horizontalLayout)
        InfoWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(InfoWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 896, 26))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        InfoWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(InfoWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        InfoWindow.setStatusBar(self.statusbar)

        self.retranslateUi(InfoWindow)
        QtCore.QMetaObject.connectSlotsByName(InfoWindow)

    def retranslateUi(self, InfoWindow):
        InfoWindow.setWindowTitle(_translate("InfoWindow", "Phantom Viewer Information", None))

