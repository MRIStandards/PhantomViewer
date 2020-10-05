
"""Make code compatible with PyQt4 and PyQt5
which have components in different module
First see if you can find PyQT4, use this, else use PyQt5"""


try:    #to make compatible with Qt4 and Qt5 
  pyqtVersion=4    #flag to indicate which gui files to use
  from PyQt4.Qt import PYQT_VERSION_STR
  from PyQt4.QtCore import Qt, QPoint
  from PyQt4.QtGui import QApplication, QMainWindow, QProgressDialog, QInputDialog, QFileDialog, QColorDialog, QAction, QTextEdit, QStatusBar, QMenuBar
  from PyQt4.QtGui import  QFont, QColor, QPixmap
except:
  pyqtVersion=5    #flag to indicate which gui files to use
  from PyQt5.Qt import PYQT_VERSION_STR
  from PyQt5.QtCore import  Qt, QPoint
  from PyQt5.QtGui import  QFont, QColor, QPainter, QPixmap, QTextOption, QScreen
  from PyQt5.QtWidgets import QApplication, QMainWindow,  QWidget, QProgressDialog, QInputDialog, QColorDialog, QLineEdit, QFileDialog, QAction, QTextEdit, QToolTip, QStatusBar, QMenuBar
  QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True) #enable highdpi scaling
  QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True) #use highdpi icons