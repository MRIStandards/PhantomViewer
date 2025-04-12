
"""Import basic pyqt modules"""



pyqtVersion=5    #flag to indicate which pyqt version is in use
from PyQt5.Qt import PYQT_VERSION_STR
from PyQt5.QtCore import  Qt, QPoint
from PyQt5.QtGui import  QFont, QColor, QPainter, QPixmap, QTextOption, QScreen
from PyQt5.QtWidgets import QApplication, QMainWindow,  QWidget, QProgressDialog, QInputDialog, QColorDialog, QLineEdit, QFileDialog, QAction, QTextEdit, QToolTip, QStatusBar, QMenuBar, QMenu, QStyleFactory

QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True) #enable highdpi scaling
QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True) #use highdpi icons