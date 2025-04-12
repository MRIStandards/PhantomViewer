'''
Created on Mar 3, 2015
  execute   "designer\pyuic4 designer\InfoGui.ui -o PhantomViewer\InfoGui.py" from system shell to regenerate ROIViewGui.py from ROIViewGui.ui
@author: stephen russek
'''

from PyQt5 import QtGui, QtCore
from InfoGui import Ui_InfoWindow

class InfoWindow(QtGui.QMainWindow):
    def __init__(self ,  parent = None):
        super(InfoWindow, self).__init__()
        self.ui = Ui_InfoWindow()
        self.ui.setupUi(self)
        self.setWindowTitle('Phantom Info')
        