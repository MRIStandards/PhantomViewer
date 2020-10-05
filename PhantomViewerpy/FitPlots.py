'''
Created on Dec 30, 2014
Uses FitPlotsGui.py created from FitPlots.ui by QT4
  execute   "designer\pyuic4 designer\FitPlotsGui.ui -o PhantomViewerpy\FitPlotsGui4.py"
  "designer\pyuic5 designer\FitPlotsGui.ui -o PhantomViewerpy\FitPlotsGui5.py" 
@author: stephen russek
'''
from pyqt import *         #imports required PyQt modules, tries PyQT4 then PyQt5
import pyqtgraph as pg
#import pyqtgraph.opengl as gl
import pyqtgraph.functions as fn
if pyqtVersion==4:
  from FitPlotsGui4 import Ui_FitPlotsGui    #main window Gui
if pyqtVersion==5:
  from FitPlotsGui5 import Ui_FitPlotsGui    #main window Gui
import numpy as np

class FitPlots(QMainWindow):
    def __init__(self , x, y, fx, fy, resy, xlabel='Time(ms)', ylabel='Signal',header='', parent = None):
        super(FitPlots, self).__init__()
        self.ui = Ui_FitPlotsGui()
        self.ui.setupUi(self)
        self.setWindowTitle('Fits')
        self.fitPlot=self.ui.gvFitPlot
        self.ui.vsROI.setMaximum(y.shape[0])  #set vertical slider maximum to number of curves
        self.ui.vsROI.valueChanged.connect(self.plotROIdata)
        self.ui.actionSave.triggered.connect(self.saveData)
        self.ui.actionData.triggered.connect(self.dataTrue)
        self.ui.actionResiduals.triggered.connect(self.residualsTrue) 
        self.x=x
        self.y=y
        self.fx=fx
        self.fy=fy
        self.resy=resy
        self.header=header
        self.data=True
        self.residuals=False
        self.fitPlot.setLabel('bottom', xlabel)
        self.fitPlot.setLabel('left', ylabel)
        self.fitPlot.plot(self.x, self.y[0,:], pen=None,  symbolPen=None, symbolSize=20, symbolBrush=(255, 0, 0)) 
        self.fitPlot.plot(self.fx, self.fy[0,:],pen={'color': (255,0,0), 'width': 3}, symbol =None)
        self.inf1 = pg.InfiniteLine(movable=True, angle=90, label='x={value:0.3f}', 
                       labelOpts={'position':0.1, 'color': (200,200,100), 'fill': (200,200,200,50), 'movable': True})
        self.inf2 = pg.InfiniteLine(movable=True, angle=0, pen=(0, 0, 200),  hoverPen=(0,200,0), label='y={value:0.3f}', 
                       labelOpts={'color': (200,0,0), 'movable': True, 'fill': (0, 0, 200, 100)})
    def plotROIdata(self):
      nROI=self.ui.vsROI.value()-1
      self.ui.gvFitPlot.clear()
      if self.data:
          self.fitPlot.plot(self.x, self.y[nROI,:], pen=None,  symbolPen=None, symbolSize=20, symbolBrush=(255, 0, 0))
          self.fitPlot.plot(self.fx, self.fy[nROI,:],pen={'color': (255,0,0), 'width': 3}, symbol =None)
          self.fitPlot.addItem(self.inf1)
          self.fitPlot.addItem(self.inf2)
      if self.residuals:
          self.fitPlot.plot(self.x, self.resy[nROI,:]/np.amax(self.fy[nROI,:]), pen=None,  symbolPen=None, symbolSize=20, symbolBrush=(255, 0, 0))
          self.fitPlot.addItem(self.inf1)
          self.fitPlot.addItem(self.inf2)
          
          
          
    def dataTrue(self):
        self.data=True
        self.residuals=False
        self.plotROIdata()
        
    def residualsTrue(self):
        self.data=False
        self.residuals=True
        self.plotROIdata()
   
    def saveData(self):
      '''output data to fileName.csv and fit function to filenameFit.csv'''
      fileName = QFileDialog.getSaveFileName(parent=None, caption="Report File Name", directory = '', selectedFilter = ".csv")
      if not fileName:  #if cancel is pressed return
        return None
      file=str(fileName)
      idx=file.index('.')
      fitfile=file[:idx]+'Fit'+file[idx:]
      header=self.header+'time(ms)'
      for i in range(self.y.shape[0]):
        header+=',ROI' +str(i+1)      
      a1=self.x.reshape(self.x.shape[0],-1)
      a2=np.transpose(self.y) 
      np.savetxt(file,np.column_stack((a1,a2)) ,fmt='%.6e', header=header, delimiter=',')
      header=self.header+'time(ms)'
      for i in range(self.fy.shape[0]):
        header+=',ROI' +str(i+1) + 'fit'       
      a1=self.fx.reshape(self.fx.shape[0],-1)
      a2=np.transpose(self.fy) 
      np.savetxt(fitfile,np.column_stack((a1,a2)) ,fmt='%.6e', header=header, delimiter=',')
      for i in range(self.fy.shape[0]):
        header+=',ROI' +str(i+1) + 'Residuals'       
      a1=self.x.reshape(self.x.shape[0],-1)
      a2=np.transpose(self.resy) 
      np.savetxt(fitfile,np.column_stack((a1,a2)) ,fmt='%.6e', header=header, delimiter=',')