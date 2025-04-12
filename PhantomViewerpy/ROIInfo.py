r'''
Created on Mar 3, 2015
  execute   "designer\pyuic5 designer\ROIInfoGui.ui -o PhantomViewerpy\ROIInfoGui5.py" from system shell to regenerate ROIViewGui.py from ROIViewGui.ui
            
@author: stephen russek

Window to plot ROI voxels, statistics, and properties
'''
from pyqt import *         #imports required PyQt modules
if pyqtVersion==5:
  from ROIInfoGui5 import Ui_ROIInfoWindow    #main window Gui
import numpy as np
import numpy.ma as ma   #masked arrays eg ROI arrays 
import pyqtgraph as pg        

class ROIInfoWindow(QMainWindow):
    def __init__(self ,  parent = None):
        super(ROIInfoWindow, self).__init__()
        self.ui = Ui_ROIInfoWindow()
        self.ui.setupUi(self)
        self.setWindowTitle('ROI Info')
        self.ui.pbUpdate.clicked.connect(self.saveChanges)
        self.imv=self.ui.imvROI
        self.imv.ui.roiBtn.setText("Line scan/ ROI")
#        self.imv.ui.histogram.plot.setLogMode(None,True)    #set the histogram y axis to a log scale    
        self.imv.vLine = pg.InfiniteLine(angle=90, movable=False)   #cross hair
        self.imv.hLine = pg.InfiniteLine(angle=0, movable=False)
        self.imv.addItem(self.imv.vLine, ignoreBounds=True)
        self.imv.addItem(self.imv.hLine, ignoreBounds=True)
        self.proxy = pg.SignalProxy(self.ui.imvROI.view.scene().sigMouseMoved, rateLimit=60, slot=self.mouseMoved)
        
    def update(self,roiset, roi, roiview):
      self.roi=roi
      self.roiview=roiview
      self.ui.lblROISet.setText(roiset.ROIName)
      self.ui.lblROIType.setText(roi.Type)
      self.ui.lblROIIndex.setText(str(roi.Index))
      self.ui.txtXCenter.setText("{:.2f}".format(roi.Xcenter))
      self.ui.txtYCenter.setText("{:.2f}".format(roi.Ycenter))
      self.ui.txtZCenter.setText("{:.2f}".format(roi.Zcenter))
      self.ui.txtd1.setText("{:.2f}".format(roi.d1))  
      self.ui.txtT1.setText("{:.2f}".format(roi.T1))
      self.ui.txtT2.setText("{:.2f}".format(roi.T2))
      self.ui.txtADC.setText("{:.2f}".format(roi.ADC*1000))
      self.ui.txtConcentration.setText("{:.2f}".format(roi.Concentration))
      self.ui.txtProtonDensity.setText("{:.2f}".format(roi.PD))
      self.ui.lblAve.setText("{:.3f}".format(roi.SignalAve))
      self.ui.lblSd.setText("{:.3f}".format(roi.SignalRMS))
      self.roiArray=roi.array.filled(0.0)
      self.ui.imvROI.setImage(self.roiArray)
      self.ui.lenVoxels.setText(str(roi.array.count()))
      
    def saveChanges(self):
      roi=self.roi
      try:
        roi.Xcenter=float(self.ui.txtXCenter.toPlainText())
        roi.Ycenter=float(self.ui.txtYCenter.toPlainText())
        roi.Zcenter=float(self.ui.txtZCenter.toPlainText())
        roi.d1=float(self.ui.txtd1.toPlainText())  
        roi.T1=float(self.ui.txtT1.toPlainText())
        roi.T2=float(self.ui.txtT2.toPlainText())
        roi.ADC=float(self.ui.txtADC.toPlainText())/1000
        roi.Concentration=float(self.ui.txtConcentration.toPlainText())
        roi.PD=float(self.ui.txtProtonDensity.toPlainText())
      except:
        pass

      self.roiview.redrawROIs()
   

    def mouseMoved(self,evt): 
        '''mouse move event to move crosshairs and display location and values'''
        data = self.imv.image  
        nRows, nCols = data.shape
        pos = evt[0]  ## using signal proxy turns original arguments into a tuple
        if self.ui.imvROI.view.sceneBoundingRect().contains(pos):
            #mousePoint = self.ui.imvROI.view.vb.mapSceneToView(pos)
            mousePoint = self.ui.imvROI.view.mapSceneToView(pos)

#             if abs(mousePoint.x()) < self.ds.FoVX[self.nCurrentImage]/2 and abs(mousePoint.y()) < self.ds.FoVY[self.nCurrentImage]/2:
#                 Xindex = int((mousePoint.x()+self.ds.FoVX[self.nCurrentImage]/2)/self.ds.PixelSpacingX[self.nCurrentImage]) #if self.ds.PixelSpacingX[self.nCurrentImage] > 0. else Xindex = int(mousePoint.x())
#                 Yindex = int((mousePoint.y()+self.ds.FoVY[self.nCurrentImage]/2)/self.ds.PixelSpacingY[self.nCurrentImage]) #if self.ds.PixelSpacingY[self.nCurrentImage] > 0. else Yindex = int(mousePoint.y())
#                 value=  self.ds.PA[self.nCurrentImage][Xindex,Yindex]      
#                 self.ui.lblValue.setText("{:.2f}".format(value))
#                 rc= self.reltoGlobal(mousePoint.x(), mousePoint.y(), self.nCurrentImage)
#                 self.ui.lblX.setText("{:.2f}".format(rc[0]))
#                 self.ui.lblY.setText("{:.2f}".format(rc[1]))
#                 self.ui.lblZ.setText("{:.2f}".format(rc[2]))
            col=int(mousePoint.x())
            row=int(mousePoint.y())
            self.ui.leH.setText("{:.2f}".format(col))
            self.ui.leV.setText("{:.2f}".format(row))
            self.imv.vLine.setPos(mousePoint.x())
            self.imv.hLine.setPos(mousePoint.y())
            if (0 <= row < nRows) and (0 <= col < nCols):
              self.ui.leValue.setText("{:.2f}".format(data[col, row]))
              #print("pos = ({:d}, {:d}), value = {!r}".format(row, col, value))
            else:
              pass
              #print("no data at cursor") 
 
