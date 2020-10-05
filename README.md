# PhantomViewer
Analysis package written in Python to analyze MRI phantoms developed by NIST and collaborators at NCI, RSNA, NIBIB

Created on Fri Oct 11 16:30:54 2013

PhantomViewer: Main class to display and analyze phantom data, mostly MRI data

Modules required (most can be installed using 'pip install package' or  'pip install package --upgrade':
verify module version using "pip show module"

   Anaconda 3 with Python 3.7    (64 bit, Python, PyQT5, Numpy, Scipy, Pillow )
   
   pydicom Version: 1.4.2:       (for DICOM file import/export)
   
   pyqtgraph Version: 0.10.0     (for plotting and image windows)
   
   lmfit Version: 1.0.0          (for nonlinear least squares fitting routines)
   
   unwrap Version: 0.1.1
   
   scikit-image  Version: 0.11.3: uses unwrap_phase from skimage.resoration
   
   PyOpenGL                      (Used for 3D plotting)

will work with Anaconda 2 with Python 2.7.11, PYQT4 (not recommended)

Uses PhantomViewerGui created from PhantomViewer.ui by Qt Designer

  convert ui file to python by executing:
  
      "designer\pyuic4 designer\PhantomViewerGui.ui -o PhantomViewerpy\PhantomViewerGui4.py"
      
  or  "designer\pyuic5 designer\PhantomViewerGui.ui -o PhantomViewerpy\PhantomViewerGui5.py"
  
  from system shell to regenerate PhantomViewerGui.py from PhantomViewerGui.ui
  
@author: Stephen Russek

Units: times in ms, distances in mm, ADC in mm2/s, Temperature in C,

VPhantom: class to describe phantoms, several virtual phantoms, including the NIST/ISMRM system phantom are available.
PV structure is shown below:

![PV Structure](https://github.com/StephenRussek/PhantomViewerNew/blob/master/icons/PVstructure.jpg)
