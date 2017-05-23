# -*- coding: utf-8 -*-
"""
Created on Thu May 11 11:03:06 2017

@author: S.Y. Agustsson
"""

from PyQt5 import QtGui as qg  # (the example applies equally well to PySide)
from PyQt5 import QtWidgets as qw
from PyQt5 import QtCore as qc
import sys
import pyqtgraph as pg
import numpy as np
from functionlibrary import transient as tr
import qdarkstyle
import os
from GUI import TransientsWidgets as tw

def main():
    app = qc.QCoreApplication.instance()
    if app is None:
        app = qg.QApplication(sys.argv)
    # Create handle prg for the Graphic Interface
    prg = tw.MainWindow()
    prg.show()
    app.exec_()



if __name__ == '__main__':
    main()