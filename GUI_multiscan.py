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
from lib import redred as rr
import qdarkstyle
import os
from GUI import rrWidgets

class MainWindow(qw.QMainWindow):  # TODO delete file
    """ Main application window

    """

    def __init__(self):
        super().__init__()


        self.title = 'Transient Analyser'
        self.left = 300
        self.top = 100
        self.width = 1400
        self.height = 900
        self.initUI()

        #set the cool dark theme
        self.setStyleSheet(qdarkstyle.load_stylesheet_pyqt5())
        pg.setConfigOption('background', 0.1)
        pg.setConfigOption('foreground', 0.7)




    def initUI(self):

        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        layout = qg.QGridLayout() # create a grid for subWidgets
        layout.setSpacing(10)
        self.setLayout(layout)


        self.test = rrWidgets.rrScanImport()
        layout.addWidget(self.test,0,0)
        self.setCentralWidget(self.test)



        self.show()


class SingleScanWidget(qw.QWidget):
    """master widget for multiple scan import/export and treatment
    """

    def __init__(self):
        """ """
        super().__init__()











if __name__ == '__main__':

    app = qc.QCoreApplication.instance()
    if app is None:
        app = qg.QApplication(sys.argv)
    # Create handle prg for the Graphic Interface
    prg = MainWindow()
    prg.show()
    app.exec_()