# -*- coding: utf-8 -*-
"""
Created on Thu May  4 16:48:43 2017

@author: S.Y. Agustsson
"""

from PyQt5 import QtGui, QtCore, QtWidgets  # (the example applies equally well to PySide)
from PyQt5.QtWidgets import QMessageBox
import sys
import pyqtgraph as pg
import numpy as np

class testGUIwindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()

        self.initUI()

    def initUI(self):
        pass


class testGUI(QtWidgets.QWidget):

    def __init__(self):
        super().__init__()

        self.initUI()

    def initUI(self):
        """ Generate GUI layout """

        # Define top level widget
        self.w = QtGui.QWidget()
        #self.showFullScreen()

        layout = QtGui.QGridLayout() # create a grid for subWidgets
        self.col = QtGui.QColor(0,0,0)
        self.setLayout(layout)


        btn = QtGui.QPushButton('do something',self)
        btn.clicked.connect(self.asdf())
        btn.setToolTip('Press to <b>do nothng</b>')
        btn.resize(btn.sizeHint())

#        qbtn = QtGui.QPushButton('Quit',self)
#        qbtn.clicked.connect(QtCore.QCoreApplication.instance().quit)
#        btn.setToolTip('Press to <b>quit</b>')
#        btn.resize(btn.sizeHint())
        # Set Layout
        layout.addWidget(btn, 0, 0)
#        layout.addWidget(qbtn, 1, 0)

    def asdf(self):
        reply = QMessageBox.question(self, 'Message', 'can you read me?', QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QtGui.QMessageBox.Yes: print('yes')



    def closeEvent(self, event):

        reply = QMessageBox.question(self, 'Message',
            "Are you sure to quit?", QMessageBox.Yes |
            QMessageBox.No, QMessageBox.No)

        if reply == QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()


if __name__ == '__main__':

    ## Always start by initializing Qt (only once per application)
    #app = QtGui.QApplication([])

#    x = np.arange(0,1000,1)
#    noise = np.random.normal(0,1,1000)/1
#    y = np.sin(x/10)+noise
     # assign previous initialization of Qt or initialize it otherwise
    app = QtCore.QCoreApplication.instance()
    if app is None:
        app = QtGui.QApplication(sys.argv)
    # Create handle prg for the Graphic Interface
    prg = testGUI()

    prg.show()
    # Start Qt event loop
    #sys.exit(app.exec_())
    app.exec_()