# -*- coding: utf-8 -*-
"""
Created on Thu May  4 16:48:43 2017

@author: S.Y. Agustsson
"""

from PyQt5 import QtGui, QtCore, QtWidgets  # (the example applies equally well to PySide)
from PyQt5.QtWidgets import QInputDialog, QMessageBox, QLineEdit, QFileDialog, QMainWindow, QTableWidget, QTableWidgetItem
from PyQt5.QtCore import pyqtSlot
import sys
import pyqtgraph as pg
import numpy as np
from functionlibrary import redred as rr


class testGUI(QtWidgets.QWidget):

    def __init__(self):
        super().__init__()
        self.title = 'the GUI test'
        self.left = 300
        self.top = 150
        self.width = 640
        self.height = 480
        self.initUI()


    def initUI(self):
        """ Generate GUI layout """
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
#        self.statusBar().showMessage('Idle')
#        # Define top level widget
#        self.w = QtGui.QWidget()
#        #self.showFullScreen()

        self.widget = QtGui.QWidget()
        layout = QtGui.QGridLayout() # create a grid for subWidgets
        self.setLayout(layout)


        self.loadButton = QtGui.QPushButton('Open File', self)
        self.loadButton.resize(self.loadButton.sizeHint())
        self.loadButton.clicked.connect(self.loadScanRaw)

        self.clrButton = QtGui.QPushButton('Clear Plot', self)
        self.clrButton.resize(self.clrButton.sizeHint())
        self.clrButton.clicked.connect(self.clearPlot)



#        self.textbox = QLineEdit(self)
#        self.textbox.setPlaceholderText('enter parameters')
#        self.textbox.resize(self.textbox.sizeHint())


        self.metadata = {'list1': [1,2,3,4,5,6, {'nested1': 'aaaaa', 'nested2': 'bbbbb'}, "seven"],
                         'dict1': {'x': 1, 'y': 2, 'z': 'three'},
                         'array1 (20x20)': np.ones((10,10))
                        }

        self.tree = pg.DataTreeWidget()

        self.createPlotWidget()



        # Set Layout
        layout.addWidget(self.loadButton, 2, 4, 1, 1)
        layout.addWidget(self.clrButton, 2, 3, 1, 1)
        #layout.addWidget(self.textbox,   2, 2, 1, 1)
        layout.addWidget(self.plotWidget, 0, 3, 2,4)
        layout.addWidget(self.tree,       0,0,2,2)



        self.show()

    def createPlotWidget(self):
        self.plotWidget = pg.PlotWidget()
        self.data = []

        pg.setConfigOption('background', 'w')
        pg.setConfigOption('foreground', 'k')


    def loadScanCSV(self):

        filename = self.openFileNameDialog()
        self.scanData = rr.rrScan()
        self.scanData.importCSV(filename)
        self.scanData.initParameters()
        self.tree.addChildren(self.scanData.parameters)
#        self.plotData()
        x = self.scanData.time
        y = self.scanData.trace
        self.plot = self.plotWidget.plot(x,y, color='r')
        print(self.scanData.trace)

    def loadScanRaw(self):

        filename = self.openFileNameDialog()
        self.scanData = rr.rrScan()
        self.scanData.importRawFile(filename)
        self.scanData.initParameters()
        print(self.scanData.parameters)
        #self.tree.addChild(self.scanData.parameters)
#        self.plotData()
        x = self.scanData.time
        y = self.scanData.trace
        self.plot = self.plotWidget.plot(x,y, color='r')



    def clearPlot(self):
        """clears all graphs from plot, after asking confermation"""
        reply = QMessageBox.question(self, 'Message',
            "Are you sure you want to clear the graph completely?", QMessageBox.Yes |
            QMessageBox.No, QMessageBox.No)

        if reply == QMessageBox.Yes:
            self.plotWidget.clear()


    def plotDataTest(self):
        """ plot curves in the plot widget"""
        x = np.arange(0,1000,1)
        noise = np.random.normal(0,1,1000)/1
        y = np.sin(x/10)+noise

        self.plot = self.plotWidget.plot(x,y, color='g')



    def plotData(self):
        """ plot curves in the plot widget"""
        x = self.scanData.time
        y = self.scanData.trace
        self.plot = self.plotWidget.plot(x,y, color='b')

    def openFileNameDialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;Python Files (*.py)", options=options)
        if fileName:
            return(fileName)

    def set_statusbar(self):
        buttonReply =  QMessageBox.question(self, 'PyQt5 message', "Do you like PyQt5?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if buttonReply == QMessageBox.Yes:
            string = 'yes'
        else:
            string = 'no'
        self.statusBar().showMessage(string)

    def closeEvent(self, event):
        """ to ask confermation when quitting program"""

        reply = QMessageBox.question(self, 'Message',
            "Are you sure to quit?", QMessageBox.Yes |
            QMessageBox.No, QMessageBox.No)

        if reply == QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()


#%%other functions

    @pyqtSlot()
    def on_click(self):
        textboxValue = self.textbox.text()
        QMessageBox.question(self, 'Message', "You typed: " + textboxValue, QMessageBox.Ok, QMessageBox.Ok)
        self.textbox.setText("")

    def getText(self):
        text, okPressed = QInputDialog.getText(self, "Get text","Your name:", QLineEdit.Normal, "")
        if okPressed and text != '':
            print(text)


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