# -*- coding: utf-8 -*-
"""
Created on Thu May  4 16:48:43 2017

@author: S.Y. Agustsson
"""

from PyQt5 import * # QtGui, QtCore, QtWidgets  # (the example applies equally well to PySide)
from PyQt5.QtWidgets import * #QInputDialog, QMessageBox, QLineEdit, QFileDialog, QMainWindow, QTableWidget, QTableWidgetItem
from PyQt5.QtCore import * #pyqtSlot
import sys
import pyqtgraph as pg
import numpy as np
from functionlibrary import redred as rr
import qdarkstyle
import os

class testGUI(QtWidgets.QWidget):

    def __init__(self):
        super().__init__()
        self.title = 'the GUI test'
        self.left = 300
        self.top = 100
        self.width = 1400
        self.height = 900
        self.initUI()
        #set the cool dark theme
        self.setStyleSheet(qdarkstyle.load_stylesheet_pyqt5())


    def initUI(self):
        """ Generate GUI layout """
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        self.setLayout_ImportScan()
        self.show()


    def setLayout_ImportScan(self):
        """ Generate the GUI layout """

        layout = QtGui.QGridLayout() # create a grid for subWidgets
        layout.setSpacing(10)
        self.setLayout(layout)

        font = QtGui.QFont()
        #font.setFamily(_fromUtf8("FreeMono"))
        font.setBold(True)
        font.setPixelSize(15)


        # Define items
        self.importFileBtn = QPushButton('Import File', self)
        self.importFileBtn.clicked.connect(self.importFile)
        layout.addWidget(self.importFileBtn, 18, 1, 1, 7)

        self.saveasCSVBtn = QPushButton('Save as CSV', self)
        self.saveasCSVBtn.clicked.connect(self.saveasCSV)
        layout.addWidget(self.saveasCSVBtn, 16, 24, 1, 5)

        self.saveFigure = QPushButton('Save Figure', self)
        #self.saveFigure.resize(self.clrButton.sizeHint())
        self.saveFigure.clicked.connect(self.nofunctionyet)
        layout.addWidget(self.saveFigure, 18, 24, 1, 5)

        self.nameTxtbox = QLabel('Series Name:', self)
        self.nameTxtbox.setFont(font)
        layout.addWidget(self.nameTxtbox, 1, 1, 1, 7)
        self.nameTxtbox = QLineEdit(self)
        self.nameTxtbox.setPlaceholderText('file name')
        self.nameTxtbox.editingFinished.connect(self.nofunctionyet)
        layout.addWidget(self.nameTxtbox, 2, 1, 1, 7)

        self.metadataTree_name = QLabel('Metadata:', self)
        self.metadataTree_name.setFont(font)
        layout.addWidget(self.metadataTree_name, 3, 1)
        self.metadataTree = pg.DataTreeWidget()
        #self.metadataTree.setHeaderItem()
        layout.addWidget(self.metadataTree, 4, 1, 13, 7)

        self.plotWidget_name = QLabel('Plot', self)
        self.plotWidget_name.setFont(font)
        layout.addWidget(self.plotWidget_name, 1, 9)
        self.setPlotWidget()
        layout.addWidget(self.plotWidget, 2, 9, 13, 20)

        # plot modification buttons
        self.modifyBox_name = QLabel('Modify', self)
        layout.addWidget(self.modifyBox_name, 15, 10)
        self.modifyBox_name.setFont(font)

        self.flipXcb = QPushButton('Flip x', self)
        self.flipXcb.clicked.connect(self.flipX)
        layout.addWidget(self.flipXcb, 16, 10)

        self.remDCcb = QPushButton('Remove DC', self)
        self.remDCcb.clicked.connect(self.removeDC)
        layout.addWidget(self.remDCcb, 17, 10)

        self.shiftXcb = QPushButton('Set Time Zero', self)
        self.shiftXcb.clicked.connect(self.setTimeZero)
        layout.addWidget(self.shiftXcb, 18, 10)

        self.timeZero = 0
        self.shiftTimeZero_input = QLineEdit(self)
        self.shiftTimeZero_input.setValidator(QtGui.QDoubleValidator())
        #self.shiftTimeZero_input.textChanged.connect(self.setShiftTimeZero)
        self.shiftTimeZero_input.returnPressed.connect(self.setTimeZero)
        layout.addWidget(self.shiftTimeZero_input, 18, 11, 1, 1)

        # Filter

#        self.filterBox = QGroupBox('filter', self)
#        layout.addWidget(self.filterBox, 15, 17, 3, 3)

        self.filterBox_name = QLabel('Filter', self)
        self.filterBox_name.setFont(font)
        layout.addWidget(self.filterBox_name, 15, 17)
        self.filterLowPassFreq_name = QLabel('Frequency [THz]', self)
        layout.addWidget(self.filterLowPassFreq_name, 16,17)
        self.filterLowPassFreq = QLineEdit(self)
        self.filterLowPassFreq.setPlaceholderText('freq')
        self.filterLowPassFreq.returnPressed.connect(self.applyFilter)
        layout.addWidget(self.filterLowPassFreq, 16, 18, 1, 3)
        self.applyFilterBtn = QPushButton('Apply', self)
        self.applyFilterBtn.clicked.connect(self.applyFilter)
        layout.addWidget(self.applyFilterBtn, 18, 17, 1, 5)


    def nofunctionyet(self):
        self.msg = QMessageBox()
        self.msg.setIcon(QMessageBox.Warning)
        self.msg.setText("No, this does nothing yet")
        self.msg.show()

    def nameChanged(self, name):
        print(name)

    def test(self):
        txt = self.shiftTimeZero_input.text()
        num = float(txt)
        print(num)

    def setPlotWidget(self):
        ''' Create the widget for plotting scans and define it's properties'''
        pg.setConfigOptions(antialias=True)
        self.plotWidget = pg.PlotWidget()
        self.data = []

#        pg.setConfigOption('background', 'w')
#        pg.setConfigOption('foreground', 'k')





    def fetchMetadata(self):
        '''retrieve metadata from rrScan() and write it in metadataTree'''
        self.metadataTree.setData(self.scanData.fetchMetadata(), hideRoot=True)
        print(self.metadataTree.headerItem())
#%% Scan Modifcations

    def flipX(self):
        self.scanData.flipTime()
        self.plotScanData()

    def flipY(self):
        self.scanData.flipTrace()
        self.plotScanData()

#    def setShiftTimeZero(self, text):
#        self.timeZero = float(text)
#        print(self.timeZero)


    def shiftTimeZero(self):
        '''shift time of scan by timeZero, value given in the QLineEdit shiftTimeZero'''
        txt = self.shiftTimeZero_input.text()
        num = float(txt)

        self.timeZero = num
        self.scanData.shiftTime(self.timeZero)

        self.plotScanData()

    def setTimeZero(self):
        '''set value given in shiftTimeZero_input() as new time zero'''
        txt = self.shiftTimeZero_input.text()
        num = float(txt)
        self.scanData.shiftTime(-self.timeZero)
        self.timeZero = num
        self.scanData.shiftTime(self.timeZero)
        self.plotScanData()



    def applyFilter(self):
        '''get filter frequency from textbox and apply the filter to a single scan'''
        freq = float(self.filterLowPassFreq.text())
        self.scanData.trace = self.scanData.rawtrace
        nyqfreq = self.scanData.nyqistFreq()

        if freq != 0 and freq < nyqfreq:

            cutfactor = freq / nyqfreq
            self.scanData.filterit(cutHigh=cutfactor)

        self.plotScanData()

    def removeDC(self):
        self.scanData.removeDC()
        self.plotScanData()




#%% Plots
    def plotScanData(self):
        """ clears the graph and plots a fresh graph from scanData"""
        self.plotWidget.clear()
        x = self.scanData.time
        y = self.scanData.trace
        self.plot = self.plotWidget.plot(x,y, pen=(255,0,0))

    def clearPlot(self):
        """clears all graphs from plot, after asking confermation"""
        reply = QMessageBox.question(self, 'Message',
            "Are you sure you want to clear the graph completely?", QMessageBox.Yes |
            QMessageBox.No, QMessageBox.No)

        if reply == QMessageBox.Yes:
            self.plotWidget.clear()

    def plotDataTest(self):
        """ plot a test curve in the plot widget"""
        x = np.arange(0,1000,1)
        noise = np.random.normal(0,1,1000)/1
        y = np.sin(x/10)+noise

        self.plot = self.plotWidget.plot(x,y, color='g')


#%% import export

    def importFile(self):
        '''import a single file form either .mat or .txt (csv) format'''
        self.scanData = rr.rrScan()
        filename = self.openFileNameDialog()
        self.scanData.importFile(filename)
        self.scanData.initParameters()
        self.fetchMetadata()
        self.plotScanData()

#        ext = os.path.splitext(filename)[-1].lower()
#        if ext == ".mat":
#            self.loadScanRaw(filename)
#        elif ext == ".txt":
#            self.loadScanCSV(filename)
#        else:
#            print("wrong file type, please try again")

    def loadScanCSV(self, filename):

        #filename = self.openFileNameDialog()
        self.scanData = rr.rrScan()
        self.scanData.importCSV(filename)
        self.scanData.initParameters()
        self.fetchMetadata()
        self.plotScanData()


    def loadScanRaw(self, filename):

        #filename = self.openFileNameDialog()
        self.scanData = rr.rrScan()
        self.scanData.importRawFile(filename)
        self.scanData.initParameters()
        print(self.scanData.parameters)
        #self.tree.addChild(self.scanData.parameters)
#        self.plotData()
        self.fetchMetadata()
        self.plotScanData()


    def saveasCSV(self):
        '''save object rrScan() to csv'''
        savedir = rr.getFolder()
        print(savedir)
        self.scanData.exportCSV(savedir)

#%% other
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