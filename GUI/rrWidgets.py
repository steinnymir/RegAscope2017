# -*- coding: utf-8 -*-
"""
Created on Thu May  4 16:48:43 2017

@author: S.Y. Agustsson
"""

from PyQt5 import QtGui as qg  # (the example applies equally well to PySide)
from PyQt5 import QtWidgets as qw
from PyQt5 import QtCore as qc
import sys
import pyqtgraph as pg
import numpy as np
from functionlibrary import redred as rr
import qdarkstyle
import os

run_local = True



class rrScanImport(qw.QWidget):
    """ """
    def __init__(self):
        super().__init__()
        self.title = 'the GUI test'
        self.left = 300
        self.top = 100
        self.width = 1400
        self.height = 900
        self.initUI()


    def initUI(self):
        """ Generate GUI layout """
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        self.makeLayout()
        self.show()


    def makeLayout(self):
        """ Generate the GUI layout """

        layout = qg.QGridLayout() # create a grid for subWidgets
        layout.setSpacing(10)
        self.setLayout(layout)

        font = qg.QFont()
        font.setBold(True)
        font.setPixelSize(15)


        # Define items
        self.importFileBtn = qw.QPushButton('Import File', self)
        self.importFileBtn.clicked.connect(self.importFile)
        layout.addWidget(self.importFileBtn, 18, 1, 1, 7)

        self.saveasCSVBtn = qw.QPushButton('Save as CSV', self)
        self.saveasCSVBtn.clicked.connect(self.saveasCSV)
        layout.addWidget(self.saveasCSVBtn, 16, 24, 1, 5)

        self.saveFigure = qw.QPushButton('Save Figure', self)
        #self.saveFigure.resize(self.clrButton.sizeHint())
        self.saveFigure.clicked.connect(self.nofunctionyet)
        layout.addWidget(self.saveFigure, 18, 24, 1, 5)

        self.nameTxtbox = qw.QLabel('Series Name:', self)
        self.nameTxtbox.setFont(font)
        layout.addWidget(self.nameTxtbox, 1, 1, 1, 7)
        self.nameTxtbox = qw.QLineEdit(self)
        self.nameTxtbox.setPlaceholderText('file name')
        self.nameTxtbox.editingFinished.connect(self.nofunctionyet)
        layout.addWidget(self.nameTxtbox, 2, 1, 1, 7)

        self.metadataTree_name = qw.QLabel('Metadata:', self)
#        self.metadataTree_name.setFont(font)
#        layout.addWidget(self.metadataTree_name, 3, 1)
        self.metadataTree = pg.DataTreeWidget()
        #self.metadataTree.setHeaderItem()
#        layout.addWidget(self.metadataTree, 4, 1, 13, 7)

        self.scanListTree_label = qw.QLabel('Scans:', self)
        self.scanListTree_label.setFont(font)
        layout.addWidget(self.scanListTree_label, 3, 1)
        self.scanListTree = pg.parametertree.ParameterTree()
        #self.metadataTree.setHeaderItem()
        layout.addWidget(self.scanListTree, 4, 1, 13, 7)



        self.plotWidget_name = qw.QLabel('Plot', self)
        self.plotWidget_name.setFont(font)
        layout.addWidget(self.plotWidget_name, 1, 9)
        self.setPlotWidget()
        layout.addWidget(self.plotWidget, 2, 9, 13, 20)

        # plot modification buttons
        self.modifyBox_name = qw.QLabel('Modify', self)
        layout.addWidget(self.modifyBox_name, 15, 10)
        self.modifyBox_name.setFont(font)

        self.flipXcb = qw.QPushButton('Flip x', self)
        self.flipXcb.clicked.connect(self.flipX)
        layout.addWidget(self.flipXcb, 16, 10)

        self.remDCcb = qw.QPushButton('Remove DC', self)
        self.remDCcb.clicked.connect(self.removeDC)
        layout.addWidget(self.remDCcb, 17, 10)

        self.shiftXcb = qw.QPushButton('Set Time Zero', self)
        self.shiftXcb.clicked.connect(self.setTimeZero)
        layout.addWidget(self.shiftXcb, 18, 10)

        self.timeZero = 0
        self.shiftTimeZero_input = pg.SpinBox(self)
        self.shiftTimeZero_input.setMinimumSize(1,25)
#        self.shiftTimeZero_input.setValidator(QtGui.QDoubleValidator())
        #self.shiftTimeZero_input.textChanged.connect(self.setShiftTimeZero)
        self.shiftTimeZero_input.valueChanged.connect(self.setTimeZero)
        layout.addWidget(self.shiftTimeZero_input, 18, 11, 1, 2)

        # Filter

#        self.filterBox = QGroupBox('filter', self)
#        layout.addWidget(self.filterBox, 15, 17, 3, 3)
#        self.makeFilterBox(layout=layout)
#
#    def makeFilterBox(self, layout):

        self.filterBox_name = qw.QLabel('Filter', self)
#        self.filterBox_name.setFont(font)
        layout.addWidget(self.filterBox_name, 15, 17)
        self.filterLowPassFreq_name = qw.QLabel('Frequency [THz]', self)
        layout.addWidget(self.filterLowPassFreq_name, 16,17)
        self.filterLowPassFreq = pg.SpinBox(self, dec=True)
        self.filterLowPassFreq.setMinimumSize(1,25)
        #self.filterLowPassFreq.setPlaceholderText('freq')
        self.filterLowPassFreq.valueChanged.connect(self.applyFilter)
        layout.addWidget(self.filterLowPassFreq, 16, 18, 1, 3)
#        self.applyFilterBtn = QPushButton('Apply', self)
#        self.applyFilterBtn.clicked.connect(self.filter_data_lowpass)
#        layout.addWidget(self.applyFilterBtn, 18, 17, 1, 5)

#%% slots
    @qc.pyqtSlot()
    def nofunctionyet(self):
        self.msg = qw.QMessageBox()
        self.msg.setIcon(qw.QMessageBox.Warning)
        self.msg.setText("No, this does nothing yet")
        self.msg.show()


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

    @qc.pyqtSlot()
    def shiftTimeZero(self):
        '''shift time of scan by timeZero, value given in the QLineEdit shift_time_scale'''
        txt = self.shiftTimeZero_input.text()
        num = float(txt)

        self.timeZero = num
        self.scanData.shiftTime(self.timeZero)

        self.plotScanData()


    @qc.pyqtSlot()
    def setTimeZero(self):
        '''set value given in shiftTimeZero_input() as new time zero'''
        txt = self.shiftTimeZero_input.text()
        num = float(txt)
        self.scanData.shiftTime(-self.timeZero)
        self.timeZero = num
        self.scanData.shiftTime(self.timeZero)
        self.plotScanData()


    @qc.pyqtSlot()
    def applyFilter(self):
        '''get filter frequency from textbox and apply the filter to a single scan'''
        freq = float(self.filterLowPassFreq.text())
        self.scanData.trace = self.scanData.rawtrace
        nyqfreq = self.scanData.nyqistFreq()

        if freq != 0 and freq < nyqfreq:

            cutfactor = freq / nyqfreq
            self.scanData.filterit(cutHigh=cutfactor)

        self.plotScanData()

    @qc.pyqtSlot()
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
        reply = qw.QMessageBox.question(self, 'Message',
            "Are you sure you want to clear the graph completely?", qw.QMessageBox.Yes |
            qw.QMessageBox.No, qw.QMessageBox.No)

        if reply == qw.QMessageBox.Yes:
            self.plotWidget.clear()

    def plotDataTest(self):
        """ plot a test curve in the plot widget"""
        x = np.arange(0,1000,1)
        noise = np.random.normal(0,1,1000)/1
        y = np.sin(x/10)+noise

        self.plot = self.plotWidget.plot(x,y, color='g')


#%% import export

    @qc.pyqtSlot()
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
        options = qw.QFileDialog.Options()
        options |= qw.QFileDialog.DontUseNativeDialog
        fileName, _ = qw.QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;Python Files (*.py)", options=options)
        if fileName:
            return(fileName)

    def set_statusbar(self):
        buttonReply =  qw.QMessageBox.question(self, 'PyQt5 message', "Do you like PyQt5?",
                                               qw.QMessageBox.Yes | qw.QMessageBox.No,
                                               qw.QMessageBox.No)
        if buttonReply == qw.QMessageBox.Yes:
            string = 'yes'
        else:
            string = 'no'
        self.statusBar().showMessage(string)



if __name__ == '__main__':

    if run_local:
    #run local script
        app = qc.QCoreApplication.instance()
        if app is None:
            app = qg.QApplication(sys.argv)

        prg = rrScanImport()
        prg.setStyleSheet(qdarkstyle.load_stylesheet_pyqt5())
        pg.setConfigOption('background', 0.1)
        pg.setConfigOption('foreground', 0.7)
        prg.show()
        app.exec_()

    else:
        #run master file script
        exec(open("../GUI_multiscan.py").read())