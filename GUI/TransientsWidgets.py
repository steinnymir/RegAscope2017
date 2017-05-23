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
from functionlibrary import transient as tr
import qdarkstyle
import os

run_local = True

class MainWindow(qw.QMainWindow):
    """ Main application window """

    def __init__(self):
        '''Initialize by setting window title, size and graphic options'''
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
        '''Create the layout, adding central widget, layout style and status
        bar. '''
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        layout = qg.QGridLayout() # create a grid for subWidgets
        layout.setSpacing(10)
        self.setLayout(layout)

        self.centralWidget = TransientAnalysisWidget()
        layout.addWidget(self.centralWidget,0,0)
        self.setCentralWidget(self.centralWidget)
        self.statusBar().showMessage('Message in statusbar.')

        self.show()


class TransientAnalysisWidget(qw.QWidget):
    """ """
    def __init__(self):
        super().__init__()
        self.title = 'the GUI test'
#        self.left = 300
#        self.top = 100
#        self.width = 1400
#        self.height = 900
        self.initUI()

        self.data = tr.Transients()
        self.data_memory = {}


    def initUI(self):
        """ Generate GUI layout """
        self.setWindowTitle(self.title)
#        self.setGeometry(self.left, self.top, self.width, self.height)

        self.makeLayout()
        self.show()


    def makeLayout(self):
        """ Generate the GUI layout """

        layout = qg.QGridLayout() # create a grid for subWidgets
        layout.setSpacing(10)
        self.setLayout(layout)

        self.setFontStyles()


        # -------- Define items ----------
        # import button
        self.importFileBtn = qw.QPushButton('Import File(s)', self)
        self.importFileBtn.clicked.connect(self.importTransients_fromFile)
        layout.addWidget(self.importFileBtn, 18, 1, 1, 7)

        # Save block
        self.saveasCSVBtn = qw.QPushButton('Save as CSV', self)
        self.saveasCSVBtn.clicked.connect(self.saveasCSV)
        layout.addWidget(self.saveasCSVBtn, 16, 24, 1, 5)
        self.saveFigure = qw.QPushButton('Save Figure', self)
        #self.saveFigure.resize(self.clrButton.sizeHint())
        self.saveFigure.clicked.connect(self.nofunctionyet)
        layout.addWidget(self.saveFigure, 18, 24, 1, 5)

        # File name block
        self.nameTxtbox = qw.QLabel('Series Name:', self)
        self.nameTxtbox.setFont(self.title_font)
        layout.addWidget(self.nameTxtbox, 1, 1, 1, 7)
        self.nameTxtbox = qw.QLineEdit(self)
        self.nameTxtbox.setPlaceholderText('file name')
        self.nameTxtbox.editingFinished.connect(self.nofunctionyet)
        self.nameTxtbox.returnPressed.connect(self.nofunctionyet)
        layout.addWidget(self.nameTxtbox, 2, 1, 1, 7)

        # Metadata tree
#        self.metadataTree_name = qw.QLabel('Metadata:', self)
#        self.metadataTree_name.setFont(font)
#        layout.addWidget(self.metadataTree_name, 3, 1)
#        self.metadataTree = pg.DataTreeWidget()
        #self.metadataTree.setHeaderItem()
#        layout.addWidget(self.metadataTree, 4, 1, 13, 7)

        # Scan list tree
        self.scanListTree_label = qw.QLabel('Scans:', self)
        self.scanListTree_label.setFont(self.title_font)
        layout.addWidget(self.scanListTree_label, 3, 1)
        self.transientData_list = qw.QListWidget()
#        self.scanListTree = pg.parametertree.ParameterTree()
        #self.metadataTree.setHeaderItem()
        layout.addWidget(self.transientData_list, 4, 1, 13, 7)

        # Plot widget
        self.plotWidget_name = qw.QLabel('Plot', self)
        self.plotWidget_name.setFont(self.title_font)
        layout.addWidget(self.plotWidget_name, 1, 9)

        self.setPlotWidget()
        layout.addWidget(self.plotWidget, 2, 9, 13, 20)


        # plot modification buttons
#        self.DataAnalysisBox_label = qw.QLabel('Modify', self)
#        self.DataAnalysisBox_label.setFont(font)
#        layout.addWidget(self.DataAnalysisBox_label, 15, 10)
        self.DataAnalysisBox = qw.QGroupBox()
        self.setDataAnalysisBox()
        layout.addWidget(self.DataAnalysisBox, 16,10)

#########################################################
#########################################################
#########################################################
#########################################################




#%% slots
    @qc.pyqtSlot()
    def nofunctionyet(self):
        self.msg = qw.QMessageBox()
        self.msg.setIcon(qw.QMessageBox.Warning)
        self.msg.setText("No, this does nothing yet")
        #self.statusBar().showMessage('Message in statusbar.')
        # MainWindow.statusBar().showMessage('pushed the wrong button')
        self.msg.show()

    def setPlotWidget(self):
        ''' Create the widget for plotting scans and define it's
        properties'''
        pg.setConfigOptions(antialias=True)
        self.plotWidget = pg.PlotWidget()

#        pg.setConfigOption('background', 'w')
#        pg.setConfigOption('foreground', 'k')

    def setDataAnalysisBox(self):
        ''' '''
        self.DataAnalysisBox_layout = qw.QGridLayout()
        self.DataAnalysisBox.setLayout(self.DataAnalysisBox_layout)

        self.DataModification_label = qw.QLabel('Modifications:', self)
        self.DataModification_label.setFont(self.subtitle_font)
        self.DataAnalysisBox_layout.addWidget(self.DataModification_label,
                                              0, 0)


        self.flipX_cb = qw.QCheckBox(self.DataAnalysisBox)
        self.flipX_cb.setText('Flip x')
        self.flipX_cb.setChecked(True)
#        self.flipX_checkBox.clicked.connect(self.flipX)
        self.DataAnalysisBox_layout.addWidget(self.flipX_cb, 1, 0)

        self.removeDC_cb = qw.QCheckBox(self.DataAnalysisBox)
        self.removeDC_cb.setText('Remove DC Offset')
        # self.removeDC_cb.clicked.connect(self.removeDC)
        self.removeDC_cb.setChecked(True)
        self.DataAnalysisBox_layout.addWidget(self.removeDC_cb, 2, 0)

        self.setTimeZero_cb = qw.QCheckBox(self.DataAnalysisBox)
        self.setTimeZero_cb.setText('Set time Zero')
        # self.setTimeZero_cb.clicked.connect(self.setTimeZero)
        self.setTimeZero_cb.setChecked(True)
        self.DataAnalysisBox_layout.addWidget(self.setTimeZero_cb, 3, 0)

        self.setTimeZero_sb = pg.SpinBox(self.DataAnalysisBox)
        self.setTimeZero_sb.setMinimumSize(1,25)
        #self.shiftTimeZero_input.valueChanged.connect(self.setTimeZero)
        self.DataAnalysisBox_layout.addWidget(self.setTimeZero_sb, 2, 1)
        #self.shiftTimeZero_input.setValidator(QtGui.QDoubleValidator())
        #self.shiftTimeZero_input.textChanged.connect(self.setShiftTimeZero)

        # Filter
        self.timeZero = 0
        self.filter_label = qw.QLabel('Filter [THz]', self.DataAnalysisBox)
        self.filter_label.setFont(self.subtitle_font)
#        self.filterBox_name.setFont(font)
        self.DataAnalysisBox_layout.addWidget(self.filter_label, 0,2)

#        self.filterLowPass_label = qw.QLabel('Low Pass',
#                                                self.DataAnalysisBox)

        self.filterLowPass_cb = qw.QCheckBox(self.DataAnalysisBox)
        self.filterLowPass_cb.setText('Low Pass Frequency')
        # self.setTimeZero_cb.clicked.connect(self.setTimeZero)
        self.filterLowPass_cb.setChecked(True)
        self.DataAnalysisBox_layout.addWidget(self.filterLowPass_cb, 1,2)

        self.filterLowPass_sb = pg.SpinBox(self, dec=True)
        self.filterLowPass_sb.setMinimumSize(1,25)
        #self.filterLowPassFreq.setPlaceholderText('freq')
        self.filterLowPass_sb.valueChanged.connect(self.applyFilter)
        self.DataAnalysisBox_layout.addWidget(self.filterLowPass_sb, 1,3)


    def setFontStyles(self):
        ''' Give settings for fonts to use in widget'''
        self.title_font = qg.QFont()
        self.subtitle_font = qg.QFont()
        self.text_font = qg.QFont()

        self.title_font.setBold(True)
        self.title_font.setPixelSize(15)

        self.subtitle_font.setPixelSize(12)
        self.subtitle_font.setBold(True)

        self.text_font.setPixelSize(10)


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

    def importTransients_fromFile(self):
        ''' Import multiple files to analysis program.
        for now only overwrites, adding append function soon.
        '''
        filename = self.openFileNameDialog()
        append = False
        self.data.importFiles(filename, append)
        print(self.data)
        self.refresh_transient_list()
        self.plot_transients_all()

    def refresh_transient_list(self):
        ''' refresh list of transients reported in list_widget'''
#        for n, transient in enumerate(self.data):
#            self.transientData_list.clear()
        self.transientData_list.addItem('test')


    def plot_transients_all(self):
        ''' '''
        for scan in self.data:
            x= scan.time
            y= scan.trace
            self.plot = self.plotWidget.plot(x,y, pen=(255,0,0))


    def plot_transients_selected(self):
        ''' '''
        pass

    def plotScanData(self):
        """ clears the graph and plots a fresh graph from scanData"""
        self.plotWidget.clear()
        x = self.scanData.time
        y = self.scanData.trace
        self.plot = self.plotWidget.plot(x,y, pen=(255,0,0))

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






















    def fetchMetadata(self):
        '''retrieve metadata from rrScan() and write it in metadataTree
                                        Deprecated
        '''
        self.metadataTree.setData(self.scanData.fetchMetadata(), hideRoot=True)
        print(self.metadataTree.headerItem())



#%% Scan Modifcations - old

    @qc.pyqtSlot()
    def flipX(self):
        self.scanData.flipTime()
        self.plotScanData()

    @qc.pyqtSlot()
    def flipY(self):
        self.scanData.flipTrace()
        self.plotScanData()

    @qc.pyqtSlot()
    def shiftTimeZero(self):
        '''shift time of scan by timeZero, value given in the QLineEdit shiftTimeZero'''
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

        prg = MainWindow()
        prg.setStyleSheet(qdarkstyle.load_stylesheet_pyqt5())
        pg.setConfigOption('background', 0.1)
        pg.setConfigOption('foreground', 0.7)
        prg.show()
        app.exec_()

    else:
        #run master file script
        exec(open("../TransientsAnalysis.py").read())