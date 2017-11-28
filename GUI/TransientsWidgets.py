# -*- coding: utf-8 -*-
"""
Created on Thu May  4 16:48:43 2017

@author: S.Y. Agustsson
"""

from PyQt5 import QtGui as QG
from PyQt5 import QtWidgets as QW
from PyQt5 import QtCore as QC
import sys
import pyqtgraph as pg
import numpy as np
from lib import transient as tr
import qdarkstyle
import os


class MainWindow(QW.QMainWindow):
    """ Main application window """

    def __init__(self):
        """Initialize by setting window title, size and graphic options"""
        super().__init__()

        self.title = 'Transient Analyser'
        self.left = 300
        self.top = 100
        self.width = 1400
        self.height = 900
        self.initUI()

        # set the cool dark theme
        self.setStyleSheet(qdarkstyle.load_stylesheet_pyqt5())
        pg.setConfigOption('background', 0.1)
        pg.setConfigOption('foreground', 0.7)

    def initUI(self):
        """Create the layout, adding central widget, layout style and status
        bar. """
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        layout = QG.QGridLayout()  # create a grid for subWidgets
        layout.setSpacing(10)
        self.setLayout(layout)


        self.centralWidget = TransientAnalysisWidget()
        layout.addWidget(self.centralWidget, 0, 0)
        self.setCentralWidget(self.centralWidget)
        self.statusBar().showMessage('Message in statusbar.')

        self.show()


class TransientAnalysisWidget(QW.QWidget):
    """ """

    def __init__(self):
        super().__init__()
        self.title = 'the GUI test'
        #        self.left = 300
        #        self.top = 100
        #        self.width = 1400
        #        self.height = 900
        self.initUI()

        self.data = tr.MultiTransients()
        self.data_memory = {}

    def initUI(self):
        """ Generate GUI layout """
        self.setWindowTitle(self.title)
        #        self.setGeometry(self.left, self.top, self.width, self.height)

        self.make_layout()
        self.show()

    def make_layout(self):
        """ Generate the GUI layout """

        layout = QG.QGridLayout()  # create a grid for subWidgets
        layout.setSpacing(10)
        self.setLayout(layout)

        self.setup_font_styles()

        # -------- Define items ----------
        # import button
        self.importFileBtn = QW.QPushButton('Import File(s)', self)
        self.importFileBtn.clicked.connect(self.import_multiple_files)
        layout.addWidget(self.importFileBtn, 18, 1, 1, 7)

        # Save block
        self.saveasCSVBtn = QW.QPushButton('Save as CSV', self)
        self.saveasCSVBtn.clicked.connect(self.saveasCSV)
        layout.addWidget(self.saveasCSVBtn, 16, 24, 1, 5)
        self.saveFigure = QW.QPushButton('Save Figure', self)
        # self.saveFigure.resize(self.clrButton.sizeHint())
        self.saveFigure.clicked.connect(self.no_function_yet)
        layout.addWidget(self.saveFigure, 18, 24, 1, 5)

        # File name block
        self.nameTxtbox = QW.QLabel('Series Name:', self)
        self.nameTxtbox.setFont(self.title_font)
        layout.addWidget(self.nameTxtbox, 1, 1, 1, 7)
        self.nameTxtbox = QW.QLineEdit(self)
        self.nameTxtbox.setPlaceholderText('file name')
        self.nameTxtbox.editingFinished.connect(self.no_function_yet)
        self.nameTxtbox.returnPressed.connect(self.no_function_yet)
        layout.addWidget(self.nameTxtbox, 2, 1, 1, 7)

        # Metadata tree
        #        self.metadataTree_name = qw.QLabel('Metadata:', self)
        #        self.metadataTree_name.setFont(font)
        #        layout.addWidget(self.metadataTree_name, 3, 1)
        #        self.metadataTree = pg.DataTreeWidget()
        # self.metadataTree.setHeaderItem()
        #        layout.addWidget(self.metadataTree, 4, 1, 13, 7)

        # Scan list tree
        self.scanListTree_label = QW.QLabel('Scans:', self)
        self.scanListTree_label.setFont(self.title_font)
        layout.addWidget(self.scanListTree_label, 3, 1)
        self.transientData_list = QW.QListWidget()
        #        self.scanListTree = pg.parametertree.ParameterTree()
        # self.metadataTree.setHeaderItem()
        layout.addWidget(self.transientData_list, 4, 1, 13, 7)

        # Plot widget
        self.plotWidget_name = QW.QLabel('Plot', self)
        self.plotWidget_name.setFont(self.title_font)
        layout.addWidget(self.plotWidget_name, 1, 9)

        self.setup_plot_widget()
        layout.addWidget(self.plotWidget, 2, 9, 13, 20)

        # plot modification buttons
        #        self.DataAnalysisBox_label = qw.QLabel('Modify', self)
        #        self.DataAnalysisBox_label.setFont(font)
        #        layout.addWidget(self.DataAnalysisBox_label, 15, 10)
        self.DataAnalysisBox = QW.QGroupBox()
        self.setup_data_analysis_box()
        layout.addWidget(self.DataAnalysisBox, 16, 10)

    #########################################################
    #########################################################
    #########################################################
    #########################################################

    # %% slots
    @QC.pyqtSlot()
    def no_function_yet(self):
        self.msg = QW.QMessageBox()
        self.msg.setIcon(QW.QMessageBox.Warning)
        self.msg.setText("No, this does nothing yet")
        # self.statusBar().showMessage('Message in statusbar.')
        # MainWindow.statusBar().showMessage('pushed the wrong button')
        self.msg.show()

    def setup_plot_widget(self):
        """ Create the widget for plotting scans and define it's
        properties"""
        pg.setConfigOptions(antialias=True)
        self.plotWidget = pg.PlotWidget()

    #        pg.setConfigOption('background', 'w')
    #        pg.setConfigOption('foreground', 'k')

    def setup_data_analysis_box(self):
        ''' '''
        self.DataAnalysisBox_layout = QW.QGridLayout()
        self.DataAnalysisBox.setLayout(self.DataAnalysisBox_layout)

        self.DataModification_label = QW.QLabel('Modifications:', self)
        self.DataModification_label.setFont(self.subtitle_font)
        self.DataAnalysisBox_layout.addWidget(self.DataModification_label,
                                              0, 0)

        self.flipX_cb = QW.QCheckBox(self.DataAnalysisBox)
        self.flipX_cb.setText('Flip x')
        self.flipX_cb.setChecked(True)
        #        self.flipX_checkBox.clicked.connect(self.flip_time_scale)
        self.DataAnalysisBox_layout.addWidget(self.flipX_cb, 1, 0)

        self.removeDC_cb = QW.QCheckBox(self.DataAnalysisBox)
        self.removeDC_cb.setText('Remove DC Offset')
        # self.removeDC_cb.clicked.connect(self.remove_DC_offset)
        self.removeDC_cb.setChecked(True)
        self.DataAnalysisBox_layout.addWidget(self.removeDC_cb, 2, 0)

        self.setTimeZero_cb = QW.QCheckBox(self.DataAnalysisBox)
        self.setTimeZero_cb.setText('Set time Zero')
        # self.setTimeZero_cb.clicked.connect(self.set_time_zero)
        self.setTimeZero_cb.setChecked(True)
        self.DataAnalysisBox_layout.addWidget(self.setTimeZero_cb, 3, 0)

        self.setTimeZero_sb = pg.SpinBox(self.DataAnalysisBox)
        self.setTimeZero_sb.setMinimumSize(1, 25)
        # self.shiftTimeZero_input.valueChanged.connect(self.set_time_zero)
        self.DataAnalysisBox_layout.addWidget(self.setTimeZero_sb, 2, 1)
        # self.shiftTimeZero_input.setValidator(QtGui.QDoubleValidator())
        # self.shiftTimeZero_input.textChanged.connect(self.setShiftTimeZero)

        # Filter
        # self.timeZero = 0
        self.filter_label = QW.QLabel('Filter [THz]', self.DataAnalysisBox)
        self.filter_label.setFont(self.subtitle_font)
        #        self.filterBox_name.setFont(font)
        self.DataAnalysisBox_layout.addWidget(self.filter_label, 0, 2)

        #        self.filterLowPass_label = qw.QLabel('Low Pass', self.DataAnalysisBox)

        self.filterLowPass_cb = QW.QCheckBox(self.DataAnalysisBox)
        self.filterLowPass_cb.setText('Low Pass Frequency')
        # self.setTimeZero_cb.clicked.connect(self.set_time_zero)
        self.filterLowPass_cb.setChecked(True)
        self.DataAnalysisBox_layout.addWidget(self.filterLowPass_cb, 1, 2)

        self.filterLowPass_sb = pg.SpinBox(self, dec=True)
        self.filterLowPass_sb.setMinimumSize(1, 25)
        # self.filterLowPassFreq.setPlaceholderText('freq')
        self.filterLowPass_sb.valueChanged.connect(self.filter_data_lowpass)
        self.DataAnalysisBox_layout.addWidget(self.filterLowPass_sb, 1, 3)

    def setup_font_styles(self):
        """ Give settings for fonts to use in widget"""
        self.title_font = QG.QFont()
        self.subtitle_font = QG.QFont()
        self.text_font = QG.QFont()

        self.title_font.setBold(True)
        self.title_font.setPixelSize(15)

        self.subtitle_font.setPixelSize(12)
        self.subtitle_font.setBold(True)

        self.text_font.setPixelSize(10)

    # %% import export

    @QC.pyqtSlot()
    def import_single_file(self):  # todo: translate for new Transient() class
        """ Import a single file form either .mat or .txt (csv) file."""
        self.scanData = tr.Transient()
        filename = self.openFileNameDialog()  # choose the file to import
        self.scanData.import_file(filename)
        self.scanData.initParameters()
        self.plotScanData()

    #        ext = os.path.splitext(filename)[-1].lower()
    #        if ext == ".mat":
    #            self.loadScanRaw(filename)
    #        elif ext == ".txt":
    #            self.loadScanCSV(filename)
    #        else:
    #            print("wrong file type, please try again")

    def import_multiple_files(self):  # todo: Improve! Rethink importing method
        """
        Import multiple files to analysis program.
        for now only overwrites, adding append function soon.
        """
        filename = self.openFileNameDialog()
        append = False
        self.data.import_files(filename, append)
        print(self.data)
        self.refresh_transient_list()
        self.plot_all_transients()

    def refresh_transient_list(self):
        """ refresh list of transients reported in list_widget"""
        #        for n, transient in enumerate(self.data):
        #            self.transientData_list.clear()
        self.transientData_list.addItem('test')

    def plot_all_transients(self):
        """ """
        for scan in self.data:
            x = scan.time
            y = scan.trace
            self.plot = self.plotWidget.plot(x, y, pen=(255, 0, 0))

    def plot_selected_transients(self):  # todo: write method
        """ """
        pass

    def plotScanData(self):  # todo: translate for TransientsSet()
        """ clears the graph and plots a fresh graph from scanData"""
        self.plotWidget.clear()
        x = self.scanData.time
        y = self.scanData.trace
        self.plot = self.plotWidget.plot(x, y, pen=(255, 0, 0))

    def loadScanCSV(self, filename):  # todo: translate for TransientsSet()

        # filename = self.openFileNameDialog()
        self.scanData = rr.rrScan()
        self.scanData.import_file_csv(filename)
        self.scanData.initParameters()
        self.plotScanData()

    def loadScanRaw(self, filename):  # todo: translate for TransientsSet()

        # filename = self.openFileNameDialog()
        self.scanData = rr.rrScan()
        self.scanData.importRawFile(filename)
        self.scanData.initParameters()
        print(self.scanData.parameters)
        # self.tree.addChild(self.scanData.parameters)
        #        self.plotData()
        self.plotScanData()

    def saveasCSV(self):  # todo: translate for TransientsSet()
        """save object rrScan() to csv"""
        savedir = rr.getFolder()
        print(savedir)
        self.scanData.export_file_csv(savedir)

    # %% Scan Modifcations - old

    @QC.pyqtSlot()
    def flip_time_scale(self):  # todo: translate for TransientsSet(), is it really useful?
        self.scanData.flip_time()
        self.plotScanData()

    @QC.pyqtSlot()
    def flip_trace(self):  # todo: translate for TransientsSet(), is it really useful?
        """ Flip data on y scale: f(x) = -f(x).
        Replot data"""
        self.scanData.flip_trace()
        self.plotScanData()

    @QC.pyqtSlot()
    def shift_time_scale(self):  # todo: translate for TransientsSet(), is it really useful?
        """shift time of scan by timeZero, value given in the QLineEdit shift_time_scale"""
        txt = self.shiftTimeZero_input.text()
        num = float(txt)
        self.timeZero = num
        self.scanData.shift_time(self.timeZero)
        self.plotScanData()

    @QC.pyqtSlot()
    def set_time_zero(self):  # todo: translate for TransientsSet(), is it really useful?
        """set value given in shiftTimeZero_input() as new time zero"""
        txt = self.shiftTimeZero_input.text()
        num = float(txt)
        self.scanData.shift_time(-self.timeZero)
        self.timeZero = num
        self.scanData.shift_time(self.timeZero)
        self.plotScanData()

    @QC.pyqtSlot()
    def filter_data_lowpass(self):  # todo: translate for TransientsSet(), is it really useful?
        """get filter frequency from textbox and apply the filter to a single scan"""
        freq = float(self.filterLowPassFreq.text())
        self.scanData.trace = self.scanData.rawtrace
        nyqfreq = self.scanData.nyqistFreq()

        if freq != 0 and freq < nyqfreq:
            cutfactor = freq / nyqfreq
            self.scanData.filter_low_pass(cutHigh=cutfactor)

        self.plotScanData()

    @QC.pyqtSlot()
    def remove_DC_offset(self):  # todo: translate for TransientsSet(), is it really useful?
        self.scanData.remove_DC_offset()
        self.plotScanData()

    # %% Plots

    def plotScanData(self):  # todo: translate for TransientsSet(), is it really useful?
        """ clears the graph and plots a fresh graph from scanData"""
        self.plotWidget.clear()
        x = self.scanData.time
        y = self.scanData.trace
        self.plot = self.plotWidget.plot(x, y, pen=(255, 0, 0))

    def clearPlot(self):  # todo: translate for TransientsSet(), is it really useful?
        """clears all graphs from plot, after asking confermation"""
        reply = QW.QMessageBox.question(self, 'Message',
                                        "Are you sure you want to clear the graph completely?", QW.QMessageBox.Yes |
                                        QW.QMessageBox.No, QW.QMessageBox.No)

        if reply == QW.QMessageBox.Yes:
            self.plotWidget.clear()

    def plotDataTest(self):  # todo: translate for TransientsSet(), is it really useful?
        """ plot a test curve in the plot widget"""
        x = np.arange(0, 1000, 1)
        noise = np.random.normal(0, 1, 1000) / 1
        y = np.sin(x / 10) + noise

        self.plot = self.plotWidget.plot(x, y, color='g')

    # %% import export


    def loadScanCSV(self, filename):

        # filename = self.openFileNameDialog()
        self.scanData = rr.rrScan()
        self.scanData.import_file_csv(filename)
        self.scanData.initParameters()
        self.plotScanData()

    def loadScanRaw(self, filename):

        # filename = self.openFileNameDialog()
        self.scanData = rr.rrScan()
        self.scanData.importRawFile(filename)
        self.scanData.initParameters()
        print(self.scanData.parameters)
        # self.tree.addChild(self.scanData.parameters)
        #        self.plotData()
        self.plotScanData()

    def saveasCSV(self):
        """save object rrScan() to csv"""
        savedir = rr.getFolder()
        print(savedir)
        self.scanData.export_file_csv(savedir)

    # %% other
    def openFileNameDialog(self):
        options = QW.QFileDialog.Options()
        options |= QW.QFileDialog.DontUseNativeDialog
        fileName, _ = QW.QFileDialog.getOpenFileName(self, "QFileDialog.getOpenFileName()", "",
                                                     "All Files (*);;Python Files (*.py)", options=options)
        if fileName:
            return (fileName)

    def set_statusbar(self):
        buttonReply = QW.QMessageBox.question(self, 'PyQt5 message', "Do you like PyQt5?",
                                              QW.QMessageBox.Yes | QW.QMessageBox.No,
                                              QW.QMessageBox.No)
        if buttonReply == QW.QMessageBox.Yes:
            string = 'yes'
        else:
            string = 'no'
        self.statusBar().showMessage(string)


if __name__ == '__main__':



    app = QC.QCoreApplication.instance()
    if app is None:
        app = QG.QApplication(sys.argv)

    prg = MainWindow()
    prg.setStyleSheet(qdarkstyle.load_stylesheet_pyqt5())
    pg.setConfigOption('background', 0.1)
    pg.setConfigOption('foreground', 0.7)
    prg.show()
    app.exec_()

