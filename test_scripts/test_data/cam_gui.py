import numpy
import cv2
from PyQt5 import QtGui as qg
from PyQt5 import QtWidgets as qw
from PyQt5 import QtCore as qc

import sys

class CamView(qw.QWidget):
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
        layout = qg.QGridLayout()  # create a grid for subWidgets
        layout.setSpacing(10)
        self.setLayout(layout)

        font = qg.QFont()
        font.setBold(True)
        font.setPixelSize(15)

        self.videoScreen = pg.

def plot_image():
    cam = cv2.VideoCapture(0)
    # print(type(cam.read()))
    ret_val, img = cam.read()


if __name__=="__main__":
    app = qw.QApplication(sys.argv)
    w = CamView()
    w.resize(600, 400)
    w.show()
    app.exec_()