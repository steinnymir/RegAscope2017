import numpy
import cv2
from PyQt5 import QtGui as qg, QtWidgets as qw, QtCore as qc, uic
import pyqtgraph as pg
from PIL import Image

import sys

def main():
    app = qw.QApplication(sys.argv)
    w = CamView()
    w.resize(600, 400)
    w.show()
    app.exec_()

def main_():
    show_webcam()


class CamView(qw.QMainWindow):
    CAMERA = 0

    def __init__(self):
        super().__init__()

        uic.loadUi('cam_gui.ui', self)
        self.cam = cv2.VideoCapture(self.CAMERA)

    def run_video(self):

        vid = cv2.VideoCapture(self.CAMERA)
        ret, frame = vid.read()
        while ret:
            # Qimg = convert(frame)
            self.label.setpixmap(Qimg)
            self.label.update()
            ret, frame = vid.read()
            qw.QApplication.processEvents()


        # vb = pg.ViewBox()
        # self.graphicsView.setCentralItem(vb)
        # vb.setAspectLocked()
        # img = pg.ImageItem()
        # vb.addItem(img)
        # vb.setRange(qc.QRectF(0, 0, 512, 512))

def plot_image():
    cam = cv2.VideoCapture(0)
    # print(type(cam.read()))
    ret_val, img = cam.read()

def show_webcam():
    cam = cv2.VideoCapture(1)
    while True:
        ret_val, img = cam.read()

        print(ret_val)
        cv2.imshow('logitech webcam', img)
        if cv2.waitKey(1) == 27:
            break



if __name__=="__main__":

    main()