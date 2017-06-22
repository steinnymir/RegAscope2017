import numpy
import cv2
from PyQt5 import QtGui as qg, QtWidgets as qw, QtCore as qc, uic
import pyqtgraph as pg
from PIL import Image
from PIL.ImageQt import ImageQt

import sys

def main():
    app = qw.QApplication(sys.argv)
    w = CamView()
    w.resize(600, 400)
    w.show()
    app.exec_()



class CamView(qw.QMainWindow):

    def __init__(self):
        super().__init__()

        # uic.loadUi('cam_gui.ui', self)

        self.layout = qw.QGridLayout()  # create a grid for subWidgets
        self.layout.setSpacing(10)
        self.setLayout(self.layout)
        self.centralWidget = videoWidget()
        self.layout.addWidget(self.centralWidget, 0, 0)
        self.setCentralWidget(self.centralWidget)
        self.statusBar().showMessage('Message in statusbar.')

        self.show()

class videoWidget(qw.QWidget):
    CAMERA = 1

    def __init__(self):
        super(videoWidget, self).__init__()

        layout = qw.QGridLayout()  # create a grid for subWidgets
        layout.setSpacing(10)
        self.setLayout(layout)

        self.camWindow = pg.ImageView()
        layout.addWidget(self.camWindow, 0,0,1,1)

        self.cam = cv2.VideoCapture(self.CAMERA)
        ret, frame = self.cam.read()
        frame = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)

        self.testImage = cv2.imread('mini.jpg',0)
        print((self.testImage))


        img = self.testImage
        img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
        self.show()

        self.run_video()
        # self.camWindow.setImage(frame)




    def run_video(self):

        vid = cv2.VideoCapture(self.CAMERA)
        ret, frame = vid.read()
        while ret:
            Qimg = self.convert_frame(frame)
            self.camWindow.setImage(Qimg)
            ret, frame = vid.read()
            qw.QApplication.processEvents()
            if cv2.waitKey(1) == 27:
                break
    def convert_frame(self,img):
        return cv2.cvtColor(img, cv2.COLOR_BGR2RGB)


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