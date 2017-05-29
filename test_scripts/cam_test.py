# -*- coding: utf-8 -*-
"""
Created on Fri May 19 15:55:57 2017

@author: S.Y. Agustsson
"""

from __future__ import print_function

from VideoCapture import Device
import VideoCapture as vc
import cv2
import sys
from PIL import Image
from PIL.ImageQt import ImageQt
import pyqtgraph as pg
from PyQt5.QtGui import QImage

from PyQt5 import QtGui as qg
from PyQt5 import QtWidgets as qw
from PyQt5 import QtCore as qc

# app = qc.QCoreApplication.instance()
# if app is None:
#     app = qw.QApplication(sys.argv)
# # added initialization after first suggestion below
# qg.QImageReader.supportedImageFormats()
#
#
#
# im = Image.open('test.gif')
# image = ImageQt(im)
# pixmap = qg.QPixmap.fromImage(image)
# # pixmap = QtGui.QPixmap('test.gif')
#
# widget = qw.QWidget()
# hbox = qw.QHBoxLayout(widget)
# lbl = qw.QLabel(widget)
# lbl.setPixmap(pixmap)
# hbox.addWidget(lbl)
# widget.setLayout(hbox)
# widget.show()
# sys.exit(app.exec_())



def show_webcam(mirror=False):
    cam = cv2.VideoCapture(0)
    # print(type(cam.read()))
    while True:
        ret_val, img = cam.read()
        if mirror:
            img = cv2.flip(img, 1)
        pg.image(img)
        cv2.imshow('my webcam', img)
        if cv2.waitKey(1) == 27:
            break  # esc to quit

def plot_image():
    cam = cv2.VideoCapture(0)
    # print(type(cam.read()))
    ret_val, img = cam.read()

    pg.image(img)



plot_image()
while True:
    if cv2.waitKey(1) == 27:
        break  # esc to quit
cv2.destroyAllWindows()
# >>> print(im.format, im.size, im.mode)
#
# vc = VC('Test.avi')
# if vc.isOpened():
#     rval,frame = vc.read()
#
# image = Image.fromarray(frame)
#
# #
# image.show()
#
#
# del cam

