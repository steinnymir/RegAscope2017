# -*- coding: utf-8 -*-
"""
Created on Fri May 19 15:55:57 2017

@author: S.Y. Agustsson
"""


from VideoCapture import Device
import pyqtgraph as pg
from PIL import Image, ImageQt




cam = Device(1)
ImageQt(cam.getImage())

#im = Image.open(image)
image.show()


del cam

