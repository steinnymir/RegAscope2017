# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 13:45:12 2017

@author: S.Y. Agustsson
"""

from functionlibrary import redred as rr

#import os
#import numpy as np
#import scipy as sp
#GUI
#from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg
#plot
import matplotlib.pyplot as plt
#from matplotlib import cm

#%% Parameters

#Choose files
testfile = 'RuCl3-Pr-0.5mW-Pu-1.5mW-T-007.0k-1kAVG.mat'
testpath = 'c://Users//sagustss//py_code//test_data'
fileRange = [0,0] # select files to use, [0,0] stands for all

# Filter:
filterOrder = 2
NyqCutFreq = 0.01

#save = False

win = pg.GraphicsWindow(title="Plotting test")
win.resize(1000,600)
win.setWindowTitle('Plotting test window title')
gui1 = win.addPlot(title='Temperature Dependence')



fig = plt.figure("test")
plt.clf()
plt1 = fig.add_subplot(111)


dataDict = rr.dir_to_dict(testpath, fileRange = [0,5])

stack = 0
stackstep = 0.000 # to separate vertically the scans
for key in dataDict:
    stack += 1
    stackval = stack * stackstep
    print(dataDict[key]['Temperature'])
    x = rr.timezero_shift(dataDict[key]['data'][0], reverse = True)
    y = rr.quick_filter(rr.remove_DC_offset(dataDict[key]['data'][1])+stackval, 
                        order = filterOrder, cutfreq = NyqCutFreq)
    gui1.plot(x,y)
    plt1.plot(x,y)  
    

    
    
    
    