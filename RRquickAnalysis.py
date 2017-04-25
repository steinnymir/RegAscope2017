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
from matplotlib import cm

#%% Parameters

#Choose files
testfile = 'RuCl3-Pr-0.5mW-Pu-1.5mW-T-007.0k-1kAVG.mat'
testpath = '..//test_data'

fileRange = [0,0] # select files to use, [0,0] stands for all

# Filter:
filterOrder = 2
NyqCutFreq = 0.03

# Choose parameter being scanned
KeyDependence = 'Temperature'

Title = KeyDependence + ' Dependence'

#save = False
 
#import data
dataDict = rr.dir_to_dict(testpath, fileRange = fileRange)

#%% GUI plot
#win = pg.GraphicsWindow(title="Plotting test")
#win.resize(1000,600)
#win.setWindowTitle('Plotting test window title')
#gui1 = win.addPlot(title='Temperature Dependence')
##

#%% initialise regualr plot
fig = plt.figure(Title, figsize = (19,10))
plt.clf()
plt1 = fig.add_subplot(111)
plt1.set_xlabel('Time, ps', fontsize=18)
plt1.set_ylabel('Differential Reflectivity', fontsize=18)
plt1.set_title(Title,fontsize=26)
plt1.tick_params(axis='x', labelsize=12)
plt1.tick_params(axis='y', labelsize=12)
color=iter(cm.rainbow(np.linspace(0,1,len(dataDict))))


#%% analysis program
#just plots the curves you tell it, for now

dataDict = rr.dir_to_dict(testpath, fileRange = fileRange)

stack = 0
stackstep = 0.000 # to separate vertically the scans

# Plot a "key" dependence from 
for key in dataDict:
    stack += 1
    stackval = stack * stackstep
    x = rr.timezero_shift(dataDict[key]['data'][0], reverse = True)
    y = rr.quick_filter(rr.remove_DC_offset(dataDict[key]['data'][1])+stackval, 
                        order = filterOrder, cutfreq = NyqCutFreq)
#    gui1.plot(x,y)
    col = next(color)
    plt1.plot(x,y, label = dataDict[key][KeyDependence], c=col)

    
plt.show()
plt.legend(bbox_to_anchor=(1.005, 1), loc='best', borderaxespad=0., ncol=1)
    
    
    
    