# -*- coding: utf-8 -*-
"""
Created on Wed May  3 21:38:01 2017

@author: sagustss
"""

from lib import redred as rr
import os
import numpy as np
import scipy
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog
from matplotlib import cm


root = tk.Tk()
root.withdraw()

baseDir = 'E://DATA//_RAW//RuCl3 - Copy//Power dependence Low - Low Temperature'

filename = filedialog.askopenfilename(initialdir = baseDir)


scan = rr.rrScan() 
scan.importRawFile(filename)

scan.filterit(cutHigh = 0.01)
scan.flipTime()
scan.shiftTime(-50)
scan.removeDC()
#scan.flip_trace()


Title = os.path.basename(baseDir)
plt.ion()
plt.figure(Title)
#plt1 = fig.add_subplot(111)
plt.plot(scan.time,scan.trace)
plt.show()
print(scan.parameters)

save = false
if save:
    saveDir = 'E://DATA//RuCl3//Analysis//' + Title + '//'
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)
    scan.exportCSV(saveDir)
    print('saved')