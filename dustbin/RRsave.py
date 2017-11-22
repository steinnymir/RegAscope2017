# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 15:25:46 2017

@author: sagustss
"""

from lib import redred as rr
import os
import tkinter as tk
from tkinter import filedialog

""" Choose save settings """
#saveas = 'csv'
saveas = 'npy'
filtermultiplier = 0.01





#choose a folder

root = tk.Tk()
root.withdraw()
dataDir = filedialog.askdirectory(initialdir = 'E://DATA//_RAW')
seriesType = input('what tipe of dependence is this?')

name = os.path.basename(dataDir)

filtermultiplier = 0.01

rr.quickPlot(dataDir, cutfreq = filtermultiplier, KeyDependence = seriesType)
#    cmd = input('"y" for save or new filter multiplier')
#    if cmd == 'y':
#        save = True
#    else:
#        filtermultiplier == cmd




saveDir = filedialog.askdirectory(initialdir = 'E://DATA//RuCl3//')
os.chdir(saveDir)
if saveas == 'csv':
    for item in os.listdir(dataDir):
        filepath = dataDir + '//' + item
        if item != 't-cal.mat' and not os.path.isdir(filepath):
            filename = os.path.splitext(item)[0] + '.txt'
            dataDict = rr.file_to_dict(filepath)
            rr.save_filtered_trace(dataDict, saveDir, filename, filtermultiplier)
elif saveas == 'npy':
    dataDict = rr.dir_to_dict(dataDir)
    rr.save_obj(dataDict, name)
else:
    print('Format Error')

