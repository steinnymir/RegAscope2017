# -*- coding: utf-8 -*-
"""
Created on Thu May  4 16:35:21 2017

@author: vgrigore
"""

from functionlibrary import redred
import scipy
import os
import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog








sourceDir='E:\\2017-02-24\\'
#Filenames=os.listdir(sourceDir)
#Filenames = sortnames(Filenames)

#Names=[]

#i=0
#while i<37 :
   # Names.append('Mn2Au-MA125-30mWhor-20mWvert-RT-1kAV-MOKE-samplerot-'+str(0+i*10)+'degree.mat')
   # i+=1

Data=redred.RRscans()
#for item in Names:
 #   Data.addfilename(sourceDir + item)
Data.addfilenamesfromfolder(sourceDir)
Data.sortnames()
Data.importselectedfiles()
Data.definesampleorientforoldmokedata()
Data.filteritall()
Data.fliptimeall()
Data.fliptraceall()
Data.shiftTimeall(-80)
Data.initParametersall()

Data.rrPlot3d()