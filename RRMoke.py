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







class RRscans(object):
    '''Contains multiple scans and methods to sort them and plot different grafs'''
    
    def __init__(self):
        '''initilize list of scans and atributes'''
        self.scans=[]#list of redred scans
        self.samplename='samplename'
        self.filenames=[]#list of the files to read
        
#%%import matlab files functions        
    def addfilename(self, fullname):
        '''add single filename to the list of the files to read'''
        self.filenames.append(fullname)
    
    def addfilenamesfromfolder(self, directory):
        '''add filenames from selected folder to the list of the files to read'''
        newnames= os.listdir(directory)
        for item in newnames:
            self.addfilename(directory + item)
    
    def choosefile(self):
        '''open dialog window to choose file to be added to the list of the files to read'''
        root = tk.Tk()
        root.withdraw()
        filenames = filedialog.askopenfilenames()
        for item in filenames:
            self.addfilename(item)
        
    def choosefilesfromfolder(self):
        '''open dialog window to choose folder with file to be added to the list of the files to read'''
        root = tk.Tk()
        root.withdraw()
        dataDir = filedialog.askdirectory(initialdir = 'E://')
        self.addfilenamesfromfolder(dataDir)
    
    def sortnames(self):
        '''dump files that can't be read(not .mat files and t-cal.mat)'''
        filenamesout=[]
        for item in self.filenames:
            ext = os.path.splitext(item)[-1].lower()
            if ext == ".mat":
                if item != 't-cal.mat':
                    filenamesout.append(item)
                else: print(item + ' was dumped')
            lse: print(item + ' was dumped')
            self.filenames = filenamesout

    def importselectedfiles(self):
        '''import all files from the list of names '''
        for item in self.filenames:
            scan=redred.MOKEScan()
            scan.importRawFile(item)
            scan.filterit()
            scan.removeDC()
            scan.flipTime()
            scan.flipTrace()
            scan.shiftTime(-80)#should add findtimeshift function
            scan.sampleOrient=float(item[item.find('rot-')+4:item.find('degree')])
            scan.initParameters()
            scan.quickplot(clear=True)
            self.scans.append(scan)
#%%functions applying redred functions to every scan in the list
            
    def filterall(self):
        for item in self.scans:
            item=item.filterit()
    
    def removeDCall(self):
        for item in self.scans:
            item=item.removeDC()
            
    def fliptimeall(self):
        for item in self.scans:
            item=item.flipTime()
    
    def fliptraceall(self):
        for item in self.scans:
            item=item.flipTrace
            
    def shiftTimeall(self, tshift):
        for item in self.scans:
            item=item.shiftTime(tshift)
    
    def initParametersall(self):
        for item in self.scans:
            item=item.initParameters()
#%% plot functions
    def rrPlot3d(self, Yparameter='Sample Orientation', title='3dplot', Xlabel= 'Time, ps', Zlabel='Kerr rotation (mrad)', colormap='viridis'):
        '''plot 3d graf with time on X trace on Z and selected parametr on Y '''
        #create 3 lists of X Y Z data
        time=[]
        trace=[]
        ypar=[]
        #for every scan object takes values
        for item in self.scans:
            time.append(item.time)
            trace.append(item.trace)
            #on Y axis will be chosen parameter which exist in scan object
            ypar.append(item.parameters[Yparameter][0])
        #Make proper arrays from lists with data
        Ypar=[]
        
        for item in range(len(self.scans[0].time)):
            
            Ypar.append(ypar)
            
        X=np.array(time)
        Y=np.transpose(np.array(Ypar))
        Z=np.array(trace)
        
        fig = plt.figure(num=2)
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(X, Y, Z, rstride=1, cstride=100, cmap=colormap)
        ax.set_xlabel(Xlabel, fontsize=20, labelpad=20)
        ax.tick_params(axis='x', labelsize=20)
        ax.tick_params(axis='y', labelsize=20)
        ax.tick_params(axis='z', labelsize=20)
    
        ax.set_ylabel(Yparameter, fontsize=20, labelpad=20)
        ax.set_zlabel(Zlabel, fontsize=20, labelpad=20)
        ax.set_title(title, fontsize=40)
        plt.show()
    
sourceDir='E:\\2017-02-24\\'
#Filenames=os.listdir(sourceDir)
#Filenames = sortnames(Filenames)

#Names=[]

#i=0
#while i<37 :
   # Names.append('Mn2Au-MA125-30mWhor-20mWvert-RT-1kAV-MOKE-samplerot-'+str(0+i*10)+'degree.mat')
   # i+=1

Data=RRscans()
#for item in Names:
 #   Data.addfilename(sourceDir + item)
Data.addfilenamesfromfolder(sourceDir)
Data.sortnames()
Data.importselectedfiles()

#Data.rrPlot3d()