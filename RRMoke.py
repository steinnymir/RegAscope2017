# -*- coding: utf-8 -*-
"""
Created on Thu May  4 16:35:21 2017

@author: vgrigore
"""

from functionlibrary import redred
import scipy
import os



def sortnames(filenames):
    filenamesout=[]
    for item in filenames:
        ext = os.path.splitext(item)[-1].lower()
        if ext == ".mat":
            if item != 't-cal.mat':
                filenamesout.append(item)
            else: print(item + ' was dumped')
        else: print(item + ' was dumped')
    return(filenamesout)

def importrrm(filenames, sourceDir):
    Data=[]
    for item in filenames:
        scan=redred.MOKEScan()
        scan.importRawFile(sourceDir+'\\'+item)
        scan.initParameters()
        scan.filterit()
        scan.removeDC()
        scan.flipTime()
        scan.shiftTime(-80)#should add findtimeshift function
        scan.quickplot()
        Data.append(scan)
    return(Data)




sourceDir='E:\\2017-02-24'
Filenames=os.listdir(sourceDir)
Filenames = sortnames(Filenames)

Data=importrrm(Filenames, sourceDir)
#Data[9].quickplot()