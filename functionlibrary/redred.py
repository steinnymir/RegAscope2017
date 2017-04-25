# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 14:11:19 2017

@author: sagustss
"""

#%% Import modules

import os
import numpy as np
import scipy as sp
import scipy.signal as spsig
#GUI
#from pyqtgraph.Qt import QtGui, QtCore
#import pyqtgraph as pg
#plot
import matplotlib.pyplot as plt
#from matplotlib import cm
#other
from datetime import datetime
import re


#%%

def import_file(filename, content = 'Daten'):
    """ 
    Import data aquired with RedRed software
    
    returns data as [time,trace]    
    
    """    
    MData = sp.io.loadmat(filename)    #load matlab file
    output = []
    
    try :
        output = MData[content]
    except KeyError:
        print('KeyError: \nNo key "' + content + '" found in ' + filename)
    
    if content == 'Daten':
        x=[[],[]]
        x[0] = output[2] #assign time axis
        x[1] = output[0] #assgin data axis
        output = x
    
    return(output)


def remove_DC_offset(trace, zerofrom = 7460, zeroto = 7500):
    """remove DC offset from signal"""
    
    shift=np.average(trace[zerofrom:zeroto:1])
    
    newtrace=trace-shift
    
    return(newtrace)


def timezero_shift(timeData, timeZero = 0, reverse = 'False'):
    """set time zero"""
    newtimedata = timeData + max(timeData) - timeZero
    if reverse:
        newtimedata = -newtimedata
    return(newtimedata)


def quick_filter(trace, order = 2, cutfreq = 0.1):
    """ apply low pass filter to data"""
    b, a = spsig.butter(order, cutfreq, 'low', analog= False)
    filtered_trace = spsig.lfilter(b,a,trace)
    return(filtered_trace)


def dir_to_dict(sourceDirectory, fileRange = [0,0]):
    """ Generate a dictionary containing info from file name and data"""
    
    #select all files if range is [0,0]
    if fileRange == [0,0]:
        fileRange[1] = len(sourceDirectory)
        
        # pick scans to work on
    fileNames = os.listdir(sourceDirectory)[fileRange[0]:fileRange[1]]
    
    DataDict = {}
    nGood, nBad = 0,0
    for item in fileNames:
        ext = os.path.splitext(item)[-1].lower()
        if ext == ".mat":
            nGood += 1
            filepath = sourceDirectory + '//' + item                                    
            DataDict[item] = name_to_info(filepath)
            DataDict[item]['data'] = import_file(filepath)
        else:
            nBad += 1
            
    print('Imported '+ str(nGood) + ' Files \n discarded '+ str(nBad) + ' Files')
    
    return(DataDict)

#%% filename interpreter

def name_to_info(file):
    """ 
    Interprets measurement parameters out of the file name
    
    it understands the following:
    'Pump Power' : 'pu','pump',
    'Probe Power': 'pr','probe',
    'Temperature': 't','temp',
    'Destruction Power': 'd','dest',
       
    """   
    #lowercase file name without .* extension
    try:
        FileName = ('.').join(os.path.basename(file).split('.')[:-1]) 
        filename = FileName.lower()
        
        #define parameter set
        parameterDict = {'Scan Date' : file_creation_date(file),
                         'Material'   : '-',
                         'Pump Power' : '-',
                         'Probe Power': '-',
                         'Temperature': '-',
                         'Destruction Power': '-',
                         'Other': '-',
                         }
        #set of possible indicators and corresponding key 
        parInd = {'Pump Power' : ['pu','pump'],
                  'Probe Power': ['pr','probe'],
                  'Temperature': ['t','temp'],
                  'Destruction Power': ['d','dest'],
                  }
        # Iterate over possible indicators and save value to corresponding key in 
        parameterIndicatorTrue = []
        for key, indicators in parInd.items():
            for string in indicators:
                if string in filename:
                    parameterIndicatorTrue.append(string)
                    # Partition string around parameter delimiter and pick the  
                    # first following float number
                    value = float(re.findall(r"\d+\.\d*", 
                                             filename.partition(string)[2])[0]) 
                    #append in Dictionary of parameters
                    parameterDict[key] = value
        pos = 100
        for string in parameterIndicatorTrue:
           x = filename.find(string)
           if x < pos: pos=x
           
        if pos != 0:
            parameterDict['Material'] = FileName[0:pos]
        #print(FileName[pos])
        
        if FileName[pos-1] in ['_','-']:
            parameterDict['Material'] = FileName[0:pos-1]
    except IndexError:
        print(FileName + ' contains no parameter info')
        pass
        
    return(parameterDict)

#%% generic utility functions

def file_creation_date(file):
    return datetime.fromtimestamp(int(os.path.getctime(file))).strftime('%Y-%m-%d %H:%M:%S')

#%% run main
if __name__ == "__main__":
    
    testfile = 'RuCl3-Pr-0.5mW-Pu-1.5mW-T-007.0k-1kAVG.mat'
    testpath = 'c://Users//sagustss//py_code//test_data'
    
    dataDict = dir_to_dict(testpath)
    print(dataDict[testfile]) # test dir_to_dict
    
    x = dataDict[testfile]['data'][0]
    y = dataDict[testfile]['data'][1]
    X = timezero_shift(x, timeZero = 50, reverse = True)
    Y = quick_filter(y, order = 3, cutfreq = 0.01)
    
    fig = plt.figure("test")
    plt.clf()
    
    ax1 = fig.add_subplot(111)
    ax1.plot(X,Y)
    
    
    
    