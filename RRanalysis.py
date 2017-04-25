# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 13:19:22 2017

RedRed Analysis, a first try!

@author: S.Y.Agustsson, V.Yu. Grigorev
"""
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import re
from datetime import datetime

#%%  function definitions
def rr_import(Filename, fltrOrder = 2, fltrCut = 0.1, timeZero = 0):
    
    """import redred data and apply lowpass filter"""
    
    MData = sp.io.loadmat(Filename)    #load matlab file
    Data = MData['Daten']           # Pick data with dark control
    d=Data[0]                       
    shift=np.average(d[7460:7500:1])    #shift zero to pre-pump value
    trace=(Data[0]-shift)

    timeShift = max(Data[2]) - timeZero
    Time = -Data[2] + timeShift #+ Data[2,np.argmin(-Data[0])]

    b, a = sp.signal.butter(fltrOrder, fltrCut, 'low', analog= False)
    filtered_trace = sp.signal.lfilter(b,a,trace)
   # plt.plot(Time,filtered_trace)
    return(Time,filtered_trace)


def rr_dictionary(sourceDirectory, fileRange = [0,0], save = True, 
                  fltrOrder = 2, fltrCut = 0.1, timeZero = 0):
    """Generate a Dictionary and save transformed data"""
    
    if fileRange == [0,0]:
        fileRange[1] = len(sourceDirectory)
    saveName = os.path.basename(sourceDirectory)
        # pick scans to work on
    fileNames=os.listdir(sourceDirectory)[fileRange[0]:fileRange[1]] 
    data = []
    countGood = 0
    countBad = 0
    for item in fileNames:
        ext = os.path.splitext(item)[-1].lower()
        if ext == ".mat":
            countGood += 1
            D = rr_import(sourceDirectory + '//' + item, 
                          fltrOrder, fltrCut, timeZero) # import scan data
            d = np.array(D) #data to array
             # add label and data for dictionary - Matlab export
            data.append(['T_' + item[0:-4], d])
        else: countBad += 1
    dataDict = dict(data) # Transform to dictionary
    print(str(countGood) + " files used\n" + str(countBad) + " files were ignored.")
    if save:  # save as matlab file in targetDirectory
        sp.io.savemat(saveName + '_FilteredData', mdict=dataDict)
        print('Saved as ' + saveName + '_FilteredData')
    else: print('not Saved')
    return(dataDict)


def rr_plot(dataDict, saveName):
    """Plot all curves in the dataDict dictionary and use keys as labels"""
    fig = plt.figure(saveName, figsize = (10,20))
    plt.clf()
    ax1 = fig.add_subplot(111)
    
    color=iter(cm.rainbow(np.linspace(0,1,len(dataDict))))
    
    for key, value in dataDict.items():
            
        time = value[0]
        trace = value[1]
        col=next(color)
        
        ax1.plot(time, trace, label=key[26:-2], c=col)
    
    
    ax1.set_xlabel('Time, ps', fontsize=20)
    ax1.set_ylabel('Differential Reflectivity', fontsize=20)
    ax1.set_title(saveName,fontsize=30)
    ax1.tick_params(axis='x', labelsize=15)
    ax1.tick_params(axis='y', labelsize=15)
    plt.show()
    plt.legend(bbox_to_anchor=(1.005, 1.1), loc=0, borderaxespad=-1., ncol=1)
    print('-  Check new Plot!!  -')


def name_interpreter_delimiter(file):
    
    #change: use most common simbol as delimiter
    filename = os.path.basename(file)
    print(filename)
    if filename.rfind('-'): 
        delim = '-'
    elif filename.rfind('_'):
        delim = '_'
    else:
        print('-- No delimiter Found --')
        
    filenameList = filename[0:-4].split(delim)
    print(filenameList)
    #parameterDict = {nameList[i]: nameList[i+1] for i in range(0, len(nameList), 2)}
    #append name
    if filenameList[0].lower() in ['pu','pr','pump','probe','temp','t']:
        parameterDict = {'Material': 'Unknown'}
    else:
        parameterDict = {'Material': filenameList[0]}
        
    #add measure date and time
    parameterDict['Measure date'] = datetime.fromtimestamp(int(os.path.getctime(file))).strftime('%Y-%m-%d %H:%M:%S')
    
    for i in range(len(filenameList)):
        if 'Pu'or'Pump' in filenameList:
            parameterDict['Pump Power'] = filenameList[filenameList.index('Pu')+1]
        if 'Pr'or'Probe' in filenameList:
            parameterDict['Probe Power'] = filenameList[filenameList.index('Pr')+1]
        if 'T'or 'Temp' or 'temp' in filenameList:
            parameterDict['Temperature'] = filenameList[filenameList.index('T')+1]

    return(parameterDict)




def name_to_info(file):
    """get measurement parameters from file name"""
    #lowercase file name without .mat
    FileName = os.path.basename(file)[0:-4]
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
    #iterate over possible indicators and save value to corresponding key in 
    parIndTrue = []
    for key, indicators in parInd.items():
        for string in indicators:
            if string in filename:
                parIndTrue.append(string)
                #partition string around parameter delimiter and pick the first following float number
                value = float(re.findall(r"\d+\.\d*", filename.partition(string)[2])[0]) 
                #append in Dictionary of parameters
                parameterDict[key] = value
    pos = 100
    for string in parIndTrue:
       x = filename.find(string)
       if x < pos: pos=x
    if pos != 0:
        parameterDict['Material'] = FileName[0:pos]
    print(FileName[pos])
    
    if FileName[pos-1] in ['_','-']:
        
        parameterDict['Material'] = FileName[0:pos-1]
    
        
        
    return(parameterDict)

def file_creation_date(file):
    return datetime.fromtimestamp(int(os.path.getctime(file))).strftime('%Y-%m-%d %H:%M:%S')




def fitFunction(t,taupulse,taudecay,linearconst,A,C,t_zero):
    return  (A*(-np.exp(taupulse**2/(8*taudecay**2)) * 
            (np.exp(-t/taudecay))+linearconst*t + C) * 
            (1+sp.special.erf(((2**0.5) * 
            (t+t_zero-(taupulse**2/(4*taudecay))))/taupulse)))


def fitRedRed(trace1):
    trace=trace1[1]
    Time=trace1[0]
    guess=( 0.01,1.68,0.1,0.4,-0.27,0.2)
    #       t,taupulse,taudecay,linearconst,A,C,t_zero
    plt.plot(Time,trace)
    popt, pcov = sp.optimize.curve_fit(fitFunction,Time,trace,p0=guess)
    popt=1
    #fig=plt.figure('fittings')
    #plt.clf()
    #plt.plot(Time,fitFunction(Time, *popt))
    #plt.plot(trace[1],fitFunction(trace[1], 0.16,1.63,0.0011,0.4,-0.27,0.2))
    #plt.plot(Time,trace - fitFunction(Time,  *popt))
    #plt.show()
    return(popt)


def fftRedRed(trace, Time):
    freq=sp.fftpack.fftfreq(trace.size,np.abs((Time[0]-Time[:1])/Time.size))
    FFT = abs(sp.fft(trace))
    plt.plot(freq, FFT)
    
#%% 

#Plot = True
#Save = False
#
sourceDir = 'E://RuCl3//2017-04-11//TemperatureDependence'
sourceDir = 'E://RuCl3//2017-04-11//no_r0'
sourceDir = 'E://RuCl3//2017-04-19//TempDependence04-19'
#saveName = os.path.basename(sourceDir) 
#n_of_scans = len(os.listdir(sourceDir))
#
#dataDict = rr_dictionary(sourceDir, save=Save, fltrOrder = 2, fltrCut = 0.05, timeZero=4)
#
#rr_plot(dataDict, saveName)

#for i in range(5):
#    scanTrace = dataDict[list(dataDict.keys())[i]]
#    fitRedRed(scanTrace)
files = os.listdir(sourceDir)
file=sourceDir + '//' + files[15]
print(name_to_info(file))


#%% Plot Settings





