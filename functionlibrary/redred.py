# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 14:11:19 2017

@author: sagustss
"""

#%% Import modules
import pickle
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
from matplotlib import cm

#%%
def main():
    
    """
    use a test file to test most funnctions in this file
    
    now set for norm_to_pump
    
    """
    
    testfile = 'RuCl3-Pr-0.5mW-Pu-1.5mW-T-007.0k-1kAVG.mat'
    testpath = '..//test_data'
    
    dataDict = dir_to_dict(testpath)
    print(dataDict[testfile]['data'][1]) # test dir_to_dict

    dataDictNorm = norm_to_pump(dir_to_dict(testpath))
    print(dataDict[testfile]['data'][1]) # test dir_to_dict
#    diff = dataDictNorm[testfile]['data'][1][1] - dataDict[testfile]['data'][1][1]
#    print(diff)
    
    fig = plt.figure("Raw Trace")
    plt.clf()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    
    for key in dataDict:
        #print(str(dataDict[key]['Pump Power']))
        x = dataDict[key]['data'][0]
        X = timezero_shift(x, timeZero = 50, reverse = True)
        
        y = dataDict[key]['data'][1]
        Y = quick_filter(y, order = 3, cutfreq = 0.01)
        Yf = remove_DC_offset(Y)
        
        yn = dataDictNorm[key]['data'][1]
        Yn = quick_filter(yn, order = 3, cutfreq = 0.01)
        Ynf = remove_DC_offset(Yn)
        
        ax1.plot(X,Yf)
        ax2.plot(X,Ynf)
#%%
class rrScan(object):
    """ class defining a scan from redred setup
        
        Attributes
        ----------
        
        time        : time axis [ps]
        trace       : differential reflectivity axis 
        pump        : pump power [mW]
        probe       : probe power [mW]
        puspot      : pump spot diameter (FWHM gauss) [micrometer]
        prspot      : probe spot diameter (FWHM gauss) [micrometer]
        temperature : temperature [K]
        date        : scan date and time
        material    : material name
    
    """
    
    
    def __init__(self, time, trace):
        """ 
        Initialize the scan by defining time and differential reflectivity data
        
        Also define any other parameter as listed
        """
        self.time = time
        self.trace = trace
        self.pump = 0
        self.probe = 0
        self.puspot = 0
        self.prspot = 0
        self.temperature = 0
        self.date = ''
        self.material = ''

    def shiftTime(self, tshift):
        """ Shift time scale by tshift. Changes time zero"""
        shiftedTime = self.time - tshift
        return(shiftedTime)
    
    def flipTime(self):
        """ Flip time scale: t = -t """
        revtime = -self.time
        return(revtime)
    
    def removeDC(self):
        """Remove DC offset defined using last 50 points on the scan"""
        noDCtrace = self.trace - np.average(self.trace[7460:7500:1])
        return(noDCtrace)

    def filterit(self, cutHigh = 0.1, cutLow = 0, order = 2):
        pass

#%%
def import_file(filename, content = 'Daten'):
    """ 
    Import data aquired with RedRed software
    
    returns data as [time,trace]
    
    """    
    MData = sp.io.loadmat(filename)    #load matlab file
    output = []
    
    if filename == 't-cal.mat':
        pass
    else:
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
    """Shift the 0 offset of a time trace, returns [new time trace] and [time shift]
    
    Some time its better to define time shift from one trace and aply it to athers, 
    since max or min value could be different, like in my case
    
    -> You can define the time shift with a single line (timeshift = max(TimeData))
    and feed it to this function as "timeZero", so no need to have two outputs 
    from this function.
    
    
    I add a comment
    """
    
    #timeshift = max(timeData)
    #newtimedata = timeData + timeshift - timeZero
    timeData = timeData - timeZero
    if reverse:
        timeData = -timeData
        
    #return(newtimedata, timeshift) 
    return(timeData)


def quick_filter(trace, order = 2, cutfreq = 0.1):
    """ apply simple low pass filter to data"""
    b, a = spsig.butter(order, cutfreq, 'low', analog= False)
    filtered_trace = spsig.lfilter(b,a,trace)
    return(filtered_trace)





def file_to_dict(filepath):
    """
    Convert file into Dictionary containing scan info and data
    
    if file is not valid returns an empty dictionary
    """
    DataDict = {}    
    filename  = os.path.basename(filepath)
    #print(filename)    
    ext = os.path.splitext(filename)[-1].lower()
     
    if ext == ".mat" and not filename == 't-cal':
        DataDict = name_to_info(filepath)
        DataDict['data'] = import_file(filepath)
        return(DataDict)




def dir_to_dict(sourceDirectory, fileRange = [0,0]):
    """ Generate a dictionary containing info from file name and data"""
    
    #select all files if range is [0,0]
    if fileRange == [0,0] or fileRange[1]<fileRange[0]:
        fileRange[1] = len(sourceDirectory)
        
        # pick scans to work on
    fileNames = os.listdir(sourceDirectory)[fileRange[0]:fileRange[1]]
    
    DataDict = {}
    nGood, nBad = 0,0
    for item in fileNames:
        filepath = sourceDirectory + '//' + item
        filename = os.path.basename(filepath)
        ext = os.path.splitext(filename)[-1].lower()
        #if not os.path.isdir(filepath) and
        if ext == ".mat" and not filename == 't-cal':
            DataDict[item] = file_to_dict(filepath) #returns nothing if file was not datafile.mat
            if DataDict[item]: 
                nGood += 1
                DataDict[item]['data'] = import_file(filepath)
            else:
                nBad += 1
            
    print('Imported '+ str(nGood) + ' Files \n discarded '+ str(nBad) + ' Files')
    
    return(DataDict)





def save_trace(dataDict, directory, filename):
    """ Generate csv file with header"""
    
    directory += '//'
    file = open(directory + filename, "w+")
    
    parameters = ['Material',
            'Scan Date',
            'Other',
            'Probe Power',
            'Pump Power',
            'Temperature',]
    units = ['','','','mW','mW','K']
    # create Header
    u = 0
    for par in parameters:
        
        string = str(dataDict[par])
        file.write(par + '\t' + string + '\t' +units[u]+ '\n')
        u += 1
    file.write('time\ttrace\n')
    for i in range(len(dataDict['data'][1])):
        time = str(dataDict['data'][0][i])
        trace = str(dataDict['data'][1][i])
        file.write(time + ',' + trace + '\n')
    print('done')
    file.close()
        
        
def save_filtered_trace(dataDict, directory, filename, filtermultiplier):
    """ Generate csv file with header"""
    
    directory += '//'
    file = open(directory + filename, "w+")
    
    parameters = ['Material',
            'Scan Date',
            'Other',
            'Probe Power',
            'Pump Power',
            'Temperature',]
    units = ['','','','mW','mW','K']
    # create Header
    u = 0
    for par in parameters:
        
        string = str(dataDict[par])
        file.write(par + '\t' + string + '\t' +units[u]+ '\n')
        u += 1
    
    
    time = dataDict['data'][0]
    # get filter cut of frequency
    nyqFreq = 0.5 * len(time) / time[-1]-time[0]
    filterFreq = nyqFreq*filtermultiplier
    
    raw = dataDict['data'][1]
    #filter data
    filt = quick_filter(dataDict['data'][1], cutfreq = filtermultiplier)
    
    
    file.write('Filter frequency\t'+str(filterFreq)[0:4]+'THz\n')
    file.write('time\traw trace\tfilteredtrace\n')
    
  
    
    for i in range(len(dataDict['data'][1])):
        stime = str(time[i])
        sraw = str(raw[i])
        sfilt = str(filt[i])        
        file.write(stime + ',' + sraw + ',' + sfilt + '\n')
    print('file ' + filename + ' is ready')
    file.close()
    


def quickPlot(sourceDirectory, cutfreq, KeyDependence = 'Temperature'):
        
    dataDict = dir_to_dict(sourceDirectory)
    # initialize regualr plot
    
    Title = KeyDependence + ' Dependence'

    fig = plt.figure(Title, figsize = (19,10))
    plt.clf()
    plt1 = fig.add_subplot(111)
    plt1.set_xlabel('Time, ps', fontsize=18)
    plt1.set_ylabel('Differential Reflectivity', fontsize=18)
    plt1.set_title(Title,fontsize=26)
    plt1.tick_params(axis='x', labelsize=12)
    plt1.tick_params(axis='y', labelsize=12)
    color=iter(cm.rainbow(np.linspace(0,1,len(dataDict))))
      
    # Plot a "key" dependence from 
    for key in dataDict:
        x = timezero_shift(dataDict[key]['data'][0], reverse = True)
        y = quick_filter(remove_DC_offset(dataDict[key]['data'][1]), 
                            order = 2, cutfreq = cutfreq)

        col = next(color)
        plt1.plot(x,y, label = dataDict[key][KeyDependence], c=col)
    
        
    plt.show()
    plt.legend(bbox_to_anchor=(1.005, 1), loc='best', borderaxespad=0., ncol=1)
    
    
def norm_to_pump(dataDict):
    """ Divide all curves in dataDict by it's pump power value"""
    dataDictNorm = dataDict
    norm = []
    for key in dataDict:
        norm = dataDict[key]['data'][1] / dataDict[key]['Pump Power']#dataDict[key]['Pump Power']
        rest = norm-dataDict[key]['data'][1][1]
        dataDictNorm[key]['data'][1] = norm
        #print('rest:  '+str(rest))
        
    return(dataDictNorm)

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
    FileName = ('.').join(os.path.basename(file).split('.')[:-1]) 
    filenamex = FileName.lower()
    filename = filenamex.replace(',','.')
    #print(filename)
    
    #errStr = FileName + ' contains no info about: '
    
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
                try:
                    value = float(re.findall(r"\d+\.\d*", 
                                         filename.partition(string)[2])[0])
                except IndexError:
                    pass
                #append in Dictionary of parameters
                
                parameterDict[key] = value
    pos = 100
    for string in parameterIndicatorTrue:
       x = filename.find(string)
       if x < pos: 
           pos=x
    
    parameterDict['Material'] = FileName[0:pos]

#    if parameterDict['Material'][-1] in ['_','-']:
#        parameterDict['Material'] = parameterDict['Material'][:-1]

    

#    except IndexError:
#        print(FileName + ' contains no parameter info')
#        pass
        
    return(parameterDict)


#%% generic utility functions

def file_creation_date(file):
    return datetime.fromtimestamp(int(os.path.getmtime(file))).strftime('%Y-%m-%d %H:%M:%S')

def save_obj(obj, name ):
    with open('obj/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open('obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)

#%% run main
if __name__ == "__main__":
    main()    
    
    
    
    
    
