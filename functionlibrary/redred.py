# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 14:11:19 2017

@author: sagustss
"""

#%% Import modules

import numpy as np
import scipy as sp
import scipy.signal as spsignal
import pickle
import csv
import os


import matplotlib.pyplot as plt
#other
from datetime import datetime
import re
from matplotlib import cm
import tkinter as tk
from tkinter import filedialog

#%%
def main():

    """ use a test file to test most funnctions in this file
    now set for norm_to_pump
    """

    testfile = 'RuCl3-Pr-0.5mW-Pu-1.5mW-T-007.0k-1kAVG.mat'
    testpath = '..//test_data//'
    savepath = "E://DATA//RuCl3//Analysis//"


    scn = rrScan()
    scn.importRawFile(testpath + testfile)

    a = scn.fetchMetadata()
    print(a)

#%% class for Redred/MOKE data

class rrScan(object):
    """
    Class defining a scan from redred setup.

    Attributes:
    -----------

        self.time = []          # time data
        self.rawtrace = []      # raw trace data
        self.trace = []         # modified trace data

        self.pumpPw = 0         # Pump Power [mW]
        self.probePw = 0        # Probe Power [mW]
        self.destrPw = 0        # Destruction Power [mW]
        self.pumpSp = 0         # Pump beam spot size, FWHM [micrometers]
        self.probeSp = 0        # Probe beam spot size, FWHM [micrometers]
        self.temperature = 0    # Temperature [K]
        self.date = ''          # Scan date in format YYYY-MM-DD hh.mm.ss
        self.material = ''      # Material name
        self.R0 = 0             # Static reflectivity
        self.filter = []        # Filter parameters used for self.trace

        self.sampleOrient = 0
        self.pumpPol = 0
        self.probePol = 0
        self.originalfilename=''
        self.analysisHistory = [] #keeps track of analysis changes performed
        self.scanID = []        # list of parameters used in the name
        self.filename = []      # String used as file name for saving data
        self.parameters = {}    # Dictionary of all parameters, somewhat redundant?

    """


    def __init__(self):
        """
        Initialize the scan by defining time and differential reflectivity data

        Also define any other parameter as listed
        """

        self.time = []          # time data
        self.rawtrace = []      # raw trace data
        self.trace = []         # modified trace data

        self.pumpPw = 0         # Pump Power [mW]
        self.probePw = 0        # Probe Power [mW]
        self.destrPw = 0        # Destruction Power [mW]
        self.pumpSp = 0         # Pump beam spot size, FWHM [micrometers]
        self.probeSp = 0        # Probe beam spot size, FWHM [micrometers]
        self.temperature = 0    # Temperature [K]
        self.date = ''          # Scan date in format YYYY-MM-DD hh.mm.ss
        self.material = ''      # Material name
        self.R0 = 0             # Static reflectivity
        self.filter = []        # Filter parameters used for self.trace

        self.sampleOrient = 0
        self.pumpPol = 0
        self.probePol = 0
        self.originalfilename=''


        # analysis trackers
        self.analysisHistory = [] #keeps track of analysis changes performed
        self.scanID = []        # list of parameters used in the name
        self.filename = []      # String used as file name for saving data
        self.parameters = {}    # Dictionary of all parameters, somewhat redundant?

        # Dictionary containing info about each metadata:
        # units and relative class attribute name
        self.metadataInfo = {'time' : ['Time Axis', 'ps'],
                             'rawtrace' : ['Y axis', '?'],
                             'trace' : ['Y axis', '?'],
                             'pumpPw' : ['Pump Power', 'mW'],
                             'probePw' : ['Probe Power', 'mW'],
                             'destrPw' : ['Destruction Power', 'mW'],
                             'pumpSp' : ['Pump Spot Size', 'microm'],
                             'probeSp' : ['Probe Spot Size', sym('mu') + 'm'],
                             'temperature' : ['Temperature', 'K'],
                             'date' : ['Date', ''],
                             'material' : ['Material', ''],
                             'R0' : ['R0', ''],
                             'filter' : ['Filter Frequency', 'NyqFreq Multiplier'],
                             'sampleOrient' : ['Sample Orientation', 'rad'],
                             'pumpPol' : ['Pump Polarization', 'rad'],
                             'probePol' : ['Probe Polarization', 'rad'],
                             'originalfilename' : ['Raw File Path', 'str'],
                           }


    def initParameters(self):
        """ Create a a dictionary of all parameters and a nameID for the scan"""
        self.parameters = {'Pump Power': [self.pumpPw,'mW'],
                           'Probe Power': [self.probePw,'mW'],
                           'Destruction Power': [self.destrPw,'mW'],
                           'Pump Spot': [self.pumpSp,'mum'],
                           'Pump polarization': [self.pumpPol,'deg'],
                           'Probe polarization': [self.probePol,'deg'],
                           'Sample Orientation': [self.sampleOrient,'deg'],
                           'Probe Spot': [self.probeSp,'mum'],
                           'Temperature': [self.temperature,'K'],
                           'R0': [self.R0,'']
                           }
        # generate scanID
        if self.material:
            self.scanID.append(self.material)
        else:
            self.scanID.append('UnknMat')
        if self.date:
            self.scanID.append(self.date)
        if self.pumpPw != 0:
            self.scanID.append('Pump' + str(self.pumpPw) + 'mW')
        if self.destrPw != 0:
            self.scanID.append('Destr' + str(self.destrPw) + 'mW')
        if self.temperature != 0:
            self.scanID.append('Temp' + str(self.temperature) + 'K')
        self.filename = ' '.join(self.scanID)
        self.filename = self.filename.replace(':','.')


    def fetchMetadata(self):
        """ Returns a dictionary containing all metadata information available"""

        metadata = {}
        var = self.__dict__
        for key in var:
            if key in self.metadataInfo:
                if not key in ['time', 'trace', 'rawtrace']:
                    try:
                        if float(var[key]) != 0:
                            metadata[key] = var[key]
                    except ValueError:
                        metadata[key] = var[key]
                    except TypeError:
                        metadata[key] = var[key]

        return(metadata)

    def pushMetadata(self, dict):
        """ Import metadata from dictionary
        Dictionary keys must be full names, ex: "Pump Power", not PumpPw
        """

        self.pumpPw = dict['PumpPower']
        self.probePw = dict['Probe Power']        # Probe Power [mW]
        self.destrPw = dict['Destruction Power']        # Destruction Power [mW]
        self.pumpSp = dict['Pump Spot Size']         # Pump beam spot size, FWHM [micrometers]
        self.probeSp = dict['Probe Spot Size']        # Probe beam spot size, FWHM [micrometers]
        self.temperature = dict['Temperature']    # Temperature [K]
        self.date = dict['Date']          # Scan date in format YYYY-MM-DD hh.mm.ss
        self.material = dict['Material']      # Material name
        self.R0 = dict['R0']             # Static reflectivity
        self.filter = dict['Filter Frequency']        # Filter parameters used for self.trace

        self.sampleOrient = dict['Sample Orientation']
        self.pumpPol = dict['Pump Polarization']
        self.probePol = dict['Probe Polarization']
        self.originalfilename = dict['Raw File Path']

        self.analysisHistory = dict['Analysis History'] #keeps track of analysis changes performed



    def shiftTime(self, tshift):
        """ Shift time scale by tshift. Changes time zero"""
        self.time = np.array(self.time) - tshift
        self.analysisHistory.append('shift time')


    def flipTime(self):
        """ Flip time scale: t = -t
        does not revert the list order"""
        self.time = -self.time
        #self.time = self.time[::-1]
        #self.trace = self.trace[::-1]
        self.analysisHistory.append('flip time')

    def flipTrace(self):
        """ Flip the trace, usually not needed"""
        self.trace = -self.trace
        self.analysisHistory.append('flip trace')

    def removeDC(self):
        """Remove DC offset defined by the average of the last(which are the first) 40 points on the scan"""
        self.trace = self.trace - np.average(self.trace[7460:7500:1])
        self.analysisHistory.append('removeDC')

    def filterit(self, cutHigh = 0.1, order = 2):
        """ apply simple low pass filter to data"""
        b, a = spsignal.butter(order, cutHigh, 'low', analog= False)
        self.trace = spsignal.lfilter(b,a,self.rawtrace)
        self.filter = [cutHigh, order]
        self.analysisHistory.append('filter')

    def filterFreq(self):
        """ Gives low pass filter frequency in THz """
        nyqFreq = abs(0.5 * len(self.time) / self.time[-1] - self.time[0])
        return(nyqFreq * self.filter[0])

    def normToPump(self):
        if self.pumpPw != 0:
            self.trace = self.trace / self.pumpPw
        self.analysisHistory.append('normalized to PumpPw')

    def quickplot(self, xlabel='Time, ps',
                  ylabel='Kerr rotation', fntsize=20,
                  title='Time depandance of the pump induced Kerr rotation',
                  clear=False):
        if clear: plt.clf()
        quickplotfig=plt.figure(num=1)
        ax=quickplotfig.add_subplot(111)
        ax.plot(self.time, self.trace,)
        ax.set_xlabel(xlabel, fontsize=fntsize)
        ax.set_ylabel(ylabel, fontsize=fntsize)
        ax.set_title(title,fontsize=fntsize)
        ax.tick_params(axis='x', labelsize=fntsize)
        ax.tick_params(axis='y', labelsize=fntsize)
        plt.show()


#%%file managment

    def importRawFile(self, file):
        data = sp.io.loadmat(file)
        try:
            self.time = data['Daten'][2]
            self.rawtrace = data['Daten'][0]
            self.trace = self.rawtrace
            self.R0 = data['DC'][0][0]

        except KeyError:
            print(file + ' is not a valid redred scan datafile')
        parDict = name_to_info(file)

        for key in parDict:
            if key == 'Scan Date':
                self.date = parDict[key]
            elif key == 'Pump Power':
                self.pumpPw = parDict[key]
            elif key == 'Probe Power':
                self.probePw = parDict[key]
            elif key == 'Temperature':
                self.temperature = parDict[key]
            elif key == 'Destruction Power':
                self.destrPw = parDict[key]
            elif key == 'Material':
                self.material = parDict[key]
                #print(self.material)
            elif key == 'Pump Spot':
                self.pumpSp = parDict[key]
            elif key == 'Probe Spot':
                self.probeSp = parDict[key]
            elif key == 'Other':
                self.other = parDict[key]
            else:
                print('Unidentified Key: ' + key)

        self.originalfilename=file
        self.initParameters()


    def importCSV(self,file):
        """
        Read a CSV containing rrScan() data, and assign to self all the data
        file should be a string of the whole path of the file"""
        try:
            f = open(file, 'r')
        except FileNotFoundError:
            print('ERROR 404: file not found')
        if f:
            metacounter=0
            for l in f:
                metacounter+=1
                line = l.split('\t')
                if 'Material' in line: self.material = line[1]
                elif 'Pump Power' in line: self.pumpPw = float(line[1])
                elif 'Pump Spot' in line: self.pumpSp = float(line[1])
                elif 'Probe Power' in line: self.probePw = float(line[1])
                elif 'Probe Spot' in line: self.probeSp = float(line[1])
                elif 'date' in line: self.date = line[1]
                elif 'Destruction Power' in line: self.destrPw = float(line[1])
                elif 'R0' in line: self.R0 = float(line[1])
                elif 'Temperature' in line: self.temperature = float(line[1])
                elif 'RawData' in line :
                    break
                    print(metacounter)
            f.close()
            skipline = True
            n=0
            data = []

            while skipline: # skip metadata section then import array of data
                try:
                    data = np.loadtxt(file, delimiter=',', skiprows=n)
                    n+=1
                    skipline = False
                except ValueError:
                    n+=1

            for i in range(len(data)):
                self.time.append(data[i][1])
                self.rawtrace.append(data[i][1])
                self.trace.append(data[i][1])


    def exportCSV_old(self, directory):
        """ save rrScan() to a file. it overwrites anything it finds"""

        file = open(directory + self.filename + '.txt', 'w+')

        # Material:
        file.write('Material:\t' + str(self.material) + '\n')
        # Date:
        file.write('Date:\t' + self.date + '\n')
        # Parameters
        file.write('------- Parameters -------\n\n')
        for key in self.parameters:
            if self.parameters[key][0] != 0:
                file.write(key + '\t' +
                           str(self.parameters[key][0]) + '\t' +
                           str(self.parameters[key][1]) + '\n' )
        #filter info
        if self.filter:
            file.write('\nLow pass filter frequency:\t' +
                       str(self.filter) +'\t' +
                       str(self.filterFreq()) +
                       'THz\n')
        #data
        file.write('---------- Data ----------\n\n')
        file.write('Time\tRawData\tdata\n')
        for i in range(len(self.trace)):
            file.write(str(self.time[i])     + ',' +
                       str(self.rawtrace[i]) + ',' +
                       str(self.trace[i])    + '\n')
            self.time = np.array(self.time)
        file.close()

    def exportCSV(self, directory):
        """
        save rrScan() to a file. it overwrites anything it finds

        Metadata is obtained from fetchMetadata(), resulting in all non0
        parameters available.
        """
        file = open(directory + self.filename + '.txt', 'w+')
        metadata = self.fetchMetadata()
#        for key in metadata:
#            name =  str(self.metadataInfo[key][0])
#            value = str(metadata[key])
#            unit = str(self.metadataInfo[key][1])
#            file.write(name + '/t' + value + '/t' + unit = '/n')
        #data
        file.write('++++++++++ Data ++++++++++\n\n')
        file.write('Time\tRawData\tdata\n')
        for i in range(len(self.time)):
            file.write(str(self.time[i])     + ',' +
                       str(self.rawtrace[i]) + ',' +
                       str(self.trace[i])    + '\n')
            self.time = np.array(self.time)
        file.close()



 #%%  class for multiple scans
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

    def importselectedfiles(self, plot=False):
        '''import all files from the list of names '''
        for item in self.filenames:
            scan=rrScan()
            scan.importRawFile(item)
            if plot: scan.quickplot()
            self.scans.append(scan)

    def definesampleorientforoldmokedata(self):
        for item in self.scans:
            item.sampleOrient=float(item.originalfilename[item.originalfilename.find('rot-')+4:item.originalfilename.find('degree')])
#%%functions applying redred functions to every scan in the list

    def filteritall(self, cutHigh = 0.1, order = 2):
        for item in self.scans:
            item=item.filterit(cutHigh, order)

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
#%%
#%% Functions
def import_file(filename, content = 'Daten'):
    """Import data aquired with RedRed software
    returns data as [time,trace]"""
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
    b, a = sp.signal.butter(order, cutfreq, 'low', analog= False)
    filtered_trace = sp.signal.lfilter(b,a,trace)
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
        #rest = norm-dataDict[key]['data'][1][1]
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
    'Destruction Power': 'd','dest','destr'

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
              'Destruction Power': ['d','dest','destr'],
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
                    parameterDict[key] = value
                except IndexError:
                    pass
                #append in Dictionary of parameters


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

def sym(letter):
    greek_alphabet = {
    u'\u0391': 'Alpha',
    u'\u0392': 'Beta',
    u'\u0393': 'Gamma',
    u'\u0394': 'Delta',
    u'\u0395': 'Epsilon',
    u'\u0396': 'Zeta',
    u'\u0397': 'Eta',
    u'\u0398': 'Theta',
    u'\u0399': 'Iota',
    u'\u039A': 'Kappa',
    u'\u039B': 'Lamda',
    u'\u039C': 'Mu',
    u'\u039D': 'Nu',
    u'\u039E': 'Xi',
    u'\u039F': 'Omicron',
    u'\u03A0': 'Pi',
    u'\u03A1': 'Rho',
    u'\u03A3': 'Sigma',
    u'\u03A4': 'Tau',
    u'\u03A5': 'Upsilon',
    u'\u03A6': 'Phi',
    u'\u03A7': 'Chi',
    u'\u03A8': 'Psi',
    u'\u03A9': 'Omega',
    u'\u03B1': 'alpha',
    u'\u03B2': 'beta',
    u'\u03B3': 'gamma',
    u'\u03B4': 'delta',
    u'\u03B5': 'epsilon',
    u'\u03B6': 'zeta',
    u'\u03B7': 'eta',
    u'\u03B8': 'theta',
    u'\u03B9': 'iota',
    u'\u03BA': 'kappa',
    u'\u03BB': 'lamda',
    u'\u03BC': 'mu',
    u'\u03BD': 'nu',
    u'\u03BE': 'xi',
    u'\u03BF': 'omicron',
    u'\u03C0': 'pi',
    u'\u03C1': 'rho',
    u'\u03C3': 'sigma',
    u'\u03C4': 'tau',
    u'\u03C5': 'upsilon',
    u'\u03C6': 'phi',
    u'\u03C7': 'chi',
    u'\u03C8': 'psi',
    u'\u03C9': 'omega',
}
    inv_map = {v: k for k, v in greek_alphabet.items()}
    for key in inv_map:
        if letter == key:
            return(inv_map[key])
#%% run main
if __name__ == "__main__":
    main()





