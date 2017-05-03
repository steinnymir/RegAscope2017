# -*- coding: utf-8 -*-
"""
Created on Wed May  3 10:37:19 2017

@author: sagustss
"""

from functionlibrary import redred as rr
import os
import numpy as np
import scipy
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog
from matplotlib import cm

file = 'test_data//RuCl3-Pr-0.5mW-Pu-1.5mW-T-005.0k-1kAVG.mat'

root = tk.Tk()
root.withdraw()
dataDir = filedialog.askdirectory(initialdir = 'E://DATA//_RAW')


scanlist = os.listdir(dataDir)

Title = 'Low Fluence Temperature Dependence 2'
Title = os.path.basename(dataDir)


fig = plt.figure(Title, figsize = (19,10))
plt.clf()
plt1 = fig.add_subplot(121)
plt2 = fig.add_subplot(122)
plt1.set_xlabel('Time [ps]', fontsize=18)
plt1.set_ylabel('Differential Reflectivity', fontsize=18)
plt1.set_title('Fitted Scans',fontsize=26)
plt1.tick_params(axis='x', labelsize=12)
plt1.tick_params(axis='y', labelsize=12)
plt2.set_xlabel('Temperature [K]', fontsize=18)
plt2.set_ylabel('Decay Time [ps]', fontsize=18)
plt2.set_title('Decay time vs Temperature',fontsize=26)
plt2.tick_params(axis='x', labelsize=12)
plt2.tick_params(axis='y', labelsize=12)
colorlist = cm.rainbow(np.linspace(0,1,len(os.listdir(dataDir))))
color=iter(colorlist)




#fig2 = plt.figure('decay times')
#plt2 = fig2.add_subplot(111)


scn = []
#define fitting function
def func(x, A, t0, c, d):
    return A * np.exp(- x / t0) + c*x +d
guess = [0.0005, 0.1, -1, -1]
fitparameters = []


#scan through folder and create a list of rrScan objects
for i in range(len(scanlist)):
    file = os.listdir(dataDir)[i]
#    scn[i] = rr.rrScan()
    scn.append(rr.rrScan())
    scn[i] .importRawFile(dataDir + '//' + file)
    scn[i] .flipTime()
    scn[i] .shiftTime(-125)
    scn[i] .removeDC()
    scn[i] .filterit(cutHigh=0.05)

#
scn = sorted(scn, key=lambda scn: scn.temperature)

for i in range(len(scanlist)-7):

    
    #l = str(scn[i].temperature) + 'K'
    xdata = scn[i].time[0:6230:1]
    ydata = scn[i].ftrace[0:6230:1]
    c = next(color)
    try:
        popt, pcov = scipy.optimize.curve_fit(func, xdata, ydata, p0 = guess)
        #guess = popt
        fitparameters.append(popt)
        
        plt1.plot(xdata, func(xdata, *popt), '--', c=c)
        plt1.plot(scn[i].time,scn[i].ftrace, 
                  c=c, 
                  label=str(scn[i].temperature) + 'K', 
                  alpha=0.5)
#        plt2.plot(par, popt[1])
    except RuntimeError:
        print('no fit parameters found for ' + str(scn[i].temperature) + 'K scan')

    
    

temp = []
tau = []
for i in range(len(scanlist)):
    try:
        tau.append(fitparameters[i][1])
        temp.append(scn[i].temperature)
        #print(temp[i], tau[i])
    except(IndexError):
        print('index error')
    
#plt2.scatter(temp[0:-7:1],tau[0:-7:1],c=colorlist)
plt2.scatter(temp,tau,c=colorlist)
#plt2.set_xscale('log')
def expfunc(x, A, t0, c):
    return A * np.exp(-x / t0) + c
poptRes, pcovRes = scipy.optimize.curve_fit(expfunc, temp[0:-7:1], tau[0:-7:1])
print(poptRes)
xdata = np.linspace(3,70)
plt2.plot(xdata, expfunc(xdata, *poptRes), 'r--')
plt2.text(40,5, str(poptRes[1]))

plt.legend()
fig.show()

#saveDir = filedialog.askdirectory(initialdir = 'E://DATA//RuCl3//Analysis//')
saveDir = 'E://DATA//RuCl3//Analysis//'
savename = saveDir + Title
fig.savefig(savename+'.png', format='png')

