# -*- coding: utf-8 -*-
"""
Created on Wed May  3 22:03:55 2017

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


root = tk.Tk()
root.withdraw()
dataDir = filedialog.askdirectory(initialdir = 'E://DATA//RuCl3//Analysis//')

scanlist = os.listdir(dataDir)

Title = os.path.basename(dataDir)

fig = plt.figure(Title, figsize = (19,10))
plt.clf()
plt1 = fig.add_subplot(121)
plt2 = fig.add_subplot(222)
plt3 = fig.add_subplot(224)
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
plt3.set_xlabel('Temperature [K]', fontsize=18)
plt3.set_ylabel('Amplitude', fontsize=18)
plt3.set_title('Decay time vs Temperature',fontsize=26)
plt3.tick_params(axis='x', labelsize=12)
plt3.tick_params(axis='y', labelsize=12)
colorlist = cm.rainbow(np.linspace(0,1,5))
color=iter(colorlist)




#fig2 = plt.figure('decay times')
#plt2 = fig2.add_subplot(111)


scn = []
#define fitting function
def func(x, A, t0, c, d):
    return A * np.exp(- x / t0) + c*x +d
def func1(x, A, t0, d):
    return A * np.exp(- x / t0) + d
guess = [0.0005, 0.1, -1, -1]
guess1 = [0.0005, 0.1, -1]
fitparameters = []


#scan through folder and create a list of rrScan objects
for i in range(len(scanlist)):
    file = os.listdir(dataDir)[i]
#    scn[i] = rr.rrScan()
    scn.append(rr.rrScan())
    scn[i] .importCSV(dataDir + '//' + file)

#    scn[i] .filterit(cutHigh=0.05)
#    scn[i] .flipTime()

    scn[i].shiftTime(-28)
#    scn[i] .removeDC()

#
scn = sorted(scn, key=lambda scn: scn.pumpPw)

for i in range(len(scanlist)-3):

    
    #l = str(scn[i].temperature) + 'K'
    xdata = scn[i].time[0:4985:1]
    ydata = scn[i].trace[0:4985:1]
    c = next(color)
    try:
        popt, pcov = scipy.optimize.curve_fit(func, xdata, ydata, p0 = guess)
        #guess = popt
        fitparameters.append(popt)
        
        plt1.plot(xdata, func(xdata, *popt), '--', c=c)
        plt1.plot(scn[i].time[250::],scn[i].trace[250::], 
                  c=c, 
                  label=str(scn[i].temperature) + 'K', 
                  alpha=0.5)
#        plt2.plot(par, popt[1])
    except RuntimeError:
        print('no fit parameters found for ' + str(scn[i].temperature) + 'K scan')

    
    

power = []
tau = []
amp= []
for i in range(len(scanlist)):
    try:
        tau.append(fitparameters[i][1])
        amp.append(fitparameters[i][0])
        power.append(scn[i].pumpPw)

        #print(temp[i], tau[i])
    except(IndexError):
        print('index error')
    
#plt2.scatter(temp[0:-7:1],tau[0:-7:1],c=colorlist)
plt2.scatter(power,tau,c=colorlist)
plt3.scatter(power,amp,c=colorlist)
plt3.set_xscale('log')
plt3.set_ylim([0,0.002])



#fit temperature
def expfunc(t, A, t0, c):
    return A * np.exp(- t / t0) + c

poptDT, pcovDT = scipy.optimize.curve_fit(expfunc, power[0:-7:1], tau[0:-7:1], p0 = [1, 10,1])

xdata = np.linspace(3,70)

plt2.plot(xdata, expfunc(xdata, *poptDT), 'r--')
#plt2.text(40,5, str(poptDT[1]))

#poptDT, pcovDT = scipy.optimize.curve_fit(expfunc, temp[0:-7:1], tau[0:-7:1], p0= [0.001,)


#fit amplitude



plt.legend()
fig.show()
#
##saveDir = filedialog.askdirectory(initialdir = 'E://DATA//RuCl3//Analysis//')
#saveDir = 'E://DATA//RuCl3//Analysis//' + Title + '//'
#savename = saveDir + Title
#if not os.path.exists(saveDir):
#    os.makedirs(saveDir)


#fig.savefig(savename +'.png', format='png')

