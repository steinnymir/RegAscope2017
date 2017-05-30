# -*- coding: utf-8 -*-
"""
Created on Wed May  3 10:37:19 2017

@author: sagustss
"""

from functionlibrary import transient as tr
from functionlibrary import genericfunctions as gfs
from functionlibrary import redred as rr

import os
import numpy as np
import scipy
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog
from matplotlib import cm


def main():
    # select files from folder:
    dir = 'C:/Users/sagustss/py_code/DATA/_RAW/RuCl3/2017-04-19/'
    dir_s = 'C:/Users/sagustss/py_code/DATA/'
    file_list = gfs.choose_filenames(initialdir=dir)

    title = 'Low Fluence'
    dependence = 'temperature'
    pump_spot = 97
    probe_spot = 64


    transient_list = get_data_from_files(file_list)  # import data
    for item in transient_list:  # iterate and correct scans
        item.description = title
        item.key_dependence = dependence
        item.pump_spot = pump_spot
        item.probe_spot = probe_spot
        item.calc_energy_densities()
        item.give_name()  # todo: make the naming procedure simpler and more 'working'

        item.clean_data(flipTime=True, shiftTime=124.5, filterLowPass=0.005)
    print_series_parameters(transient_list)






    transient_list = sort_scan_list_by_parameter(transient_list, dependence)

    quickplot_list(transient_list, title, dependence)
    plt.show()

    save_dir = gfs.choose_folder(initialdir=dir)
    save_dir += '/' + title + ' - ' + dependence + 'dependence'

    for item in transient_list:
        item.export_file_csv(save_dir)

def print_series_parameters(transient_list):
    print('\t\t\t\t Pump\tProbe\n----------------------------')
    print('energy density:\t {0:3f}\t{1:3f}'.format(transient_list[0].pump_energy, transient_list[0].probe_energy))
    print('power:\t\t\t {0}\t{1}'.format(transient_list[0].pump_power, transient_list[0].probe_power))


# define fitting function
def func(x, A, t0, c, d):
    return A * np.exp(- x / t0) + c * x + d


def quickplot_list(transient_list, title, dependence):
    """ simple plot of a list of transients """  # todo: move to transients.py -> under multitransients()
    fig = plt.figure(num=1)
    plt.clf()
    ax = fig.add_subplot(111)
    ax.set_xlabel('Time [ps]', fontsize=18)
    ax.set_ylabel('Differential Reflectivity', fontsize=18)
    ax.set_title('Fitted Scans', fontsize=26)
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)

    # todo: make nice color iteration, that follows parameter value
        # using gcd didnt work, or the function for it is buggy
    # colorlist_length, color_step = get_parameter_min_max_minstep(transient_list, dependence)
    colorlist_length =  len(transient_list)
    colorlist = cm.rainbow(np.linspace(0, 1, colorlist_length))
    color = iter(colorlist)

    for curve in transient_list:
        # l = str(scn[i].temperature) + 'K'
        xdata = curve.time
        ydata = curve.trace
        parameter = getattr(curve, dependence)
        col = next(color)
        ax.plot(xdata, ydata, c=col, label=str(curve.temperature) + 'K', alpha=0.5)
    return fig


def get_parameter_min_max_minstep(transient_list, dependence):
    par_values = []
    for item in transient_list:  # error
        par_values.append(getattr(item, dependence))
    par_values, size = sorted(par_values), len(par_values)
    res = [par_values[i + 1] - par_values[i] for i in range(size) if i + 1 < size]
    step = min(res)
    minPar = min(par_values)
    maxPar = max(par_values)

    gcd = gcd_list(par_values)
    print(par_values)
    print(gcd)
    even_list_len = (maxPar - minPar) / gcd
    print(even_list_len)
    return even_list_len, step


def gcd(a, b):
    if (b == 0):
        return a
    else:
        return gcd(b, a % b)


def gcd_list(A):
    res = A[0]
    for c in A[1::]:
        res = gcd(res, c)
    return res



def get_data_from_files(file_list):
    temp_list = []
    transients_list = []
    for i, file in enumerate(file_list):
        temp_list.append(tr.Transient())
        temp_list[i].import_file(file)
        if len(temp_list[i].time) != 0:  # important to not create empty scans in the list of scans
            transients_list.append(temp_list[i])
    print('imported {0} files'.format(len(transients_list)))

    return transients_list


def sort_scan_list_by_parameter(transients_list, parameter):
    sorted_list = sorted(transients_list, key=lambda transients_list: getattr(transients_list, parameter))
    return sorted_list


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################


def crap():
    scn = sorted(scn, key=lambda scn: scn.temperature)

    for i in range(len(scanlist) - 7):

        # l = str(scn[i].temperature) + 'K'
        xdata = scn[i].time[0:6230:1]
        ydata = scn[i].trace[0:6230:1]
        c = next(color)
        try:
            popt, pcov = scipy.optimize.curve_fit(func, xdata, ydata, p0=guess)
            # guess = popt
            fitparameters.append(popt)

            plt1.plot(xdata, func(xdata, *popt), '--', c=c)
            plt1.plot(scn[i].time[50::], scn[i].trace[50::],
                      c=c,
                      label=str(scn[i].temperature) + 'K',
                      alpha=0.5)
        # plt2.plot(par, popt[1])
        except RuntimeError:
            print('no fit parameters found for ' + str(scn[i].temperature) + 'K scan')

    temp = []
    tau = []
    amp = []
    for i in range(len(scanlist)):
        try:
            tau.append(fitparameters[i][1])
            amp.append(fitparameters[i][0])
            temp.append(scn[i].temperature)

            # print(temp[i], tau[i])
        except(IndexError):
            print('index error')

    # plt2.scatter(temp[0:-7:1],tau[0:-7:1],c=colorlist)
    plt2.scatter(temp, tau, c=colorlist)
    plt3.scatter(temp, amp, c=colorlist)
    plt3.set_xscale('log')
    plt3.set_ylim([0, 0.002])

    # fit temperature
    def expfunc(t, A, t0, c):
        return A * np.exp(- t / t0) + c

    poptDT, pcovDT = scipy.optimize.curve_fit(expfunc, temp[0:-7:1], tau[0:-7:1], p0=[1, 10, 1])

    xdata = np.linspace(3, 70)

    plt2.plot(xdata, expfunc(xdata, *poptDT), 'r--')
    # plt2.text(40,5, str(poptDT[1]))

    # poptDT, pcovDT = scipy.optimize.curve_fit(expfunc, temp[0:-7:1], tau[0:-7:1], p0= [0.001,)


    # fit amplitude



    plt.legend()
    fig.show()

    # saveDir = filedialog.askdirectory(initialdir = 'E://DATA//RuCl3//Analysis//')
    saveDir = 'E://DATA//RuCl3//Analysis//' + Title + '//'
    savename = saveDir + Title
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)

    for i in range(len(scn)):
        scn[i].export_file_csv(saveDir)
    fig.savefig(savename + '.png', format='png')

    # fig2 = plt.figure('decay times')
    # plt2 = fig2.add_subplot(111)


    scn = []

    # define fitting function
    def func(x, A, t0, c, d):
        return A * np.exp(- x / t0) + c * x + d

    guess = [0.0005, 0.1, -1, -1]
    fitparameters = []

    # scan through folder and create a list of rrScan objects
    for i in range(len(scanlist)):
        file = os.listdir(dataDir)[i]
        #    scn[i] = rr.rrScan()
        scn.append(rr.rrScan())
        scn[i].importRawFile(dataDir + '//' + file)
        scn[i].filterit(cutHigh=0.05)
        scn[i].flipTime()
        scn[i].shift_time(-125)
        scn[i].remove_DC_offset()

    #
    scn = sorted(scn, key=lambda scn: scn.temperature)

    for i in range(len(scanlist) - 7):

        # l = str(scn[i].temperature) + 'K'
        xdata = scn[i].time[0:6230:1]
        ydata = scn[i].trace[0:6230:1]
        c = next(color)
        try:
            popt, pcov = scipy.optimize.curve_fit(func, xdata, ydata, p0=guess)
            # guess = popt
            fitparameters.append(popt)

            plt1.plot(xdata, func(xdata, *popt), '--', c=c)
            plt1.plot(scn[i].time[50::], scn[i].trace[50::],
                      c=c,
                      label=str(scn[i].temperature) + 'K',
                      alpha=0.5)
        # plt2.plot(par, popt[1])
        except RuntimeError:
            print('no fit parameters found for ' + str(scn[i].temperature) + 'K scan')

    temp = []
    tau = []
    amp = []
    for i in range(len(scanlist)):
        try:
            tau.append(fitparameters[i][1])
            amp.append(fitparameters[i][0])
            temp.append(scn[i].temperature)

            # print(temp[i], tau[i])
        except(IndexError):
            print('index error')

    # plt2.scatter(temp[0:-7:1],tau[0:-7:1],c=colorlist)
    plt2.scatter(temp, tau, c=colorlist)
    plt3.scatter(temp, amp, c=colorlist)
    plt3.set_xscale('log')
    plt3.set_ylim([0, 0.002])

    # fit temperature
    def expfunc(t, A, t0, c):
        return A * np.exp(- t / t0) + c

    poptDT, pcovDT = scipy.optimize.curve_fit(expfunc, temp[0:-7:1], tau[0:-7:1], p0=[1, 10, 1])

    xdata = np.linspace(3, 70)

    plt2.plot(xdata, expfunc(xdata, *poptDT), 'r--')
    # plt2.text(40,5, str(poptDT[1]))

    # poptDT, pcovDT = scipy.optimize.curve_fit(expfunc, temp[0:-7:1], tau[0:-7:1], p0= [0.001,)


    # fit amplitude



    plt.legend()
    fig.show()

    # saveDir = filedialog.askdirectory(initialdir = 'E://DATA//RuCl3//Analysis//')
    saveDir = 'E://DATA//RuCl3//Analysis//' + Title + '//'
    savename = saveDir + Title
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)

    for i in range(len(scn)):
        scn[i].export_file_csv(saveDir)
    fig.savefig(savename + '.png', format='png')


def generate_threeplot_window(Title, dependence):
    fig = plt.figure(Title, figsize=(19, 10))
    plt.clf()
    plt1 = fig.add_subplot(121)
    plt2 = fig.add_subplot(222)
    plt3 = fig.add_subplot(224)
    plt1.set_xlabel('Time [ps]', fontsize=18)
    plt1.set_ylabel('Differential Reflectivity', fontsize=18)
    plt1.set_title('Fitted Scans', fontsize=26)
    plt1.tick_params(axis='x', labelsize=12)
    plt1.tick_params(axis='y', labelsize=12)
    plt2.set_xlabel(dependence, fontsize=18)
    plt2.set_ylabel('Decay Time [ps]', fontsize=18)
    plt2.set_title('Decay time vs Temperature', fontsize=26)
    plt2.tick_params(axis='x', labelsize=12)
    plt2.tick_params(axis='y', labelsize=12)
    plt3.set_xlabel('Temperature [K]', fontsize=18)
    plt3.set_ylabel('Amplitude', fontsize=18)
    plt3.set_title('Decay time vs Temperature', fontsize=26)
    plt3.tick_params(axis='x', labelsize=12)
    plt3.tick_params(axis='y', labelsize=12)
    colorlist = cm.rainbow(np.linspace(0, 1, len(os.listdir(dataDir))))
    color = iter(colorlist)

    return (plt1, plt2, plt3, color)


if __name__ == '__main__':
    main()





    # print('exited')
    # guess = [0.0005, 0.1, -1, -1]
    # fitparameters = []
