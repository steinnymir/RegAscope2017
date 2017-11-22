# -*- coding: utf-8 -*-
"""
Created on Wed May  3 10:37:19 2017

@author: Steinn Ymir Agustsson
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.optimize import curve_fit, fmin
from matplotlib import cm, colorbar

from lib import genericfunctions as gfs
from lib import redred as rr
from lib.transient import Transient, MultiTransients


def main():
    plt.ioff()  # turn of the damn interactive mode that blocks everything!!

    # select files from folder:
    dir = 'C:/Users/sagustss/py_code/DATA/RuCl3/'
    dir_s = dir
    file_list = gfs.choose_filenames(initialdir=dir)
    # file_list = ('C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_4.8_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_4.32_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_5.0_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_5.5_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_6.0_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_6.5_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_7.0_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_7.5_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_8.0_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_8.5_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_9.0_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_9.5_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_10.0_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_12.0_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_14.0_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_16.0_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_18.0_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_20.0_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_22.0_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_24.0_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_26.0_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_28.0_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_30.0_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_35.0_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_40.0_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_45.0_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_50.0_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_55.0_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_60.0_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_65.0_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_70.0_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_80.0_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_90.0_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_100.0_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_150.0_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_200.0_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_250.0_K.txt',
    #              'C:/data/RuCl3/Low Fluence_temperature/RuCl3_LowFluence_300.0_K.txt')

    print(file_list)

    data = MultiTransients(file_list)
    # data.import_files(file_list)
    data.shift_time(-124.5)
    print(data.series_name, '   ', data.key_parameter)

    # data.quickplot()
    # plt.show()
    #

    ax0, ax1, ax2, colorlist = generate_threeplot_window(data.series_name, data.key_parameter, len(data.transients))

    fig = plt.figure(num=23423, figsize=(16, 12), dpi=80)
    fig.set_figwidth(1024)
    fig.set_figheight(768)
    plt.clf()
    ax = fig.add_subplot(111)
    ax.set_xlabel('Time [ps]', fontsize=18)
    ax.set_ylabel('Differential Reflectivity', fontsize=18)
    ax.set_title('Fitted Scans', fontsize=26)
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)

    initial_values = [0.00005, 0.05, -1, -1]

    popt, pcov, fit_data = data.fit_transients(single_exponential_activation, initial_values, 750, 1, ext_plot=ax, colorlist=colorlist)


    print(fit_data.keys())
    i=0
    for key, item in fit_data.items():

        fig = plt.figure(key)
        ax = fig.add_subplot(111)
        i += 1
        item.quickplot(plt_handle=ax, show=False)
    print(data.transients[0].pump_energy)

    #
    # X = []
    # A = []
    # t = []
    # for key, value in popt.items():
    #
    #     temperature = float(key.split(' ')[0])
    #     if temperature > 70:
    #         pass
    #     else:
    #         X.append(temperature)
    #         A.append(-value[0])
    #         t.append(value[1])
    # # color = iter(colorlist)

    # ax1.plot(X, t, 'o', c='b')
    # ax2.plot(X, A, 'o', c='r')


    # .quickplot(plt_handle=ax, show=False, 'o', c='b')
    # item.quickplot(plt_handle=ax, show=False, 'o', c='b')

    plt.show()
    # ax1.set_xlim([0, 100])
    # ax2.set_xlim([0, 100])
  #  ax1.set_xscale('log')
  #  ax1.set_ylim([0, 0.002])


# define fitting function
def single_exponential_activation(x, A, t0, c, d):
    return A * (1 - np.exp(- x / t0)) + c * x + d


def double_exponential_pos_neg(x, A1, t1, A2, t2, c, d):
    return A1 * (1 - np.exp(- x / t1)) + A2 * (1 - np.exp(- x / t2)) + c * x + d


def print_series_parameters(transient_list):
    print('{0:16} {1:5} {2:5}'.format(' ', 'Pump', 'Probe'))
    print('{0:16}:{1:.3f} {2:.3f}'.format('Energy Density', transient_list[0].pump_energy,
                                          transient_list[0].probe_energy))
    print('{0:16}:{1:.3f} {2:.3f}'.format('Power', transient_list[0].pump_power, transient_list[0].probe_power))


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
    colorlist_length = len(transient_list)
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


def get_data_from_files(file_list, key_parameter=None, description=None):
    temp_list = []
    transients_list = []
    for i, file in enumerate(file_list):
        temp_list.append(tr.Transient(key_parameter=key_parameter, description=description))
        temp_list[i].import_file(file)
        if len(temp_list[i].time) != 0:  # important to not create empty scans in the list of scans
            transients_list.append(temp_list[i])
    print('imported {0} files'.format(len(transients_list)))

    return transients_list


def sort_scan_list_by_parameter(transients_list, parameter):
    sorted_list = sorted(transients_list, key=lambda transients_list: getattr(transients_list, parameter))
    return sorted_list


def generate_threeplot_window(Title, dependence, colorlist_length):
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
    plt2.set_xlabel('Temperature [K]', fontsize=18)
    plt2.set_ylabel('Decay Time [ps]', fontsize=18)
    plt2.set_title('Fit Parameters as function of Temperature', fontsize=20)
    plt2.tick_params(axis='x', labelsize=12)
    plt2.tick_params(axis='y', labelsize=12)
    plt3.set_xlabel('Temperature [K]', fontsize=18)
    plt3.set_ylabel('Amplitude', fontsize=18)
    # plt3.set_title('Amplitude vs Temperature', fontsize=20)

    plt3.tick_params(axis='x', labelsize=12)
    plt3.tick_params(axis='y', labelsize=12)
    colorlist = cm.rainbow(np.linspace(0, 1, colorlist_length))

    return plt1, plt2, plt3, colorlist


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


if __name__ == '__main__':
    main()