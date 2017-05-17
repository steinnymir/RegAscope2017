# -*- coding: utf-8 -*-
"""
Created on Mon May 15 10:18:49 2017

@author: S.Y. Agustsson

module containing general purpouse function

"""
#%% imports
import tkinter as tk
from tkinter import filedialog
from datetime import datetime
import pickle
import os
import re
import scipy.constants as spconst


def main():

    pass


#%% File/Folder Popup dialogs

def chooseFolder(initialdir = 'E://'):
    root = tk.Tk()
    root.withdraw()
    dataDir = filedialog.askdirectory(initialdir = initialdir)
    return(dataDir)

def chooseFile():
    pass

def chooseFiles():
    pass

#%% Generic Utilities

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







#%% naming

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
    parameterDict = {'date' : file_creation_date(file),
                     'material'   : '',
                     'pump_power' : 0,
                     'probe_power': 0,
                     'temperature': 0,
                     'destruction_power': 0,
                     'other': '',
                     }
    #set of possible indicators and corresponding key
    parInd = {'pump_power' : ['pu','pump'],
              'probe_power': ['pr','probe'],
              'temperature': ['t','temp'],
              'destruction_power': ['d','dest','destr'],
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

    parameterDict['material'] = FileName[0:pos]

    return(parameterDict)

class NameString(object):
    """ gets names of variables from whithin a class"""
    def __init__(self,*values):
        for v in values:
            self.__dict__[v] = v

#%% Transformations

def getEnergyDensity(spot_size = 100, power = 1, rep_rate = 283):
    """ Return energy density in microJoule /cm2
        spotsize as FWHM diameter of gaussian profile in micrometers
        power in mW
        reprate in Hz
    """
    area = ((spot_size / 2)**2 *spconst.pi )/ 10**8
    print(power / (rep_rate * area))

def nyqistFreq(timedata):
    """returns the Nyquist frequency from time data"""
    return(abs(0.5 * len(timedata) / timedata[-1] - timedata[0]))

if __name__ == "__main__":
    main()
