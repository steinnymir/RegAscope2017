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



class NameString(object):
    """ gets names of variables from whithin a class"""
    def __init__(self,*values):
        for v in values:
            self.__dict__[v] = v



def getMetadataFromName(filepath):
    ''' interprets name of file and returns a dictionary with all contained metadata
    it requires a name with elements separated by ['-','_',',',' '].
    example: RuCl3_pu_15mW_pr_5mW_t_4.5K.mat
             RuCl3-pu-15mW-pr-5mW-t-4.5K.mat

    replacement of old name_to_info'''
    # transform to more usable string form

    metadataDict = {}
    metadataDict['date'] = file_creation_date(filepath) # get file creation date

    FileName = ('.').join(os.path.basename(filepath).split('.')[:-1])
    filename = FileName.lower().replace(',','.')

    #identify most recurring possible separator from separators list in filename
    separators = ['-','_',',',' ']
    sepCount = {}
    for item in separators:
        sepCount[item] = filename.count(item)

    inverse = [(value, key) for key, value in sepCount.items()]
    sep = max(inverse)[1]
    #generate a list of elements as separated by sep


    metadata_names = {'pump_power' : ['pu','pump'],
                      'probe_power': ['pr','probe'],
                      'destruction_power': ['d','dest','destr'],
                      'temperature': ['t','temp'],
                      'pump_polarization': ['pupol'],
                      'probe_polarization': ['prpol'],
                      'destruction_polarization': ['dpol'],
                      'sample_orientation': ['sor'],
              }
#    metadata_names_inv =  [(value, key) for key, value in metadata_names.items()]
#    print(metadata_names_inv)
    parameter_list = []

    if sepCount[sep] > 1: # use separator based interpreter
        filename_list = filename.split(sep)
        for key, item in metadata_names.items():
            for parameter in item:

                try:
                    # look for an identifier as defined in metadata_names and
                    # if there, get the value in the next item, that should be
                    # the value corresponding to such parameter
                    parameter_index = filename_list.index(parameter)

                    value = re.findall("\d+\.\d+", filename_list[parameter_index + 1])
                    parameter_list.append(parameter)
                    # save retrieved parameter value
                    metadataDict[key] = float(value[0])
                    filename_list.pop(parameter_index) # twice because first removes
                    filename_list.pop(parameter_index) # string, second the value
                except :
                    pass
        metadataDict['material'] = filename_list.pop(0)
        metadataDict['other'] = filename_list

    else: # use 'parameter name search' based interpreter
        parameter_list = []
        for key, item in metadata_names.items():
            for parameter in item:

                if parameter in filename:
                    parameter_list.append(parameter)
                    try:
                        # Partition string around parameter delimiter and pick the
                        # first following float number
                        parameter_list.append(parameter)
                        value = float(re.findall(r"\d+\.\d*",
                                             filename.partition(parameter)[2])[0])
                        metadataDict[key] = value
                    except IndexError:
                        pass
        pos = 100
        for string in parameter_list:
            x = filename.find(string)
            if x < pos:
                pos=x

        metadataDict['material'] = FileName[0:pos]
    print(filename)

    return(metadataDict)





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
