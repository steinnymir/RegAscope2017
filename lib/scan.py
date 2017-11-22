# -*- coding: utf-8 -*-
"""
Created on Nov 22 09:57:29 2017

@author: S.Y. Agustsson
"""

from lib import genericfunctions as gfs
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit, fmin

from matplotlib import cm, pyplot as plt, colorbar
import scipy.io as spio
import scipy.signal as spsignal
import os
import re
import pandas as pd


def main():

    pass



class Scan(object):
    """ Creates an object(OBJ) that contains data and metadata of a single time resolved scan.
    This is used to store all data and metadata while performing a pump-probe measurement.


    Data is contained in:
    - raw_time: raw measured time scale
    - raw_trace: measured data points
    - time: time trace modified by analysis functions applied on this object
    - trace: data modified by analysis functions applied on this object.

    Metadata entries are described in __init__


    """

    def __init__(self):
        """
        """
        self.timeScale = [],
        self.stagePositions = [],
        self.data = {},

        self.currentScan = pd.DataFrame(
            {
            'time': None,
            'currentStagePosition': None,
            'x':None,
            'y': None,
            'Theta': None,
            'r': None,
            'aux1': None,
            'aux2': None,
            'temperature': None
            },
            index = None,
            )

        self.metadata = {
            'date': None,
            'notes': None,
            'sample':{ # information about sample under examination
                'name': None,
                'material': None,
                'notes' : None,
            },
            'pump':{ # parameters of pump beam
                'power': None,
                'spotDiameter': None,
                'energyDensity': None,
                'frequency': None,
                'polarization': None,
            },
            'probe': { # parameters of probe beam
                'power': None,
                'spotDiameter': None,
                'energyDensity': None,
                'frequency': None,
                'polarization': None,
                },
            'destructionPulse': { # parameters of destruction beam
                'power': None,
                'spotDiameter': None,
                'energyDensity': None,
                'frequency': None,
                'polarization': None,
                'delay': None,
            },
            'temperature': None,
            'Notes': None,
            }

    def init_scan(self):
        self.currentScan['stagePosition']

    def build_stage_position(self):
        for i in self.timeScale:
            self.stagePoisitions[i] = self.timeScale[i] * 299.792458

    def store_scan(self, name):

        self.data[name] = self.currentScan


    def reset_current_scan(self):
        del self.currentScan
        self.currentScan = pd.DataFrame(
            {
            'stagePosition': None,
            'time': None,
            'currentStagePosition': None,
            'x':None,
            'y': None,
            'Theta': None,
            'r': None,
            'aux1': None,
            'aux2': None,
            'temperature': None
            },
            index = None,
            )



if __name__ == "__main__":
    main()
