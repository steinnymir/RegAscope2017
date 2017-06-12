# -*- coding: utf-8 -*-
"""
Created on Mon May 15 09:57:29 2017

@author: S.Y. Agustsson
"""

from functionlibrary import genericfunctions as gfs
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit, fmin

from matplotlib import cm, pyplot as plt, colorbar
import scipy.io as spio
import scipy.signal as spsignal
import os
import re
import pandas


def main():
    files = os.listdir('C:/data/RuCl3/Low Fluence_temperature')
    filepaths = []
    for name in files:
        filepaths.append('C:/data/RuCl3/Low Fluence_temperature/' + name)

    series = MultiTransients(transients_list=filepaths, series_name='Low Fluence - temperature dependence',
                             key_parameter='temperature')
    print(series.key_parameter)
    series.quickplot()

    plt.show()


class Transient(object):
    """ Creates an object(OBJ) that contains data and metadata of a single time resolved scan. The standard scan is that
    obtained from a single color pump probe setup, aka 'RedRed'.
    Data can be imported through use of OBJ.import_file(filepath), which imports either the raw .mat format outputted by RegaScope2012,
    or the .txt format outputted by this class (called csv here due to the csv format of contained data). To output such
    .txt use OBJ.export_csv(directory).

    Data is contained in:
    - raw_time: raw measured time scale
    - raw_trace: measured data points
    - time: time trace modified by analysis functions applied on this object
    - trace: data modified by analysis functions applied on this object.

    Metadata entries are described in __init__


    """

    def __init__(self, key_parameter=None, series_name=None, description=None):
        """
        :param key_parameter: str, name of parameter iterated in the series this scan belongs to.
        :param description: description of the scan/series.
        """

        ######################################
        #                Data                #
        ######################################

        self.raw_time = np.array([])  # time data
        self.raw_trace = np.array([])  # raw trace data

        self.time = np.array([])  # cleaned time axis
        self.trace = np.array([])  # cleaned and modified data trace
        self.description = description
        self.key_parameter = key_parameter
        if self.key_parameter is not None:
            self.key_parameter_value = getattr(self, self.key_parameter)
        else:
            self.key_parameter_value = None
        self.series_name = series_name
        # ignore list for metadata export. Add here any further non-metadata attributes created in this class.
        self.DATA_ATTRIBUTES = ('raw_time', 'raw_trace', 'time', 'trace', 'DATA_ATTRIBUTES')

        ######################################
        #              Metadata              #
        ######################################

        self.name = None  # String used as file name for saving data
        self.material = None  # Material name
        self.date = None  # Scan date in format YYYY-MM-DD hh.mm.ss
        self.original_filepath = None  # Path to original raw file

        # parameters

        self.pump_power = None  # Pump Power [mW]
        self.probe_power = None  # Probe Power [mW]
        self.destruction_power = None  # Destruction Power [mW]

        # Spot size represents the FWHM diameter from Gaussian fit of the beam profile
        self.pump_spot = None  # Pump beam spot size [micrometers]
        self.probe_spot = None  # Probe beam spot size [micrometers]
        self.destruction_spot = None  # Destruction beam spot size [micrometers]

        # Excitation densities calculated from power, spot size and repetition rate
        self.pump_energy = None
        self.probe_energy = None
        self.destruction_energy = None

        # destruction pulse parameters:
        self.t12 = None  # delay between pump and destruction pulses, [ps]

        # Polarization are measured clockwise in propagation direction of the beam, 0 = 12o'clock
        self.pump_polarization = None  # Pump beam polarization [deg]
        self.probe_polarization = None  # Probe beam polarization [deg]
        self.destruction_polarization = None  # Destruction beam polariz. [deg]
        self.sample_orientation = None  # Sample orientation [deg]
        self.temperature = None  # Temperature [K]
        self.R0 = None  # Static reflectivity

        ######################################
        #              Analysis              #
        ######################################

        self.analysis_log = {}  # Keeps track of analysis changes performed
        self.fit_function = None
        self.fit_parameters = None





        ######################################
        #             input info             #
        ######################################

    # %% metadata management
    def key_parameter_value(self):
        try:
            self.key_parameter_value = getattr(self, str(self.key_parameter))
            return self.key_parameter_value
        except AttributeError:
            raise AttributeError('No key parameter detected.')

    def calc_energy_densities(self, rep_rate=283000):
        """ recalculate metadata depending on given parameters.
            it calculates energy densities
        """
        beams = ['pump', 'probe', 'destruction']

        for beam in beams:
            if getattr(self, (beam + '_spot')) is None:
                pass
            else:
                power = getattr(self, (beam + '_power'))
                spot = getattr(self, (beam + '_spot'))
                if beam == 'pump':
                    rep_rate = rep_rate / 2  # pump has half reprate (darkcontrol)
                energy = round(gfs.get_energy_density(spot, power, rep_rate), 3)
                print(energy)
                setattr(self, (beam + '_energy'), energy)

    def input_attribute(self, attribute_name, value):
        """
        manually input values for metadata attributes
        :param attribute_name: name of parameter or attribute
        :param value: value to assign to parameter
        """
        setattr(self, attribute_name, value)

    def get_metadata(self):
        """ Returns a dictionary containing all metadata information available
            keys are attribute names and values the corresponding value.
            :rtype: dict
        """
        metadata = {'analysis_log': {}}  # make dict for metadata. Also create analysis_log entry, used later
        attributes = self.__dict__
        for key, value in attributes.items():
            if key not in self.DATA_ATTRIBUTES and value is not None:  # ignore time, trace and all other fields defined in data_attributes
                try:
                    # if parameter is a number != 0 append to metadata
                    metadata[key] = float(attributes[key])
                except (TypeError, ValueError):
                    metadata[key] = value
        return metadata

    def log_it(self, keyword, overwrite=False, *args, **kargs):
        """
        Generate log entry for analysis_log.
            creates a key with given key in analysis_log, making it:
                - boolean if no other args or kargs are given, flips previous
                values written in log
                - list if *args are passed
                - dictionary if **kargs are passed
            if overwrite is False, it appends values on previous logs,
            if True, it obviously overwrites them.
            :param key: string with name of analysis function used
            :type overwrite: bool
        """

        # key_string = 'analysis_log[' + key + ']'

        # make the right type of entry for the log
        if kargs or args:
            if kargs:  # make a dictionary of parameters passed
                entry = {}
                for key in kargs:
                    entry[key] = kargs[key]
            if args:  # make a lsit of parameter passed
                entry = []
            for arg in args:
                entry.append(arg)
        else:  # make None to trigger boolean behaviour
            entry = None
        # Check if previous logs with this key and eventually overwrite/append
        # the new entry.
        try:
            # prev = getattr(self, key_string)
            previous_value = self.analysis_log[keyword]
            # if it finds such a key in this library, it overwrites/appends new information.
            if entry == None:  # trigger boolean behaviour, flipping previous registered status if available
                self.analysis_log[keyword] = not previous_value

            elif entry is list:
                if overwrite:
                    self.analysis_log[keyword] = entry
                    # setattr(self, key_string, previous_value + entry)
                else:
                    entry = previous_value.append(entry)
                    self.analysis_log[keyword] = entry
            elif entry is dict:
                if overwrite:
                    self.analysis_log[keyword] = entry
                else:
                    new_entry = {}
                    for key in entry:
                        if key in previous_value:
                            new_entry[key] = previous_value[key].append(entry[key])
                        else:
                            new_entry[key] = entry[key]

                    self.analysis_log[keyword] = new_entry
        except KeyError:  # rises Key error when key was not previously assigned -> no previous record of this analysis
            if entry == None:  # if boolean, its it's first usage, therefore set it to true
                self.analysis_log[keyword] = True
            else:  # otherwise just append the list/dictionary.
                self.analysis_log[keyword] = entry

    def give_name(self):
        """Define name attribute as material_date."""

        if self.key_parameter is None:
            self.key_parameter = input('What is the Key parameter for basename?  ')
        if self.description is None:
            self.description = input('Add brief description for file name:  ')

        self.name = (str(self.material) + '_' +
                     str(self.description) + '_' +
                     str(getattr(self, self.key_parameter)) + '_' +
                     str(self.get_unit(self.key_parameter)))

    def get_unit(self, parameter):
        """ Returns the unit of the given parameter.
        works for
            - power
            - polarization
            - R0
            - trace
            - time
            - energy
            - tempeature
            parameter:type parameter: str
            """
        splitpar = parameter.split('_')
        if splitpar[-1] == 'power':
            return ('mW')
        elif splitpar[-1] == 'polarization' or splitpar[-1] == 'orientation':
            return ('deg')
        elif splitpar[-1] == 'R0':
            return ('V')
        elif splitpar[-1] == 'trace':
            return ('')
        elif splitpar[-1] == 'time':
            return ('ps')
        elif splitpar[-1] == 'energy':
            return ('mJ/cm^2')
        elif splitpar[-1] == 'temperature':
            return ('K')
        else:
            return ('')

    # %% import export

    def import_file(self, filepath, cleanData=True, key_parameter=None, description=None):
        """Imports a file, csv or .mat
        uses self.import_file_mat() and self.import_file_csv,
        depending on file extension"""

        try:  # if it finds the file requested...
            ext = os.path.splitext(filepath)[-1].lower()
            basename = os.path.basename(filepath)
            if ext == '.mat':
                if basename.lower() != 't-cal.mat':
                    self.import_file_mat(filepath)
                else:
                    print('Ignored t-cal.mat')

            elif ext == '.txt':
                self.import_file_csv(filepath)
            else:
                print("Invalid format. Couldn't import file: " + filepath)
            # self.importMetadata()
            if key_parameter is not None:
                self.key_parameter = key_parameter
            if description is not None:
                self.description = description
            self.give_name()
            print('Imported {0} as {1}'.format(basename, self.name))
            if cleanData and len(self.raw_time) != 0:
                self.clean_data()

        except FileNotFoundError:
            print('File ' + filepath + ' not found')

    def import_file_mat(self, filepath):
        """Import data from a raw .mat file generated by redred software.
        extracts data about raw_time raw_trace and R0.

        """
        self.original_filepath = filepath
        data = spio.loadmat(filepath)
        try:  # if it finds the right data structure
            self.raw_time = data['Daten'][2]
            self.raw_trace = data['Daten'][0]

            self.R0 = data['DC'][0][0]
            # get all metadata from name
            metadataDict = gfs.get_metadata_from_name(filepath)  # todo: add eventual non scripted parameters
            # write metadata to relative attributes
            for key in metadataDict:
                try:
                    setattr(self, key, metadataDict[key])
                except KeyError:
                    print('invalid key: ' + key)

        except KeyError:
            print(filepath + ' is not a valid redred scan datafile')

    def import_file_csv(self, filepath):
        """
        Import data from a .txt file containing metadata in the header.

        Metadata should be coded as variable names from this class:
            material, date, pump_power, temperature, probe_polarization etc...
        Data expected is 4 couloms: raw_time, raw_trace, time, trace.

        filepath should full path to file as string.

        """
        # ---------- get metadata ----------

        # dictionary of attributes where to assign parameters
        attributes = self.__dict__
        parameters = []
        for attribute in attributes:  # use only non-data attributes
            if attribute not in self.DATA_ATTRIBUTES:
                parameters.append(attribute)

        with open(filepath, 'r') as f:
            n = 0
            for l in f:
                # search for the data coulomn header
                if 'raw_time' in l:
                    dataOffset = n + 1  # used for offset of data fetching, data starts from next line -> +1
                    columnHeaders = l.replace('\n', '').replace(' ', '').split(',')
                else:
                    n += 1
                # split each line from file into a list
                word = l[:-1:].split(': ')
                # if the first word corresponds to an attribute name
                if word[0] in parameters:
                    key = word[0]
                    value_string = word[1].replace(' ', '')
                    if self.get_unit(key):
                        # if parameter expects units, get only numbers,
                        value = float(re.findall("\d+\.\d+", value_string)[0])
                    else:  # otherwise get the whole string
                        value = value_string
                    # create/assign attribute from imported parameter
                    setattr(self, key, value)
        self.key_parameter_value = getattr(self, self.key_parameter)

        # ---------- get data ---------- using pandas! :)

        data = pandas.read_csv(filepath, names=columnHeaders, skiprows=dataOffset)
        for col in data.columns:
            # make a list of float tipe data for each dataset found
            col_data = getattr(data, col).astype(float).tolist()
            data_list = []
            for i in col_data:
                if i != 'nan': data_list.append(i)

            setattr(self, col, data_list)

    def export_file_csv(self, directory):
        """
        save Transient() to a .txt file in csv format (data)
        Metadata header is in tab separated values, generated as 'name': 'value' 'unit'
        data is comma separated values, as raw_time, raw_trace, time, trace.

        Metadata is obtained from get_metadata(), resulting in all non0 parameters available.
        """
        # ----------- metadata -----------


        # self.give_name()  # creates a name for the file
        print('Exporting {0}'.format(self.name))
        metadata = self.get_metadata()
        logDict = metadata.pop('analysis_log', None)  # separate log, will be printed later
        logDict.pop('', None)  # remove the useless empty entry
        name = metadata.pop('name', None)
        original_filepath = metadata.pop('original_filepath', None)

        # open file with name self.name in overwrite mode
        file = open(directory + name + '.txt', 'w+')
        # make a title in the header
        file.write('RedRed Scan\n\nMetadata\n\n')
        # write metadata as: parameter: value unit
        for key in metadata:
            try:
                line = (key + ': ' +
                        str(metadata[key]) + ' ' +
                        self.get_unit(key) + '\n')
                file.write(line)
            except TypeError:
                print("Type error for " + key + 'when writing to file: ' +
                      self.name)
        # write analysis log as function: values
        file.write('\nAnalysis\n')
        for key in logDict:
            line = key + ': ' + str(logDict[key]) + '\n'
            file.write(line)

            # ----------- Data -----------
            # Data header followed by column heads:
        file.write('\n\nData\n\n')
        file.write('raw_time, raw_trace, time, trace\n')

        for i in range(len(self.raw_time)):
            line = str(self.raw_time[i]) + ',' + str(self.raw_trace[i])
            try:  # try appending analysied data, or skip if it is finished
                # this because of the deleting of initial and final data
                line += ',' + str(self.time[i]) + ',' + str(self.trace[i])
            except IndexError:
                pass
            finally:
                line += '\n'

            file.write(line)
        file.close()

    # %% Data manipulation

    def clean_data(self, cropTimeScale=True, shiftTime=0, flipTime=True, removeDC=True, filterLowPass=True,
                   flipTrace=False):
        """Perform a standard set of data cleaning, good for quick plotting and test purposes."""
        if cropTimeScale:
            self.crop_time_scale()
        if shiftTime:
            self.shift_time(shiftTime)
        if filterLowPass:
            self.filter_low_pass()
        if flipTrace:
            self.flip_trace()
        if removeDC:
            self.remove_DC_offset()
        if flipTime:
            self.flip_time()

    def crop_time_scale(self):  # todo: fix the overwriting issue
        """chops time scale to the monotonous central behaviour, deleting the wierd ends.
        ATTENTION: overwrites self.time and self.trace, deleting any previous changes"""

        # clear previous time and trace, and the analysis log since it goes lost
        self.analysis_log = {}  # reset log
        self.time = []  # reset time and trace
        self.trace = []
        # print('crop time scale, len: ' + str(len(self.raw_time)))
        maxT = max(self.raw_time)
        minT = min(self.raw_time)

        # pick the type of behaviour. redred scans will always be false in this statement
        if self.raw_time[0] < self.raw_time[1]:
            start = 0
            while self.raw_time[start] < maxT:
                start += 1
            end = start
            while self.raw_time[end] > minT:
                end += 1
            i = 0
            while i in range(end - start):
                self.time.append(self.raw_time[i + start])
                self.trace.append(self.raw_trace[i + start])
                i += 1
        elif self.raw_time[0] > self.raw_time[1]:
            start = 0
            while self.raw_time[start] > minT:
                start += 1
            end = start
            while self.raw_time[end] < maxT:
                end += 1
            i = 0
            while i in range(end - start):
                self.time.append(self.raw_time[i + start])
                self.trace.append(self.raw_trace[i + start])
                i += 1
        self.log_it('Crop Time Scale', maxtime=maxT, mintime=minT)

    def shift_time(self, tshift):
        """ Shift time scale by tshift. Changes time zero
        writes to analysis_log the shifted value, or increases it if already present"""
        self.time = np.array(self.time) - tshift
        self.log_it('Shift Time', tshift)

    def flip_time(self):
        """ Flip time scale: t = -t
        also reverts order in the array"""
        self.time = self.time[::-1]
        self.time = -np.array(self.time)
        self.trace = self.trace[::-1]
        self.log_it('Flip Time')

    def flip_trace(self):
        """ Flip the Y trace, usually not needed from matlab redred software"""
        self.trace = -self.trace
        self.log_it('Flip Trace')

    def remove_DC_offset(self, window=40):  # todo: change range in case of flipped scan!!!
        """Remove DC offset.
        offset is caluclated with 40 points (~700fs) taken at negative time delays.
        such delay is at the end of the scan in raw data, or at the beginning
        if scan was reverted by flip_time """
        try:
            reverted = self.analysis_log['Flip Time']
        except KeyError:
            reverted = False

        if reverted:
            shift = np.average(self.trace[0:window:1])
        else:
            tpoints = len(self.time)
            shift = np.average(self.trace[tpoints - window:tpoints:1])

        self.trace = self.trace - shift
        self.log_it('Remove DC', window=window, shift=shift)

    def filter_low_pass(self, cutHigh=0.1, order=1,
                        return_frequency=False):  # todo: add different methods between which to choose
        """ apply simple low pass filter to data. if return_frequency is True, returns the filter frequency value
        in THz ( if time data is in ps)

        This function applies a linear filter twice, once forward and once backwards.
        The combined filter has linear phase.
        To avoid spikes at edges of the scan, Gustaffson's method is used:
         F. Gustaffson, “Determining the initial states in forward-backward filtering”,
         Transactions on Signal Processing, Vol. 46, pp. 988-992, 1996.
         """

        b, a = spsignal.butter(order, cutHigh, 'low', analog=False)
        self.trace = spsignal.filtfilt(b, a, self.trace, method='gust')
        frequency = gfs.get_nyquist_frequency(self.time) * cutHigh
        self.log_it('Low Pass Filter', frequency=frequency, nyq_factor=cutHigh, order=order)
        if return_frequency:
            return frequency

    def normalize_to_parameter(self, parameter):
        """ Normalize scan by dividing by its pump power value"""
        if getattr(self, parameter):
            if getattr(self, parameter) != 0:
                self.trace = self.trace / getattr(self, parameter)
        else:
            print('Normalization failed: invalid parameter name')
        logkey = 'Normalized by ' + parameter.replace('_', ' ')
        self.log_it(logkey)

    def quickplot(self, xlabel='Time [ps]', ylabel='Trace', fntsize=15, title='Transient', clear=False, raw=False):
        """Generates a quick simple plot with matplotlib """
        if clear: plt.clf()
        quickplotfig = plt.figure(num=1)
        ax = quickplotfig.add_subplot(111)
        if raw:
            ax.plot(self.raw_time, self.raw_trace, 'o')

        else:
            ax.plot(self.time, self.trace, )
        ax.set_xlabel(xlabel, fontsize=fntsize)
        ax.set_ylabel(ylabel, fontsize=fntsize)
        ax.set_title(title, fontsize=fntsize)
        ax.tick_params(axis='x', labelsize=fntsize)
        ax.tick_params(axis='y', labelsize=fntsize)
        plt.show()


class MultiTransients(object):
    """ list of transients corresponding to a certain dependence series"""

    datadir = 'C:/Users/sagustss/py_code/DATA'

    def __init__(self, transients_list=None, series_name=None, description=None, key_parameter=None):
        """
        Initialize the transients_list can be a list of transient objects or of path strings pointing to data files
        :param transients_list: list
            can be a list of transient objects or of path strings pointing to data files or None
            if not given, will create an empty object, where to later load data.
        :param series_name: str
            title of the series, used for printing on graphs for example.
        :param description: str
            some description.. not really used.
        :param key_parameter: str
            the parameter which changes throughout each scan, making it a "key_parameter" dependence series.
        """

        if transients_list is None:
            self.transients = []
            self.key_parameter = None
            self.description = None
            self.series_name = None
            self.material = None
            self.key_parameter_list = None
            self.metadata = None

        else:
            if type(transients_list[0]) == Transient:
                self.transients = transients_list
                self.import_metadata_from_transients()
            elif type(transients_list[0]) in (str, dir):
                self.transients = []
                self.import_files(transients_list)
            self.metadata = self.get_metadata()
            if key_parameter is None:
                self.key_parameter = self.transients[0].key_parameter
                if self.key_parameter is None:
                    self.key_parameter = self.get_dependence_parameter()
            else:
                self.key_parameter = key_parameter
            if description is None:
                self.description = self.transients[0].description
                if self.description is None:
                    self.description = None
            else:
                self.description = description
            if series_name is None:
                self.series_name = self.transients[0].series_name
                if self.series_name is None:
                    self.series_name = None
            else:
                self.series_name = series_name

            self.key_parameter_list = []

        self.sort_scan_list_by_parameter()
        self.update_key_parameter_list()

    # %% Metadata

    def give_name(self, name=None):

        if name is None:
            name = input('Choose a name for this series')
        self.series_name = name
        for scan in self.transients:
            self.series_name = name

    def get_dependence_parameter(self):  # todo: fix it, doesnt work, unhashable dict in definition of valdic
        """find the variable parameter within the series of scans.
        :returns : str name of dependence parameter, list of str if multiple dependence parameters found, also prints
        message.
        """
        dependence_parameter = []
        for key in self.metadata:
            valdic = {i: self.metadata[key].count(i) for i in frozenset(self.metadata[key])}
            if len(valdic) != 1:
                dependence_parameter.append(key)
        if 'date' in dependence_parameter:
            dependence_parameter.remove('date')
        if 'originalfilename' in dependence_parameter:
            dependence_parameter.remove('originalfilename')
        if len(dependence_parameter) > 1:
            print("Warning: multiple variables change between scans")
            print('Please choose between:')
            for parameter in dependence_parameter:
                print(parameter)
            dependence_parameter = input('\nEnter chosen key parameter: ')
        if len(dependence_parameter) < len(self.transients):
            print('Warning: multiple scans with same key parameter')
        else:
            dependence_parameter = str(dependence_parameter)
        return dependence_parameter

    def get_metadata(self):
        """Create a Dictionary of all metadata from all single scans.
        Each entry of the dictionary represents a parameter. Its values are a
        list of the value corresponding to the scan."""
        # get metadata from first scan, to initialize dictionary
        metadata = self.transients[0].get_metadata()

        for key in metadata:
            metadata[key] = [metadata[key]]
        # construct a dictionary containing all metadata
        skip = True  # skip first entry, since it was already written during
        # initialization
        for transient in self.transients:
            if not skip:
                md = transient.get_metadata()
                for key in metadata:
                    metadata[key].append(md[key])
            else:
                skip = False

        return metadata

    def update_transients_metadata(self):
        """ assign metadata from multitransient object to each scan"""
        for scan in self.transients:
            scan.key_parameter = self.key_parameter
            scan.description = self.description
            scan.series_name = self.series_name
            scan.material = self.material

    # %% Import Export

    def import_metadata_from_transients(self):

        metadata = self.get_metadata()
        try:
            self.key_parameter = metadata['key_parameter'][0]
        except KeyError:
            pass
        try:
            self.material = metadata['material'][0]
        except KeyError:
            pass
        try:
            self.description = metadata['description'][0]
        except KeyError:
            pass
        try:
            self.series_name = metadata['series_name'][0]
        except KeyError:
            pass

    def import_files(self, files, append=False, key_parameter=None, description=None):
        """imports any series of data files. Files can be:
               - string of full path of a single scan
               - list of full paths of a single scan
               - folder from which all files will be imported
           - append : if true, appends new scans to object, if false overwrites.
        """
        if not append:
            self.transients = []  # clear scans in memory
            # check if 'files' is single file (str), list of files ([str,str,...]) or folder containing files.
        if isinstance(files, str):
            self.transients.append(Transient(key_parameter=key_parameter, description=description))
            self.transients[-1].import_file(files)
            print('Imported file ' + files)
        elif isinstance(files, list) or isinstance(files, tuple):
            for i in range(len(files)):
                self.transients.append(Transient(key_parameter=key_parameter, description=description))
                self.transients[-1].import_file(files[i])
            print('Imported files form list')
        elif os.path.isdir(files):
            folderlist = os.listdir(files)
            for i in range(len(folderlist)):
                fullpath = files + '//' + folderlist[i]
                self.transients.append(Transient(key_parameter=key_parameter, description=description))
                self.transients[-1].import_file(fullpath)
                print('Imported files form folder')
        self.import_metadata_from_transients()
        # self.key_parameter = self.get_dependence_parameter()
        # self.sort_scan_list_by_parameter()  # todo: uncomment when get dependence parmater is fixed

        # %% data analysis

    def saveas_csv(self, directory=None):  # todo: implement dynamic paramter choosing option
        """ creates a directory inside the given directory where it will save all data in csv format."""
        if directory is None:
            directory = gfs.choose_folder('C:/Users/sagustss/py_code/DATA')
        save_dir = directory + '/' + self.series_name + '_' + self.key_parameter + '/'
        if os.path.exists(save_dir):
            n = 1
            new_save_dir = save_dir + '_1'
            while True:
                if os.path.exists(new_save_dir):
                    n += 1
                    new_save_dir = new_save_dir.split('_')[0] + '_' + str(n)
                else:
                    save_dir = new_save_dir
        os.makedirs(save_dir)
        self.update_transients_metadata()
        for item in self.transients:
            item.export_file_csv(save_dir)

    def clean_data_all_scans(self, cropTimeScale=True, shiftTime=0, flipTime=True, removeDC=True, filterLowPass=True,
                             flipTrace=False):
        """

        :return:
        """
        for transient in self.transients:
            transient.clean_data(cropTimeScale=cropTimeScale, shiftTime=shiftTime, flipTime=flipTime, removeDC=removeDC,
                                 filterLowPass=filterLowPass, flipTrace=flipTrace)

    def input_attribute(self, attribute_name, value):
        """
        manually input values for metadata attributes
        :param attribute_name: name of parameter or attribute
        :param value: value to assign to parameter
        """
        setattr(self, attribute_name, value)

    def sort_scan_list_by_parameter(self, reverse=False):
        transients_list = self.transients
        parameter = self.key_parameter
        # self.transients = sorted(self.transients, key=lambda transients: getattr(transients, parameter))
        self.transients.sort(key=lambda x: getattr(x, parameter), reverse=reverse)
        self.update_key_parameter_list()

        # return sorted_list

    def update_key_parameter_list(self):
        self.key_parameter_list = []
        for transient in self.transients:
            self.key_parameter_list.append(getattr(transient, self.key_parameter))

    # %% analysis

    def filter_low_pass(self, cutHigh=0.1, order=2):
        for item in self.transients:
            item = item.filter_low_pass(cutHigh, order)

    def remove_DC_offset(self):
        for item in self.transients:
            item = item.remove_DC_offset()

    def flip_time(self):
        for item in self.transients:
            item = item.flip_time()

    def flip_trace(self):
        for item in self.transients:
            item = item.flip_trace

    def shift_time(self, tshift):
        for item in self.transients:
            item = item.shift_time(tshift)


            # %% plot functions

    # %% plot

    def quickplot(self, figure=1):
        """ simple plot of a list of transients """  # todo: move to transients.py -> under multitransients()
        fig = plt.figure(num=figure)
        plt.clf()
        ax = fig.add_subplot(111)
        ax.set_xlabel('Time [ps]', fontsize=18)
        ax.set_ylabel('Differential Reflectivity', fontsize=18)
        ax.set_title(self.series_name, fontsize=26)
        ax.tick_params(axis='x', labelsize=12)
        ax.tick_params(axis='y', labelsize=12)

        # todo: make nice color iteration, that follows parameter value
        colorlist_length = len(self.transients)
        colorlist = plt.cm.rainbow(np.linspace(0, 1, colorlist_length))
        color = iter(colorlist)

        for curve in self.transients:
            xdata = curve.time
            ydata = curve.trace
            label = str(getattr(curve, self.key_parameter)) + str(curve.get_unit(self.key_parameter))
            col = next(color)
            ax.plot(xdata, ydata, c=col, label=label, alpha=0.7)
        return fig

    def quickplot_OLD(self):
        """ simple plot of a list of transients """
        fig = plt.figure(num=516542)
        plt.clf()  #
        ax = fig.add_subplot(111)
        ax.set_xlabel('Time [ps]', fontsize=18)
        ax.set_ylabel('Differential Reflectivity', fontsize=18)
        ax.set_title(self.series_name, fontsize=26)
        ax.tick_params(axis='x', labelsize=12)
        ax.tick_params(axis='y', labelsize=12)

        key_parameter_values = []
        for transient in self.transients:
            key_parameter_values.append(getattr(transient, self.key_parameter))

        key_parameter_max = max(key_parameter_values)
        print(key_parameter_max)
        n = 1
        # while key_parameter_range % 10 != 0:
        #     n *=10
        #     key_parameter_range * n

        colorlist = cm.rainbow(np.linspace(0, 1, 1001))
        # colorlist = cm.rainbow(np.logspace(0,3,1000)) / 100

        for curve in self.transients:
            # l = str(scn[i].temperature) + 'K'
            xdata = curve.time
            ydata = curve.trace
            parameter = float(getattr(curve, self.key_parameter))
            parameter_label = str('{0} {1}'.format(parameter, curve.get_unit(self.key_parameter)))
            color_number = (parameter / key_parameter_max) * 999
            print(color_number, parameter)
            ax.plot(xdata, ydata, c=colorlist[color_number], label=str() + 'K', alpha=0.5)
        plt.draw()
        return fig

    def rrPlot3d(self, Yparameter='Sample Orientation', title='3dplot', Xlabel='Time, ps',
                 Zlabel='Kerr rotation (mrad)',
                 colormap='viridis'):  # todo: correct to new TransientsSet() class system
        '''plot 3d graf with time on X trace on Z and selected parametr on Y '''
        # create 3 lists of X Y Z data
        time = []
        trace = []
        ypar = []
        # for every scan object takes values
        for item in self.transients:
            time.append(item.time)
            trace.append(item.trace)
            # on Y axis will be chosen parameter which exist in scan object
            ypar.append(item.parameters[Yparameter][0])
        # Make proper arrays from lists with data
        Ypar = []

        for item in range(len(self.transients[0].time)):
            Ypar.append(ypar)

        X = np.array(time)
        Y = np.transpose(np.array(Ypar))
        Z = np.array(trace)

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

    def fit_transients(self, fit_function, parameters, fit_from=0, fit_to=0, method='curve_fit', ext_plot=None,
                       print_results=True, recursive_optimization=False, colorlist=None):
        """
            Fit given model to a series of Transients.
        :param fit_function:
            Model which will be fitted to the data
        :param parameters: list, list of lists
            Initial parameters for the given function
        :param fit_from: int
            Minimum  from which to perform fit.
        :param fit_to: int
            Maximum data point (x axis) from which to perform fit.
        :param method: function
            Fitting method used: supports 'curve_fit'
        :param recursive_optimization: bool
            If true, it uses the optimized fit from previous cycle to initialize the next fitting
        :param ext_plot: bool
            if true, plots the results in a matplotlib figure
        :return all_popt: dict
            dictionary with transient label as key and fit optimized parameters as values
        :return all_pcov: dict

        """
        if ext_plot is None:
            fig = plt.figure('Fit of transients')
            plt.clf()
            ax = fig.add_subplot(111)
            ax.set_xlabel('Time [ps]', fontsize=18)
            ax.set_ylabel('Differential Reflectivity', fontsize=18)
            ax.set_title(self.series_name, fontsize=26)
            ax.tick_params(axis='x', labelsize=12)
            ax.tick_params(axis='y', labelsize=12)
        if colorlist is None:
            colorlist_length = len(self.transients)
            colorlist = plt.cm.rainbow(np.linspace(0, 1, colorlist_length))

        else:
            ax = ext_plot
        color = iter(colorlist)


        all_popt = {}
        all_pcov = {}
        last_popt = parameters
        for i, transient in enumerate(self.transients):
            xdata = transient.time[fit_from:-fit_to]
            ydata = transient.trace[fit_from:-fit_to]
            label = '{0} {1}'.format(transient.key_parameter_value, transient.get_unit(transient.key_parameter))

            try:
                if len(parameters[0]) > 1:
                    guess = parameters[i]
            except TypeError:
                if recursive_optimization:
                    guess=last_popt
                else:
                    guess = parameters

            all_popt[label] = []
            all_pcov[label] = []
            if method == 'curve_fit':
                try:
                    popt, pcov = curve_fit(fit_function, xdata, ydata, p0=guess)
                    if recursive_optimization:
                        last_popt = popt
                    if print_results:
                        print('{0}: popt: {1}'.format(label, popt))

                    col = next(color)
                    ax.plot(xdata, fit_function(xdata, *popt), '--', c=col)
                    ax.plot(xdata, ydata, c=col, label=label, alpha=0.5)
                    all_popt[label] = popt
                    all_pcov[label] = pcov

                except RuntimeError:
                    print('no fit parameters found for transient: {}'.format(label))
            elif method == 'fmin':
                print('fmin not yet implemented')  # todo: add support for fmin

        if ext_plot:
            pass
            # plt.show()
        return all_popt, all_pcov

class Data(object):
    """ This object stores data obtained from a fit such as decay times, amplitudes etc and provides analysis tools"""
    def __init__(self,x,y,x_label,y_label):
        """ initialization"""
        self.x_data = x
        self.y_data = y
        self.x_data = np.array(self.xdata)
        self.y_data = np.array(self.ydata)

        self.x_label = x_label
        self.y_label = y_label

    def quickplot(self, fntsize=15, title='Dependence', clear=False, plt_handle=None, show=True):
        """

        :param title: str
            window title
        :param clear:
            if true clears the figure before replotting
        :param plt_handle:
            handle to which to add the plot, used to plot in a pre-existant figure
        :param show:
            if true, uses plt.show()
        :return:

        """
        if plt_handle is None:
            fig = plt.figure(num=title)
            ax = fig.add_subplot(111)
        else:
            ax = plt_handle

        ax.scatter(self.x_data, self.y_data, 'o')
        ax.set_xlabel(self.x_label, fontsize=15)
        ax.set_ylabel(self.y_label, fontsize=15)
        ax.set_title(title, fontsize=15)
        ax.tick_params(axis='x', labelsize=15)
        ax.tick_params(axis='y', labelsize=15)
        if show:
            plt.show()

    def fit_data(self, fit_function, initial_parameters, x_range=(0, 0)):
        """

        :param function: func
            function used for data fitting
        :param initial_parameters: list
            initial parameters for inizializing fitting procedure
        :param x_range: [int,int]
            limits to use for fitting [ start , stop ]
        :param ydata: lsit or np.array
            y_data, takes precedence over self.y_data
        :return popt:
            optimized parameters
        :return pcov:
        """
        if x_range == (0,0):
            x_range = (0,len(self.x_data))
        x_data = self.x_data[x_range[0]:x_range[1]]
        y_data = self.y_data
        popt, pcov = curve_fit(fit_function, x_data, y_data, p0=initial_parameters)
        return popt, pcov


class FitFunction(object):
    """ wrapper class for fit functions"""
    def __init__(self):
        """ """
    @staticmethod()
    def doubleExp_lin_const(x, A, t0, c, d):
        return A * (1 - np.exp(- x / t0)) + c * x + d

    @staticmethod()
    def double_exponential_pos_neg(x, A1, t1, A2, t2, c, d):
        return A1 * (1 - np.exp(- x / t1)) - A2 * (1 - np.exp(- x / t2)) + c * x + d

    @staticmethod()
    def expfunc(x, A, t0, c):
        return A * np.exp(- x / t0) + c

if __name__ == "__main__":
    main()
