# -*- coding: utf-8 -*-
"""
Created on Mon May 15 09:57:29 2017

@author: S.Y. Agustsson
"""

from functionlibrary import genericfunctions as gfs
import numpy as np
import scipy as sp
import scipy.io as spio
import scipy.signal as spsignal
import matplotlib.pyplot as plt
import os
import re
import pandas


def main():
    testmat = 'RuCl3-Pr-0.5mW-Pu-1.5mW-T-005.0k-1kAVG.mat'
    testcsv = 'RuCl3- 2017-04-19 17.33.14 Pump1.5mW Temp7.0K.txt'
    testpath = '..//test_scripts//test_data//'
    savepath = "E://DATA//RuCl3//"

    matfile = testpath + testmat
    csvfile = testpath + testcsv

    scanMat = Transient()
    scanCSV = Transient()
    scanMat.import_file(matfile)

    print('imported mat file\n')
    # fig = plt.figure(num=1351)
    # axis = fig.add_subplot(111)
    # axis.plot(scanMat.time, 'o')
    # print(scanMat.time[0:20])
    # print(scanMat.time[-20:-1])
    scanMat.quickplot(raw=False)
    scanMat.export_file_csv(testpath)
    # getfile = gfs.choose_filename(testpath)
    scanCSV.import_file_csv(testpath + scanMat.name + '.txt')
    print('result: ')
    print(scanCSV.time[0:20])
    scanCSV.quickplot()


#    scan2.import_single_file(csvfile)
#    scan2.quickplot


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

    def __init__(self):
        """ """

        ######################################
        #                Data                #
        ######################################

        self.raw_time = np.array([])  # time data
        self.raw_trace = np.array([])  # raw trace data

        self.time = np.array([])  # cleaned time axis
        self.trace = np.array([])  # cleaned and modified data trace
        # ignore list for metadata export. Add here any further non-metadata attributes created in this class.
        self.DATA_ATTRIBUTES = ('raw_time', 'raw_trace', 'time', 'trace', 'DATA_ATTRIBUTES')

        ######################################
        #              Metadata              #
        ######################################

        self.name = ''  # String used as file name for saving data
        self.material = ''  # Material name
        self.date = ''  # Scan date in format YYYY-MM-DD hh.mm.ss
        self.original_filepath = ''  # Path to original raw file

        # parameters

        self.pump_power = 0  # Pump Power [mW]
        self.probe_power = 0  # Probe Power [mW]
        self.destruction_power = 0  # Destruction Power [mW]

        # Spot size represents the FWHM diameter from Gaussian fit of the beam profile
        self.pump_spot = 0  # Pump beam spot size [micrometers]
        self.probe_spot = 0  # Probe beam spot size [micrometers]
        self.destruction_spot = 0  # Destruction beam spot size [micrometers]

        # Excitation densities calculated from power, spot size and repetition rate
        self.pump_energy = 0
        self.probe_energy = 0
        self.destruction_energy = 0

        # destruction pulse parameters:
        self.destruction_delay = 0  # delay between pump and destruction pulses, [ps]

        # Polarization are measured clockwise in propagation direction of the beam, 0 = 12o'clock
        self.pump_polarization = 0  # Pump beam polarization [deg]
        self.probe_polarization = 0  # Probe beam polarization [deg]
        self.destruction_polarization = 0  # Destruction beam polariz. [deg]
        self.sample_orientation = 0  # Sample orientation [deg]
        self.temperature = 0  # Temperature [K]
        self.R0 = 0  # Static reflectivity
        # self.T0 = 0             # Static Transimittivity

        ######################################
        #              Analysis              #
        ######################################

        self.analysis_log = {}  # Keeps track of analysis changes performed

    # %% metadata management

    def calc_energy_densities(self, rep_rate=283000):
        """ recalculate metadata depending on given parameters
        calculates:
               - energy densities
        """

        beams = ['pump', 'probe', 'destruction']

        for beam in beams:
            power = getattr(self, (beam + '_power'))
            spot = getattr(self, (beam + '_spot'))
            if beam == 'pump':
                rep_rate = rep_rate / 2  # pump has half reprate
                # (darkcontrol)
            energy = gfs.get_energy_density(spot, power, rep_rate)
            setattr(self, (beam + '_energy'), energy)

    def init_metadata(self):  # todo: Write this method
        """
        Initializes all metadata variables obtainable from the file name, when a file has the correct naming structure:
            MaTeRiAl_pu12mW_pr5mW_de50mW_temp4K_pupol45_prpol-45_001.xxx
        """
        print("method initMetadata still not made...")

    def get_metadata(self):
        """ Returns a dictionary containing all metadata information available
            keys are attribute names and values the corresponding value.
            :rtype: dict
        """
        metadata = {'analysis_log': {}}  # make dict for metadata. Also create analysis_log entry, used later
        var = self.__dict__

        for key in var:
            if key not in self.DATA_ATTRIBUTES:  # ignore time, trace and all other fields defined in data_attributes
                if 'analysis_log' in key:
                    metadata['analysis_log'][key[13:-1:]] = var[key]  # todo: check why it needs key[13:-1]
                else:
                    try:
                        # if parameter is a number != 0 append to metadata
                        if float(var[key]) != 0:
                            metadata[key] = var[key]
                    # if not a number, and not an empty string
                    except ValueError:
                        if len(var[key]) == 0:
                            pass
                        else:
                            metadata[key] = var[key]
                    # if not a number, and not an empty list
                    except TypeError:
                        if len(var[key]) == 0:
                            pass
                        else:
                            metadata[key] = var[key]

        return metadata

    def log_it(self, key, overwrite=False, *args, **kargs):
        """
        Generate log entry for analysis_log.
            creates a key with given key in analysis_log, making it:
                - boolean if no other args or kargs are given, flips previous
                values written in log
                - list if *args are passed
                - dictionary if **kargs are passed
            if overwrite is False, it appends values on previous logs,
            if True, it obviously overwrites them.
            :type overwrite: bool
        """

        key_string = 'analysis_log[' + key + ']'
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
            prev = self.analysis_log[key]
            # if it finds such a key in this library, it overwrites/appends new information.
            if entry == None:  # trigger boolean behaviour, flipping previous registered status if available
                self.analysis_log[key] = not prev

            elif entry is list:
                if overwrite:
                    self.analysis_log[key] = entry
                    setattr(self, key_string, prev + entry)
                else:
                    entry = prev.append(entry)
                    self.analysis_log[key] = entry
            elif entry is dict:
                if overwrite:
                    self.analysis_log[key] = entry
                else:
                    new_entry = {}
                    for key in entry:
                        if key in prev:
                            new_entry[key] = prev[key].append(entry[key])
                        else:
                            new_entry[key] = entry[key]

                    self.analysis_log[key] = new_entry
        except KeyError:  # rises Key error when key was not previously assigned
            if entry == None:  # if boolean, its it's first usage, therefore set it to true
                setattr(self, key_string, True)
            else:  # otherwise just append the list/dictionary.
                setattr(self, key_string, entry)

    def make_name(self):
        """Define name attribute as material_date."""
        self.name = str(self.material) + '_' + str(self.date)

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

    def import_file(self, filepath, cleanData=True):
        """Imports a file, csv or .mat
        uses self.import_file_mat() and self.import_file_csv,
        depending on file extension"""
        try:  # if it finds the file requested...
            ext = os.path.splitext(filepath)[-1].lower()
            if ext == '.mat':
                if os.path.basename(filepath).lower() != 't-cal.mat':
                    self.import_file_mat(filepath)
                else:
                    print('Ignored t-cal.mat')
            elif ext == '.txt':
                self.import_file_csv(filepath)
            else:
                print("Invalid format. Couldn't import file: " + filepath)
            # self.importMetadata()
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
                    value_string = word[1]
                    if self.get_unit(key):
                        # if parameter expects units, get only numbers,
                        value = float(re.findall("\d+\.\d+", value_string)[0])
                    else:  # otherwise get the whole string
                        value = value_string
                    # create/assign attribute from imported parameter
                    setattr(self, key, value)

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

        self.make_name()  # creates a name for the file
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
            if self.analysis_log['Flip Time']: pass
            shift = np.average(self.trace[0:window:1])
        except KeyError:
            tpoints = len(self.time)
            shift = np.average(self.trace[tpoints - window:tpoints:1])
        self.trace = self.trace - shift
        self.log_it('Remove DC', window=window, shift=shift)

    def filter_low_pass(self, cutHigh=0.1, order=1, return_frequency=False):
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
        self.log_it('Low Pass Filter',
                    frequency=frequency,
                    nyq_factor=cutHigh,
                    order=order)
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


class TransientsSet(object):
    def __init__(self):
        """ """
        self.transients_list = []
        self.transients = []

    # %% Metadata

    def get_dependence_parameter(self):
        """find the variable parameter within the series of scans.
        :returns : str name of dependence parameter, list of str if multiple dependence parameters found, also prints
        message.
        """
        metadata = self.get_metadata()
        dependence_parameter = []
        for key in metadata:
            valdic = {i: metadata[key].count(i) for i in metadata[key]}
            if len(valdic) != 1:
                dependence_parameter.append(key)
        if 'date' in dependence_parameter:
            dependence_parameter.remove('date')
        if 'originalfilename' in dependence_parameter:
            dependence_parameter.remove('originalfilename')
        if len(dependence_parameter) > 1:
            print("Warning: multiple variables change between scans")
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
        for scan in self.transients:
            if not skip:
                md = scan.get_metadata()
                for key in metadata:
                    metadata[key].append(md[key])
            else:
                skip = False
        return metadata

    # %% Import Export

    def import_files(self, files, append=False):
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
            self.transients.append(Transient())
            self.transients[-1].import_single_file(files)
            print('Imported file ' + files)
        elif isinstance(files, list):
            for i in range(len(files)):
                self.transients.append(Transient())
                self.transients[-1].import_single_file(files[i])
            print('Imported files form list')
        elif os.path.isdir(files):
            folderlist = os.listdir(files)
            for i in range(len(folderlist)):
                fullpath = files + '//' + folderlist[i]
                self.transients.append(Transient())
                self.transients[-1].import_single_file(fullpath)
                print('Imported files form folder')
        self.update_scan_list()

    def update_scan_list(self):
        """ Update the list of names of the single scans in self.scans.
        Name is taken from 'name' attribute of Transient() class.
        """
        self.transients_list = []
        for i in range(len(self.transients)):
            self.transients_list.append(self.transients[i].name)

            # %% data analysis

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


class TransientFit(object):
    """ class for making fits to transient data"""

    def __init__(self):
        """ initialize"""
        pass

    def make_fit_function(self):
        """ create the fitting function to use.
        """


if __name__ == "__main__":
    main()
