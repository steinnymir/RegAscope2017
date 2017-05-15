# -*- coding: utf-8 -*-
"""
Created on Mon May 15 10:23:20 2017

@author: S.Y. Agustsson
"""

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#                  old rubbish functions
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





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
        stime = str(time[i])\
        sraw = str(raw[i])
        sfilt = str(filt[i])
        file.write(stime + ',' + sraw + ',' + sfilt + '\n')
    print('file ' + filename + ' is ready')
    file.close()


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