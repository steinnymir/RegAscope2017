# -*- coding: utf-8 -*-
"""

@author: Steinn Ymir Agustsson
"""

from lib import transient
from imp import reload
reload(transient)

import os
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.colorbar as colorbar
import matplotlib.cm as cmx
from mpl_toolkits.axes_grid1 import make_axes_locatable
%matplotlib inline

import pandas as pd


###############
#   import    #
###############



filepath = 'E:/data/RuCl3/mat/kerr_rotation/fluence_3.8K/'
files = os.listdir(filepath)

cutFreq = 0.01 # THz
t0 = 85
key_parameter = 'pump_power'
description = ''

usebg = False
trs = []
k_parameters = []
data = {}
if 'background.mat' in files:
    bg_trace = transient.Transient()
    bg_trace.import_file(filepath + 'background.mat',
                           cleanData=False,
                           key_parameter=key_parameter,description=description)
    files = files.remove('background.mat')

    usebg = True

for file in files:
    try:

        tr = transient.Transient()
        tr.import_file(filepath + file,
                           cleanData=False,
                           key_parameter=key_parameter,description=description)

        if usebg:
            tr.trace = tr.trace - bg_trace.trace
        tr.crop_time_scale()
        tr.shift_time(t0)
        tr.filter_low_pass(cutFreq)
    #     tr.flip_trace()
        tr.remove_DC_offset()
        tr.flip_time()

        data[getattr(tr,key_parameter)] = tr.trace


        k_parameters.append(float(getattr(tr,key_parameter)))

#        tr.trace = np.divide(tr.trace,max(tr.trace))
    #     tr.trace = np.divide(tr.trace,tr.pump_power)

        trs.append(tr)
    except Exception as exc:
        print('skipped file: {0}\nerror: {1}'.format(file,exc))
data = pd.DataFrame(data,index=tr.time)
print('Imported {0} scan(s) as {1} dependence'.format(len(trs),trs[0].key_parameter))


#################
#   plotting    #
#################

ax = []
fig = plt.figure(figsize=(16,9))
ax = fig.add_subplot(111)
jet = cm = plt.get_cmap('jet')
cNormLog = colors.LogNorm(vmin=min(k_parameters), vmax=max(k_parameters))
cNorm  = colors.Normalize(vmin=min(k_parameters), vmax=max(k_parameters))
scalarMap = cmx.ScalarMappable(norm=cNormLog, cmap=jet)
lines = []
for trace in trs:
    colorVal = scalarMap.to_rgba(trace.pump_power)
    retLine, = ax.plot(trace.time, trace.trace, color=colorVal)
    lines.append(retLine)
# make labels
    labels = []
for i in range(len(lines)):
    labels.append((lines[i],k_parameters[i]))
labels = sorted(labels, key=lambda kpar: float(kpar[1]))
lbLines=[]
lbVals=[]
for lbl in labels:
    lbLines.append(lbl[0])
    lbVals.append(lbl[1])
ax.legend(lbLines, lbVals, loc='upper right')
# plt.xlim((-5,5))
# plt.ylim((-0.0002,0.0008))
ax.grid()

























def main():
    pass

if __name__ == '__main__':
    main()