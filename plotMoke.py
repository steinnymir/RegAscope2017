# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from scipy import signal
from scipy import io
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm

#MData = scipy.io.loadmat('Mn2Au-MA70-20mW-10mW-RT-10kAVG')
#Data = MData['Daten']

def fitfunction(t,taupulse,taudecay,linearconst,A,C,t_zero):
    return  A*(-np.exp(taupulse**2/(8*taudecay**2))*(np.exp(-t/taudecay))+linearconst*t + C)*(1+scipy.special.erf(((2**0.5)*(t+t_zero-(taupulse**2/(4*taudecay))))/taupulse))

def import_redred(Filename):
    MData = io.loadmat(Filename)
    Data = MData['Daten']
    d=Data[0]
    shift=np.average(d[7460:7500:1])
    trace=30*(-Data[0]+shift)/57
    #trace=-trace  #/trace.min(axis=0)
    Time = -Data[2] 
     #Data[2,np.argmin(-Data[0])]
    b, a = signal.butter(2, 0.1, 'low', analog= False)
    filtered_trace=signal.lfilter(b,a,trace)
    return(Time,filtered_trace)
    
    
def fitredred(trace1):
    trace=trace1[0]
    Time=trace1[1]
    guess=(0.16,1.68,0.0011,0.4,-0.27,0.2)
    plt.plot(Time,trace)
    popt, pcov = scipy.optimize.curve_fit(fitfunction,Time,trace,p0=guess)
    #matplotlib.pyplot.plot(Time,fitfunction(Time, *popt))
    #matplotlib.pyplot.plot(trace[1],fitfunction(trace[1], 0.16,1.63,0.0011,0.4,-0.27,0.2))
    
    #matplotlib.pyplot.plot(Time,trace - fitfunction(Time,  *popt))
    return(popt)
    #matplotlib.pyplot.plot(np.fft.fft(trace-fitfunction(-Data[2], *popt)))
    #matplotlib.pyplot.plot(Time,fitfunction1(Time, *popt))
    #matplotlib.pyplot.plot(Time,fitfunction2(Time, *popt))
    #matplotlib.pyplot.plot(Time,fitfunction1(Time, *popt)*fitfunction2(Time, *popt))

def fftredred(trace, Time):
    freq=scipy.fftpack.fftfreq(trace.size,np.abs((Time[0]-Time[:1])/Time.size))
    FFT = abs(scipy.fft(trace))
    plt.plot(freq, FFT)
    
Names=[]
angle=np.arange(0,370,10)
A=[]
Trace=[]
Time=[]
#Names.append('Mn2Au-MA125-20mW-10mW-RT-5kAV-MOKE-0degree')
#angle.append(0)
i=0
while i<37 :
    Names.append('E:\\2017-02-24\\Mn2Au-MA125-30mWhor-20mWvert-RT-1kAV-MOKE-samplerot-'+str(0+i*10)+'degree')
    i+=1
PS=[]
MData = io.loadmat(Names[1])
Data = MData['Daten']
timeShift=Data[2,np.argmin(-Data[0])]
for item in Names:
    data=import_redred(item)
    dt=data[1]
    PS.append(np.average(dt[5640:5692:1]))
    trace=data[1]
    time=data[0]+timeShift
    Trace.append(trace[0:5777:1])
    Time.append(time[0:5777:1])
    #plt.plot(time,trace)
    #fftredred(trace, time)
    #tr=import_redred(item)
    #q=tr[0]
    #A.append(np.sum(q[5650:5700]))
    #matplotlib.pyplot.plot(angle,A)
X=np.array(Time)
Angle=[]
i=np.arange(0,5777,1)
for item in i:
    Angle.append(angle)
Y=np.transpose(np.array(Angle))
Z=np.array(Trace)
fig = plt.figure(num=1)
ax = fig.add_subplot(111, projection='3d')
# Grab some test data.
#X, Y, Z = axes3d.get_test_data(0.05)

# Plot a basic wireframe.
ax.plot_surface(X, Y, Z, rstride=1, cstride=100, cmap='viridis')
ax.set_xlabel('Time(ps)', fontsize=20, labelpad=20)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
ax.tick_params(axis='z', labelsize=20)

ax.set_ylabel('Angle (degrees)', fontsize=20, labelpad=20)
ax.set_zlabel('Kerr rotation (mrad)', fontsize=20, labelpad=20)
ax.set_title(' vertical Probe, horizontal Pump', fontsize=40)
#plt.show()

fig2=plt.figure(num=2)
ax1 = fig2.add_subplot(111)
ax1.plot(angle,PS)
ax1.set_xlabel('Angle (degrees)', fontsize=40)
ax1.set_ylabel('Kerr rotation (mrad)', fontsize=40)
ax1.set_title('Angular dependence of Kerr rotation at 2 ps \n (horizontal Probe, vertical Pump)',fontsize=40)
ax1.tick_params(axis='x', labelsize=20)
ax1.tick_params(axis='y', labelsize=20)


iteract=np.arange(0,18,1)
ass=0
assum=[]
for item in iteract:
    ass=Trace[item]+Trace[item+18]
    assum.append(ass)



crossect=[]
idx=0

for item in i:
    kr=[]
    trrr=[]
    for item2 in assum:
        tr=item2
        kr.append(tr[item])
    av=np.average(kr)
    zerro=np.max(kr-av)
    krnorm=(kr-av)/zerro 
    crossect.append(krnorm)
x=Time[0]
y=angle[0:18:1]
X1, Y1=np.meshgrid(x,y)
Z1=np.transpose(np.array(crossect))
fig4 = plt.figure(num=3)
ax3 = fig4.add_subplot(111)
ax4 = fig4.add_subplot(211)
ax4.plot(Time[8],Trace[8])
# Grab some test data.
#X, Y, Z = axes3d.get_test_data(0.05)

# Plot a basic wireframe.
ax3.pcolor(X1, Y1, Z1, cmap='pink')
#ax3.set_xlabel('Time(ps)', fontsize=20, labelpad=20)
#ax3.tick_params(axis='x', labelsize=20)
#ax3.tick_params(axis='y', labelsize=20)
#ax3.tick_params(axis='z', labelsize=20)


#ax3.set_ylabel('Angle (degrees)', fontsize=20, labelpad=20)
#ax3.set_zlabel('Kerr rotation (mrad)', fontsize=20, labelpad=20)
#ax3.set_title(' vertical Probe, horizontal Pump', fontsize=40)
plt.show()    


    
        