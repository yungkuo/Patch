# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 20:04:19 2014

@author: KyoungWon
"""

import matplotlib.pyplot as plt
import numpy as np
from pylab import *

def DAQread(fileName, SqFrameRate, SqDAQrate, SqFreqV, frame):
    period=int(SqDAQrate/SqFreqV)
    tm=[]
    Vm=[]
    Im=[]
    with open(fileName, 'r') as f:
        next(f) #skip headings
        for line in f:
            fields= line.split('\t')
            tm.append(fields[0])
            Im.append(fields[1])
            Vm.append(fields[2])
    Vm=float16(np.array(Vm))/10
    Im=float16(np.array(Im))
    tm=float16(np.array(tm))
    downsamp=round((SqDAQrate/SqFrameRate))
    V=np.reshape(Vm, (len(Vm)/downsamp, downsamp))
    I=np.reshape(Im, (len(Im)/downsamp, downsamp))
    V=mean(V, axis=1)
    I=mean(I, axis=1)
    t=float64(arange(frame)/SqFrameRate)
    
    figure(3)
    plt.subplot(221)
    line_Current, = plt.plot(tm, Im/30, 'c', label="Current / (30 pA)" )
    line_Voltage, = plt.plot(tm, Vm, 'm', label=" Voltage (V)")
    plt.xlabel('time [s]')
    
    #
    #plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
               
    plt.subplot(222)
    line_Current, = plt.plot(tm[:period*5], Im[:period*5]/30, 'c', label="Current / (30 pA)" )
    line_Voltage, = plt.plot(tm[:period*5], Vm[:period*5], 'm', label=" Voltage (V)")
    plt.legend(bbox_to_anchor=(0.1, 1.2), loc=2, borderaxespad=0.)
    plt.xlabel('time [s]')
    
    return t, V, I