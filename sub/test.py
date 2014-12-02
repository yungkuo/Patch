# -*- coding: utf-8 -*-
"""
Created on Sun Nov 16 17:40:38 2014

@author: Philip
"""
plt.close('all')
m = 10
threshold = -4

burstdata=[]
for i in range(npoint):
    #signal = bsearch_py.burstarray(c)
    signal = barray[i]
    a = bsearch_py.burstsearch_py(signal, m, threshold)
    burstdata.append(a)

tsignal = np.arange(signal.size)*2 + 1

ref_alternation = np.zeros(meanI.shape[0])
ref_alternation[::2] = 1

for i in range(npoint):
    plt.figure()
    plt.plot(meanI[:,i], '.-')
    signal = bsearch_py.burstarray(meanI[:,i])
    #plt.plot(tsignal, signal, '.-m')
    plt.axhline(threshold)

    ax = plt.twinx()
    ax.plot(ref_alternation, 'k', alpha=0.5)

    bursts=burstdata[i]
    for start, stop, score in bursts:
        plt.axvspan(start*period/2 - 0.5 , stop*period/2 + 0.5, facecolor='g', alpha=0.4)



