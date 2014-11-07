# -*- coding: utf-8 -*-
"""
Created on Mon Nov 03 16:01:40 2014

@author: KyoungWon
"""
from scipy.ndimage.filters import gaussian_filter1d
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


LPF_G1d=np.zeros_like(spot_intensity) 
HPF=np.zeros_like(spot_intensity) 
for i in range((npoint)):
    LPF_G1d[:,i]=gaussian_filter1d(spot_intensity[:,i], sigma=10)
    HPF[:,i]=spot_intensity[:,i]-LPF_G1d[:,i]
    
    
    fig= plt.figure()    
    ax=plt.subplot(111)
    plt.plot(spot_intensity[:,i], '-o' 'c', label="Original" )
    plt.plot(LPF_G1d[:,i], 'm', label="LPF filtered" )
    plt.plot(HPF[:,i], 'g', label="HP component" )
    ax.xaxis.set_major_locator(MultipleLocator(4))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))

    ax.xaxis.grid(True,'major',linewidth=2)
    ax.yaxis.grid(True,'major',linewidth=2)
    ax.grid(True)
   # for j in range(len(spot_zip(A,B):
   #     ax.annotate('%s)' %j, xy=(i,j), xytext=(30,0), textcoords='offset points')
   #     ax.annotate('(%s,' %i, xy=(i,j))

    
    plt.show()
