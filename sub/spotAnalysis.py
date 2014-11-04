# -*- coding: utf-8 -*-
"""
Created on Wed Sep 24 14:02:08 2014

@author: KyoungWon
"""
import numpy as np
from pylab import * 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, leastsq
from scipy import asarray as ar,exp
import matplotlib.gridspec as gridspec



def avg_period(t, signal , period, threshold):
    result=np.zeros((period))
    ncycle=len(signal)/period
    frame=len(signal)
    thresh_line=np.ones(frame)*threshold    
    F=np.zeros((2))
    for i in range(len(signal)):
        k = int(mod(i,period))
        result[k] = result[k] + signal[i]
    F[0]=sum(result[:period/2])/(period/2)
    F[1]=sum(result[period/2:])/(period/2)
    dFF=(F[1]-F[0])/F[0]*100
    sortedI=np.zeros(len(signal))
    for i in range(ncycle):
        sortedI[i*(period/2) : (i+1)*(period/2)] =signal[i*period:i*period+period/2]
    for j in range(ncycle):
        sortedI[frame/2 + j *(period/2) : frame/2 + (j+1)*(period/2)] = signal[j*period+period/2:(j+1)*period]
    return result, dFF, np.array(sortedI)

    
def threshold_avgperiod(threshold, signal, period):
    nbin=np.zeros((period))
    result=np.zeros((period))
    F=np.zeros((2))
    for i in range(len(signal)):
        if signal[i] > threshold:
            k=int(mod(i,period))
            result[k]=result[k]+signal[i]
            nbin[k]=nbin[k]+1
    norm_result=result/nbin
    F[0]=sum(result[:period/2])/sum(nbin[:period/2])
    F[1]=sum(result[period/2:])/sum(nbin[period/2:])
    dFF=(F[1]-F[0])/F[0]*100
    return norm_result, dFF

def difference(curve, period):
    nperiod=len(curve)/period
    diff=np.zeros(nperiod)
    maxval=curve.max()
    for i in range(nperiod):
        phase1=sum(curve[period*i:period*i+period/2])
        phase2=sum(curve[period*i+period/2:(i+1)*period])
        diff[i]=phase1-phase2
    diff=diff/maxval
    diff=diff/(period/2)   #normalization
    return diff
        
        
def filted_diff(curve, period, threshold):
    nframe=len(curve)
    ncycle=nframe/period
    diff1=np.ones(ncycle)
    diff1[:]=np.NAN
    diff2=np.ones(ncycle)
    diff2[:]=np.NAN
    #Von=np.ones(ncycle)
    #Von[:]=np.NAN
    #Voff=np.ones(ncycle)
    #Voff[:]=np.NAN
    k=0
    l=0
    F=[]
    for i in range(ncycle):
        if threshold <= min(curve[i*period : (i+1)*period]):
            diff1[k]=sum(curve[i*period : i*period+period/2])-sum(curve[i*period+period/2 : (i+1)*period])
            F.append(curve[i*period:(i+1)*period])
            #Von[k]=sum(curve[i*period : i*period+period/2])
            #Voff[k]=sum(curve[i*period+period/2 : (i+1)*period])
            k=k+1
    for i in range(ncycle-1):
        if threshold <= min(curve[i*period+(period/2) : (i+1)*period+(period/2)]):
            diff2[l]=sum(curve[i*period+(period/2) : (i+1)*period])-sum(curve[(i+1)*period : (i+1)*period+(period/2)])
            #Von[k]=sum(curve[i*period : i*period+period/2])
            #Voff[k]=sum(curve[i*period+period/2 : (i+1)*period])
            l=l+1
    F=np.array(F)
    Favg=np.mean(F)
    dff1=diff1/(period/2)/Favg
    dff2=diff2/(period/2)/Favg
    dff_avg=(np.nanmean(dff1)-np.nanmean(dff2))/2
    return dff1, dff2, dff_avg #Von/(period/2)/np.nanmax(curve), Voff/(period/2)/np.nanmax(curve)


def evenodd(diff):       
    npoint=len(diff[0,:])
    binsize=30
    n=np.zeros((binsize, npoint))
    bins=np.zeros((binsize+1, npoint))
    even=np.zeros((binsize, npoint))
    odd=np.zeros((binsize, npoint))
    for i in range(npoint):
        abs_diff=abs(diff[:,i])
        hist_range=(-1*np.nanmax(abs_diff), np.nanmax(abs_diff))
        n[:,i], bins[:,i]  = np.histogram(diff[:,i], binsize, range=hist_range)
        even[:,i]=(n[:,i]+n[::-1,i])/2   # make an even function
        odd[:,i]=(n[:,i]-n[::-1,i])/2    # make an odd function
    return even, odd, n, bins
    
    
    
def gaussian(x, a, mu, sigma):
    return a*exp(-(x-mu)**2/(2*sigma**2))
def oddGaussian(x, a, mu, sigma):
    return a*exp(-(x-mu)**2/(2*sigma**2))*(x-mu)
    
