# -*- coding: utf-8 -*-
"""
Created on Sun Sep 21 10:31:06 2014

@author: Philip
"""
import numpy as np
import matplotlib.pyplot as plt

def movingaverage(vector, span):
    
    if span==2:
        print("Your moving average span is 2" )
        print("Moving average should be larger than 3")
        print("Adjust to 3 span as a moving average ")
        
    else:
        if span%2==0:
            print("Your moving average span is %d, which is even" % span)
            span=span-1
            print("Change span to %d" %span)
            
    hf_span=(span-1)/2
    length=len(vector)
    mvavg=np.zeros((length))
    mvavg[0]=vector[0]
    mvavg[-1]=vector[-1]
    for i in range(1,length-1):
        if i < hf_span:
            mvavg[i]=np.sum(vector[0:2*i+1])/(2*i+1)
        elif i >= length-hf_span:
            r_half=length-i-1
            mvavg[i]=np.sum(vector[-1-2*r_half:])/(r_half*2+1)
        else:
            mvavg[i]=np.sum(vector[i-hf_span:i+hf_span+1])/(hf_span*2+1)    

    return mvavg
    