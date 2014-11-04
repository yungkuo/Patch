# -*- coding: utf-8 -*-
"""
Created on Mon Sep 01 21:57:55 2014

@author: Philip
"""

def ReadBinMov(fileName, nrow, ncol):
    import numpy as np
    
    with open(fileName, "rb") as fid:
        data_array = np.fromfile(fid, np.uint16)   # read little endian format
    
    # reshape vector into appropriately oriented, 3D array
    frame = len(data_array)/(nrow*ncol);
    mov = data_array.reshape((frame, nrow, ncol))  
    # arrange row first 
    # For MATLAB, arrange col first 
    return mov, frame