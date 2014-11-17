# -*- coding: utf-8 -*-
"""
Created on Sun Nov 16 17:40:38 2014

@author: Philip
"""
burstmask = np.zeros(len(barray[0]), dtype=bool)
bursts = burstdata[0]
for start, stop, score in bursts:
    burstmask[start : stop] = True
