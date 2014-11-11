# -*- coding: utf-8 -*-
"""
Created on Sun Nov 09 10:56:49 2014

@author: Philip
"""
b=diff_th1[:,0]
a=[b[i] for i in range(len(b)) if isnan(b[i]) != 1]

for i in range(10):
    lista=len(a)-i
    