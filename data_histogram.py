# -*- coding: utf-8 -*-
"""
Created on Sat Oct 25 15:09:32 2014

@author: Philip
"""

import matplotlib.pyplot as plt
from pylab import *   #include numpy, plotting function
import numpy as np
from sub import excel
plt.close("all")


filePath = 'C:/Users/Philip/Google Drive/Python/Patch_analysis'
filename1 = 'patch.xlsx'
filename2 = 'Nopatch.xlsx'
file1=filePath+'/'+filename1
file2=filePath+'/'+filename2

patch=excel.openexcel(file1)
nopatch=excel.openexcel(file2)


fig, ax=plt.subplots()

max_diff=max(np.nanmax(patch), np.nanmax(nopatch))
hist_range=(-1*max_diff, max_diff)
n, bins, patches = ax.hist(patch, 20, range=hist_range, color='b', normed=True, alpha=0.5, label='patched')
n_shift, bins, patches = ax.hist(nopatch, 20, range=hist_range, color='r', normed=True, alpha=0.5, label='Nopatch')
ax.set_title('Histogram of patched and nopatched, $\Delta F / F$')
#ax.axis([max_diff*-1.2, max_diff*1.2, 0, n.max()*1.2])
ax.legend()