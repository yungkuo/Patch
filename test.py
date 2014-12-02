# -*- coding: utf-8 -*-
"""
Created on Sun Nov 30 11:41:23 2014

@author: Philip
"""

plt.figure()
abs_diff1=np.abs(diff_th1[:,2])
abs_diff2=np.abs(diff_th2[:,2])
max_diff=max(np.nanmax(abs_diff1), np.nanmax(abs_diff2))
hist_range=(-1*max_diff, max_diff)
n, bins, patches = plt.hist(diff_th1[:,2], binNum, range=hist_range, color='b', normed=True, alpha=0.5, label='1st-2nd')
n_shift, bins, patches = plt.hist(diff_th2[:,2], binNum, range=hist_range, color='r', normed=True, alpha=0.5, label='2nd-1st')
#plt.set_title('Histogram of $\Delta I/I_{Max}$ , thresholded')
plt.axis([max_diff*-1.2, max_diff*1.2, 0, n.max()*1.2])

plt.figure()
abs_diff1=np.abs(Ndiff1hist[1][7])
abs_diff2=np.abs(Ndiff2hist[1][7])
max_diff=max(np.nanmax(abs_diff1), np.nanmax(abs_diff2))
hist_range=(-1*max_diff, max_diff)
n, bins, patches = plt.hist(Ndiff1hist[1][7], binNum, range=hist_range, color='b', normed=True, alpha=0.5, label='1st-2nd')
n_shift, bins, patches = plt.hist(Ndiff2hist[1][7], binNum, range=hist_range, color='r', normed=True, alpha=0.5, label='2nd-1st')
#plt.set_title('Histogram of $\Delta I/I_{Max}$ , thresholded')
plt.axis([max_diff*-1.2, max_diff*1.2, 0, n.max()*1.2])
