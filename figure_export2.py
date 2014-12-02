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

curve=spot_intensity[:,0]
ncycle=nframes // period
half=period/2
record1=0
record2=0
phase1up=[]
phase2up=[]
for i in range(ncycle):
    min1= np.min(curve[i*period:i*period+half])
    max1= np.max(curve[i*period:i*period+half])
    min2= np.min(curve[i*period+half:(i+1)*period])
    max2= np.max(curve[i*period+half:(i+1)*period])
    if min1 >= max2: #phase1up
        record1=record1+1
        record2=0
    elif max1 <= min2: #phase2up
        record2=record2+1
        record1=0
    else:
        if record1 > 0:
            t1=(record1, i)
            phase1up.append(t1)
        elif record2 > 0:
            t2=(record2, i)
            phase2up.append(t2)
        record1=0
        record2=0

Maxcorr1=sorted(phase1up, reverse=True)[0]
Maxcorr2=sorted(phase2up, reverse=True)[0]
Phase1upArray= curve[(Maxcorr1[1]-Maxcorr1[0]-1) *period : (Maxcorr1[1]+1)*period]
Phase2upArray= curve[(Maxcorr2[1]-Maxcorr2[0]-1) *period : (Maxcorr2[1]+1)*period]
length = len(Phase1upArray)/period
plt.figure()
plt.plot(curve[(Maxcorr1[1]-Maxcorr1[0]-1) *period : (Maxcorr1[1]+1)*period])
for i in range(1,length):
        plt.axvspan( i*period -0.5 , (i+0.5) *period -0.5, facecolor='g', alpha=0.4)