# -*- coding: utf-8 -*-
"""
Created on Tue Oct 14 13:15:32 2014

@author: KyoungWon
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 23:19:14 2014

@author: Philip
"""
from numba import jit
from IPython import get_ipython
ip = get_ipython()
ip.magic("reset -f")
import matplotlib.pyplot as plt
from matplotlib.widgets import CheckButtons
from pylab import *   #include numpy, plotting function
import numpy as np
from sub import ReadBinMov, polygon, smooth, point, DAQ, spotAnalysis, multiplot, fit
plt.close("all")



"""
Control Pannel
"""
MaskThresh = 0.05;      # we will throw away some pixels that are too dim, this is the fraction of the peak intensity below which to discard.
SqAnalysis = 1;
TriAnalysis = 0;import numpy as np
existROIinfo = 0;       # whether or not to load ROIinfo.m file (1 = yes; 0 = no)
method = 1;
backGND_corr = 1;
Photobleaching_corr=0;
binNum=20   # Histogram bin number

"""
Define path and file name
"""
filePath = 'G:/MATLAB/patch/2014-08-11 HEK nanoparticles/nanorods 630/fov2 - good/122055_take2 100Hz'
#filePath = 'C:/patch/2014-08-11 HEK nanoparticles/nanorods 630/fov2 - good/122126_take3 100Hz preillu'
#filePath = 'G:/patch/2014-08-11 HEK nanoparticles/nanorods 630/fov2 - good/122126_take3 100Hz preillu'
filePath=filePath+'/'
filePath=filePath.replace('\\','/')

ParameterFileName='experimental_parameters.txt'
MovFile='Sq_camera.bin'


"""
Get DAQ aquisition data
"""
Value=[]
for row in open(filePath+ParameterFileName, 'rb'):  # r is read, b is binary
    ParamName, val = row.split('\t')
    Value.append(float(val))
SqRec = Value[0]           # Square wave recorded? (0 = No; 1 = Yes)
TriRec = Value[1]          # Triangle wave recorded? (0 = No; 1 = Yes)
SqDAQrate = Value[3]       # Square wave DAQ sampling rate (Hz)
SqFreqV = Value[4]         # Square wave frequency (Hz)
SqFrameRate = Value[6]     # Square frame rate (Hz)
TriDAQrate = Value[7]      # Triangle wave DAQ sampling rate (Hz)
TriFreqV = Value[8]        # Triangle wave frequency (Hz)
TriFrameRate = Value[10]   # Triangle frame rate (Hz)
ncol = int(Value[14])
nrow = int(Value[15])
DAQ_Sq = 'Sqrwave_DAQ.txt'
DAQ_Tri = 'Triwave_DAQ.txt'
Cam_Sq = 'Sq_camera.bin'
Cam_Tri = 'Tri_camera.bin'
Cam_fname= 'Sq_camera.bin'


"""
Import movie & Downsampled t, V, I
"""
mov, frame = ReadBinMov.ReadBinMov(filePath+MovFile, nrow, ncol)   # mov is [row col frame]
  # the first few frames are often screwed, replace with subsequent frame
  # !these should be excluded from analysis!
mov[:4,:,:]=np.tile(mov[4,:,:], (4,1,1))   # new axis make [row col] to [row col 1]
mov=np.double(mov)
params=[mov, SqFrameRate, MaskThresh, filePath, DAQ_Sq, SqDAQrate, method]

t, V, I = DAQ.DAQread(filePath+DAQ_Sq, SqFrameRate, SqDAQrate, SqFreqV, frame)   # Extract downsampled t, V and I
period=int(SqFrameRate/SqFreqV)


"""
Define background and bgnd correction and photobleaching correction
"""
refimg=np.mean(mov, axis=0)
plt.figure(1)
plt.subplot(211)
plt.imshow(refimg, cmap=cm.Greys_r, vmin=refimg.min(), vmax=refimg.max())
plt.show()
plt.title('Averaged Image')

bg, bg_pts = polygon.mean_polygon(mov, refimg)   # Define background by polygon input and return its mean intensity of background

plt.subplot(212)
line_bgnd, = plt.plot(bg, label="BGND")
plt.legend(bbox_to_anchor=(0.8, 1.2), loc=2, borderaxespad=0.)
plt.show()
plt.title('Averaged Background of ROI')

bg_smooth=smooth.movingaverage(bg, frame/10)

plt.plot(bg_smooth, 'r', label="Smoothen BGND")
plt.xlabel('Frame')
plt.legend(bbox_to_anchor=(0.8, 1.2), loc=2, borderaxespad=0.)

mov_bg_cr=mov-np.tile(bg_smooth[:,newaxis,newaxis], (1,nrow,ncol))  # bgnd corrected 3d movie
mov_bg_cr1 = np.sum(np.sum(mov_bg_cr, axis=1), axis=1)/(ncol*nrow)  # bgnd corrected avg intensity
mov_pb=smooth.movingaverage(mov_bg_cr1, frame/10)   # photobleaching

pb_constant=np.polyfit(t,mov_pb, 1)
pbleach=np.polyval(pb_constant,t)
pbc=pb_constant[1]/pbleach

plt.figure(2)
line_smoothen, = plt.plot(t,mov_pb, label="Smoothen I")
line_pbleaching, = plt.plot(t,pbleach, label="Photobleaching")
line_pb_correct_I, = plt.plot(t,multiply(mov_pb,pbc), label="P.B. corrected I")
plt.title('Photobleaching correction')
plt.legend([line_smoothen, line_pbleaching, line_pb_correct_I], ["Smoothen I", "Photobleaching", "P.B. corrected I"], loc=2)
plt.xlabel('time [s]')
plt.show()


if backGND_corr == 1:
    if Photobleaching_corr == 1:
        mov_f=np.zeros((frame,nrow,ncol))
        for i in range(frame):
            mov_f[i,:,:]=mov_bg_cr[i,:,:]*pbc[i] 
    else: 
        mov_f=mov_bg_cr
else:
    mov_f=mov

            
"""
Define points (QDs) of interest, and their nearest peak position
"""
plt.figure(1)
plt.subplot(211)
fig = gcf()
fig.canvas.manager.window.raise_()

pts = point.pIO(mov)
pts=np.array(pts)
pts_new = point.localmax(refimg, pts)
npoint=size(pts_new[:,0])



"""
Extracting mean intensity of 5 X 5 mask in your points of interest
"""
spot_intensity=np.zeros((frame,npoint))
for n in range(npoint):
    for i in range(frame):
        temp=point.mask(mov_f[i,:,:], pts_new[n,:])
        spot_intensity[i,n]=temp.mean()

plt.figure(3)
plt.subplot(212)
fig = gcf()
fig.canvas.manager.window.raise_()

"""
Apply threshold
"""
threshold=np.zeros((npoint))
thresh_line=np.zeros((frame,npoint))
for i in range(npoint):
    plt.text(0.1, 0.9*spot_intensity[:,i].max(), "Click the treshold point", color='red', fontsize=14, fontweight='bold')
    plot(t, spot_intensity[:,i])
    temp=np.array(ginput(1))
    threshold[i]=temp[0,1]
    cla()

"""
Data Analysis for spot
"""
dFF=np.zeros((npoint))  # Unthresholded, averaged delta (F) / F 
dFF_th=np.zeros((npoint))   # Thresholded, averaged delta (F) / F unit 
avg=np.zeros((npoint,period)) # Average Intensity over multiple cycle, unthresholded
filted_avg=np.zeros((npoint,period)) # Average Intensity over multiple cycle, thresholded
sortedI=np.zeros((frame, npoint))  # Sort spot_intensity such that 1st phase to 0:frame/2 and 2nd phase to frame/2 : frame
diff=np.zeros((frame/period,npoint))  #list of delta F (I1st - I2nd) per cycle, unthresholded
diff_th1=np.zeros((frame/period,npoint)) #list of delta F (I1st - I2nd) per cycle, thresholded
diff_th2=np.zeros((frame/period,npoint)) #list of delta F (I2nd - I1st (next)) per cycle, thresholded
dff_avg=np.zeros(npoint)
for i in range(npoint):
    filted_avg[i,:], dFF_th[i]  =spotAnalysis.threshold_avgperiod(threshold[i], spot_intensity[:,i], period)
    avg[i,:], dFF[i], sortedI[:,i] =spotAnalysis.avg_period(t, spot_intensity[:,i], period, threshold[i])
    diff[:,i]=spotAnalysis.difference(spot_intensity[:,i], period)
    diff_th1[:,i], diff_th2[:,i], dff_avg[i] = spotAnalysis.filted_diff(spot_intensity[:,i], period, threshold[i])
    
#figNum=multiplot.multiplot(filted_avg)

"""
2nd Analysis
"""
oe_ratio=np.zeros((npoint))
n_diff=np.zeros((binNum, npoint))
bins=np.zeros((binNum+1, npoint))
for i in range(npoint):
    oe_ratio[i] = multiplot.spotAnalysis(refimg, pts[i,:], sortedI[:,i], spot_intensity[:,i], threshold[i], diff_th1[:,i], diff_th2[:,i], filted_avg[i,:], t, binNum)

    
#even, odd, nhist, bins = spotAnalysis.evenodd(diff_th1)   #bins are dI/I


"""
Changeable plot of selected points
"""
#fig, ax = plt.subplots()
#plt.title('Intensity trajectory', fontsize=14, fontweight='bold')
#plt.xlabel('time [s]')
#fig.subplots_adjust(left=0.2)   # subplot starts at 50 % from the left
#lines = plt.plot(t, spot_intensity, visible=False)
#lines_thresh = plt.plot(t, thresh_line, visible=False, linestyle='--')

# Build check button axes  
#rax = plt.axes([0.02, 0.4, 0.13, 0.2], aspect='equal')  # [x start, y start, x span, y span]
#labels = range(npoint)
#labels=map(str, labels)
#nofFalse=False, 
#nofFalse=nofFalse*(npoint)
#check = CheckButtons(rax, labels, nofFalse)
#def func(label):    # button's label
#    i = labels.index(label)    # label is either '2 Hz' or '4 Hz' '6 Hz', so i is 0 or 1 or 2
#    lines[i].set_visible(not lines[i].get_visible())    #if lines[i] is not visible, make it visible
#    lines_thresh[i].set_visible(not lines_thresh[i].get_visible())  
#    fig.canvas.draw()
#check.on_clicked(func)   # when the button is clicked, call func with button label        

"""
even odd gaussian decomposition
"""
#y_even=np.zeros((len(nhist[:,0]), npoint))
#y_odd=np.zeros((len(nhist[:,0]), npoint))
#oe_ratio=np.zeros((npoint))
#for i in range(npoint):
#    y_even[:,i], y_odd[:,i], oe_ratio[i]=fit.EOGdecomposition(nhist[:,i], bins[:,i])
#

#multiplot.VonVoffplot(sortedI, threshold)

#
#plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
#plt.subplot(222)
#line_Current, = plt.plot(tm[:SqFreqV*5], Im[:SqFreqV*5]/30, 'c', label="Current / (30 pA)" )
#line_Voltage, = plt.plot(tm[:SqFreqV*5], Vm[:SqFreqV*5], 'm', label=" Voltage (V)")
#plt.legend(bbox_to_anchor=(0.1, 1.2), loc=2, borderaxespad=0.)
#plt.xlabel('time [s]')


#plt.legend([line_Current, line_Voltage], ["Current / (30 pA)", " Voltage (V)"])
#mov_f = multiply(mov_bg_cr, np.tile(bg_smooth[:,newaxis,newaxis], (1,nrow,ncol))
#data_dir = os.path.abspath(filepath1) + '/'

# Check if data_dir exists
#assert os.path.exists(data_dir), "Path '%s' does not exist." % data_dir

# Load all the filenames in a list
#from glob import glob
#file_list = sorted(glob(data_dir + 'Sq_camera.bin'))


 
    