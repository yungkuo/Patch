"""
Patch-clamp voltage sensing analysis script.

Author: KyoungWon
"""

from __future__ import division

import numpy as np
from scipy.ndimage.filters import gaussian_filter1d
import matplotlib as mpl
import matplotlib.pyplot as plt

#from IPython import get_ipython
#ip = get_ipython()
#ip.magic("reset -f")
#plt.close("all")

from sub import (ReadBinMov, polygon, smooth, point, DAQ, spotAnalysis,
                 multiplot)
from search import bsearch_py
"""
Control Pannel
"""
MaskThresh = 0.05   # we will throw away some pixels that are too dim, this is
                    # the fraction of the peak intensity below which to discard
SqAnalysis = True
TriAnalysis = 0
existROIinfo = 0    # whether or not to load ROIinfo.m file (1 = yes; 0 = no)
method = True
backGND_corr = True
Photobleaching_corr = False  # for background
spot_pbc = False        # for QDs
spot_LPF = True       # LPF for spot
binNum = 20         # Histogram bin number

"""
Define path and file name
"""
filePath = 'G:/MATLAB/patch/2014-08-11 HEK nanoparticles/nanorods 630/fov2 - good/122055_take2 100Hz'
#filePath = 'C:/patch/2014-08-11 HEK nanoparticles/nanorods 630/fov2 - good/122126_take3 100Hz preillu'
#filePath = 'G:/patch/2014-08-11 HEK nanoparticles/nanorods 630/fov2 - good/122055_take2 100Hz'
filePath=filePath+'/'
filePath=filePath.replace('\\','/')

ParameterFileName = 'experimental_parameters.txt'
MovFile = 'Sq_camera.bin'


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
 # mov is [row col frame]
mov, frame = ReadBinMov.ReadBinMov(filePath+MovFile, nrow, ncol)
  # the first few frames are often screwed, replace with subsequent frame
  # !these should be excluded from analysis!

mov[:4, :, :] = np.tile(mov[4, :, :], (4, 1, 1))
mov = mov.astype('float64')
params = [mov, SqFrameRate, MaskThresh, filePath, DAQ_Sq, SqDAQrate, method]

# Extract downsampled t, V and I
t, V, I = DAQ.DAQread(filePath+DAQ_Sq, SqFrameRate, SqDAQrate, SqFreqV, frame)
period = int(np.round(SqFrameRate / SqFreqV))
assert period == SqFrameRate / SqFreqV, "SqFrameRate / SqFreqV should be int"

"""
Define next V is whether rising (1) or falling (-1) or staying phase (0)
"""
dV = np.diff(V)
rp = [] # rising phase
sp = [] # staying phase
fp = [] # falling phase
for i in range(len(dV)):
    if dV[i] >= 0.05:
        rp.append(i)
    elif dV[i] <= -0.05:
        fp.append(i)
    else:
        sp.append(i)

"""
Define background and bgnd correction and photobleaching correction
"""
refimg = np.mean(mov, axis=0)
plt.figure(1)
plt.subplot(211)
plt.imshow(refimg, cmap=mpl.cm.Greys_r, vmin=refimg.min(), vmax=refimg.max())
plt.title('Averaged Image')
plt.show()


# Define background region by polygon input and return its mean intensity
bg, bg_pts = polygon.mean_polygon(mov, refimg)

plt.subplot(212)
line_bgnd, = plt.plot(bg, label="BGND")
plt.legend(bbox_to_anchor=(0.8, 1.2), loc=2, borderaxespad=0.)
plt.title('Averaged Background of ROI')
plt.show()


#bg_smooth=smooth.moving_average(bg, frame/10)
bg_smooth = gaussian_filter1d(bg, sigma=100)
plt.plot(bg_smooth, 'r', label="Smoothen BGND")
plt.xlabel('Frame')
plt.legend(bbox_to_anchor=(0.8, 1.2), loc=2, borderaxespad=0.)

# bgnd corrected 3d movie
mov_bg_cr = mov-np.tile(bg_smooth[:, np.newaxis, np.newaxis], (1,nrow,ncol))
# bgnd corrected avg intensity
mov_bg_cr1 = np.sum(np.sum(mov_bg_cr, axis=1), axis=1)/(ncol*nrow)


if backGND_corr == True:
    if Photobleaching_corr == True:
        mov_f=np.zeros((frame,nrow,ncol))
        mov_pb=smooth.moving_average(mov_bg_cr1, frame/10)   # photobleaching

        pb_constant=np.polyfit(t,mov_pb, 1)
        pbleach=np.polyval(pb_constant,t)
        pbc=pb_constant[1]/pbleach

        plt.figure(2)
        line_smoothen, = plt.plot(t,mov_pb, label="Smoothen I")
        line_pbleaching, = plt.plot(t,pbleach, label="Photobleaching")
        line_pb_correct_I, = plt.plot(t, np.multiply(mov_pb,pbc),
                                      label="P.B. corrected I")
        plt.title('Photobleaching correction')
        plt.legend([line_smoothen, line_pbleaching, line_pb_correct_I],
                   ["Smoothen I", "Photobleaching", "P.B. corrected I"],
                   loc=2)
        plt.xlabel('time [s]')
        plt.show()

        for i in range(frame):
            mov_f[i,:,:]=mov_bg_cr[i,:,:]*pbc[i]
    else:
        mov_f=mov_bg_cr
else:
    mov_f=mov

del mov, mov_bg_cr, mov_bg_cr1
"""
Define points (QDs) of interest, and their nearest peak position
"""
plt.figure(1)
plt.subplot(211)
fig = plt.gcf()
fig.canvas.manager.window.raise_()

pts = point.pIO(mov_f)
pts=np.array(pts)
pts_new = point.localmax(refimg, pts)
npoint = np.size(pts_new[:, 0])


"""
Extracting mean intensity of 5 X 5 mask in your points of interest
"""
spot_intensity=np.zeros((frame,npoint))
for n in range(npoint):
    for i in range(frame):
        temp = point.mask(mov_f[i,:,:], pts_new[n,:])
        spot_intensity[i,n] = temp.mean()

if spot_pbc == True:
    for n in range(npoint):
        for i in range(frame):
            temp = smooth.moving_average(spot_intensity[:,n], frame/10)
            pb_constant = np.polyfit(t, temp, 1)
            pbleach = np.polyval(pb_constant, t)
            pbc = pb_constant[1]/pbleach
            spot_intensity[:,n] = np.multiply(spot_intensity[:, n], pbc)
elif spot_LPF == True:
    for i in range(npoint):
        temp = gaussian_filter1d(spot_intensity[:,i], sigma=100)
        spot_intensity[:,i] = spot_intensity[:,i]-temp
        spot_intensity[:,i]=spot_intensity[:,i]+abs(spot_intensity[:,i].min())

plt.figure(3)
plt.subplot(212)
fig = plt.gcf()
fig.canvas.manager.window.raise_()

"""
Apply threshold
"""
threshold=np.zeros((npoint))
thresh_line=np.zeros((frame,npoint))
for i in range(npoint):
    plt.text(0.1, 0.9*spot_intensity[:,i].max(), "Click the treshold point",
             color='red', fontsize=14, fontweight='bold')
    plt.plot(t, spot_intensity[:,i])
    temp=np.array(plt.ginput(1))
    threshold[i]=temp[0,1]
    plt.cla()

"""
Data Analysis for spot
"""
meanI, mean_thI = spotAnalysis.meanspotI(spot_intensity, period, threshold)
mean_thI = [np.array(m) for m in mean_thI]

# Unthresholded, averaged delta (F) / F
dFF = np.zeros((npoint))
# Thresholded, averaged delta (F) / F unit
dFF_th = np.zeros((npoint))
# Average Intensity over multiple cycle, unthresholded
avg = np.zeros((npoint,period))
# Average Intensity over multiple cycle, thresholded
filted_avg = np.zeros((npoint,period))
# Sort spot_intensity such that 1st phase to 0:frame/2 and 2nd phase to
# frame/2 : frame
#sortedI = np.zeros((frame, npoint))
##list of delta F (I1st - I2nd) per cycle, unthresholded
#diff=np.zeros((frame/period,npoint))
#list of delta F (I1st - I2nd) per cycle, thresholded
diff_th1 = np.zeros((frame/period,npoint))
#list of delta F (I2nd - I1st (next)) per cycle, thresholded
diff_th2 = np.zeros((frame/period,npoint))
dff_avg = np.zeros(npoint)
rising = np.zeros((len(rp), npoint))
staying = np.zeros((len(sp), npoint))
falling = np.zeros((len(fp), npoint))
Rdiff = []
Sdiff = []
Fdiff = []
selected_rsf = []
Maxcorr1 = []
Maxcorr2 = []
for i in range(npoint):
    filted_avg[i,:], dFF_th[i], rtemp, stemp, ftemp  = \
            spotAnalysis.threshold_avgperiod(threshold[i], spot_intensity[:,i],
                                             period, rp, sp, fp)
    Rdiff.append(rtemp)
    Sdiff.append(stemp)
    Fdiff.append(ftemp)

    avg[i,:], dFF[i] = spotAnalysis.avg_period(t, spot_intensity[:,i],
                        period, threshold[i])
    #diff[:,i]=spotAnalysis.difference(spot_intensity[:,i], period)

    diff_th1[:,i], diff_th2[:,i], dff_avg[i] = \
        spotAnalysis.filted_diff(spot_intensity[:,i], period, threshold[i])

    rising[:,i], staying[:,i], falling[:,i] = \
        spotAnalysis.rf_selective_diff(spot_intensity[:,i], rp, sp, fp)

    #temp1, temp2 = spotAnalysis.correlating_cycle(spot_intensity[:,i], frame, period)
    temp1, temp2 = spotAnalysis.norm_corr(spot_intensity[:,i], threshold[i], period)
    Maxcorr1.append(temp1)
    print("%d QD 1st > 2nd has %d consequetive cycle" % (i, len(temp1)/period))
    Maxcorr2.append(temp2)
    print("%d QD 1st < 2nd has %d consequetive cycle" % (i, len(temp2)/period))
    #print("%d'th QD" % len(temp1)/period))
    #print(len(temp2)/period)
#figNum=multiplot.multiplot(filted_avg)



"""
2nd Analysis
"""

oe_ratio = np.zeros((npoint))
n_diff = np.zeros((binNum, npoint))
bins = np.zeros((binNum+1, npoint))
Ndiff1hist = []
Ndiff2hist = []
avg1 = []
avg2 = []
barray=[]
for i in range(npoint):
    btemp=spotAnalysis.burstarray(mean_thI[i])  # barray is (F[n+2]-2F[n+1]+F[n])/Favg for burst analysis
    barray.append(btemp)
    temp1, tavg1 = spotAnalysis.multiNdF(diff_th1[:,i])
    temp2, tavg2 = spotAnalysis.multiNdF(diff_th2[:,i])
    Ndiff1hist.append(temp1)
    Ndiff2hist.append(temp2)
    avg1.append(tavg1)
    avg2.append(tavg2)
    #Ndiff2.append(spotAnalysis.multiNdF(diff_th2[:,i]))
    oe_ratio[i] = multiplot.spotAnalysis(refimg, pts[i,:], meanI[:,i],
                                spot_intensity[:,i], threshold[i],
                                diff_th1[:,i], diff_th2[:,i], filted_avg[i,:],
                                t, binNum, rising[:,i], staying[:,i],
                                falling[:,i], Rdiff[i], Sdiff[i], Fdiff[i])



NdF = []
for i in range(npoint):
    temp_dF = multiplot.Ndiffhisto(Ndiff1hist[i], Ndiff2hist[i], binNum)
    NdF.append(temp_dF)
#even, odd, nhist, bins = spotAnalysis.evenodd(diff_th1)   #bins are dI/I

burstdata=[]
for i in range(npoint):
    a=bsearch_py.burstsearch_py(np.array(barray[i]), 10, 0.3)
    burstdata.append(a)
    #for j in a


result=np.zeros((npoint, 12))
for i in range(npoint):
    result[i, 0]=len(Maxcorr1[i])/period  # Maximum consequetive on > off cycle
    result[i, 1]=len(Maxcorr2[i])/period  # Maximum consequetive on < off cycle
    result[i, 3:11]=NdF[i]   # dF/F,  N cycle averaged
result[:,2]=dff_avg  # dF/F
result[:,11]=oe_ratio  # odd/even gaussian peak

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



