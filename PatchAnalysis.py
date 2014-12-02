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
backGND_corr = False
Photobleaching_corr = False  # for background
spot_pbc = False        # for QDs
spot_LPF = True       # LPF for spot
binNum = 20         # Histogram bin number
shift_frame = False  # control exp for
"""
Define path and file name
"""
filePath = 'G:/MATLAB/patch/2014-08-11 HEK nanoparticles/nanorods 630/fov1/114640_take1 100Hz'
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
SqFrameRate = Value[6]     # Square nframes rate (Hz)
TriDAQrate = Value[7]      # Triangle wave DAQ sampling rate (Hz)
TriFreqV = Value[8]        # Triangle wave frequency (Hz)
TriFrameRate = Value[10]   # Triangle nframes rate (Hz)
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
mov, nframes = ReadBinMov.ReadBinMov(filePath+MovFile, nrow, ncol)
  # the first few frames are often screwed, replace with subsequent nframes
  # !these should be excluded from analysis!

# Extract downsampled t, V and I
t, V, I = DAQ.DAQread(filePath+DAQ_Sq, SqFrameRate, SqDAQrate, SqFreqV, nframes)
period = int(np.round(SqFrameRate / SqFreqV))
assert period == SqFrameRate / SqFreqV, "SqFrameRate / SqFreqV should be int"

if shift_frame:
    mov = mov[1:-3, :, :]
    nframes = nframes - period
    t = t[1:-3]
    V = V[1:-3]
    I = I[1:-3]
mov[:4, :, :] = np.tile(mov[4, :, :], (4, 1, 1))
mov = mov.astype('float64')

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


#bg_smooth=smooth.moving_average(bg, nframes/10)
bg_smooth = gaussian_filter1d(bg, sigma=100)
plt.plot(bg_smooth, 'r', label="Smoothen BGND")
plt.xlabel('Frame')
plt.legend(bbox_to_anchor=(0.8, 1.2), loc=2, borderaxespad=0.)

# bgnd corrected 3d movie
mov_bg_cr = mov-np.tile(bg_smooth[:, np.newaxis, np.newaxis], (1,nrow,ncol))
# bgnd corrected avg intensity
mov_bg_cr1 = np.sum(np.sum(mov_bg_cr, axis=1), axis=1)/(ncol*nrow)   #1d intensity profile


if backGND_corr:
    if Photobleaching_corr == True:
        mov_f=np.zeros((nframes,nrow,ncol))
        mov_pb=smooth.moving_average(mov_bg_cr1, nframes/10)   # photobleaching

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

        for i in range(nframes):
            mov_f[i,:,:]=mov_bg_cr[i,:,:]*pbc[i]
    else:
        mov_f=mov_bg_cr
else:
    mov_f=mov

#del mov, mov_bg_cr, mov_bg_cr1
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
spot_intensity=np.zeros((nframes,npoint))
for n in range(npoint):
    for i in range(nframes):
        temp = point.mask(mov_f[i,:,:], pts_new[n,:])
        spot_intensity[i,n] = temp.mean()
del mov_f

if spot_pbc:
    for n in range(npoint):
        for i in range(nframes):
            temp = smooth.moving_average(spot_intensity[:,n], nframes/10)
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
thresh_line=np.zeros((nframes,npoint))
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

dFF = np.zeros((npoint)) # Unthresholded, averaged delta (F) / F
dFF_th = np.zeros((npoint)) # Thresholded, averaged delta (F) / F unit
avg = np.zeros((npoint,period)) # Avg Intensity over multiple cycle, unthresholded
avg_th = np.zeros((npoint,period)) # Average Intensity over multiple cycle, thresholded

#diff=np.zeros((nframes/period,npoint)) #list dF (over dV) for each cycle, unthresholded
diff_th1 = np.zeros((nframes/period,npoint)) #list of dF for each cycle, thresholded
#list of delta F (I2nd - I1st (next)) per cycle, thresholded
diff_th2 = np.zeros((nframes/period,npoint))

dff_avg = np.zeros(npoint)
Fstd=np.zeros(npoint)
rising = np.zeros((len(rp), npoint))
staying = np.zeros((len(sp), npoint))
falling = np.zeros((len(fp), npoint))
Rdiff = []
Sdiff = []
Fdiff = []
selected_rsf = []
Maxcorr1 = []
Maxcorr2 = []
barray=[] # input for burst search
burstdata=[] # output of burst search

for i in range(npoint):
    avg_th[i,:], dFF_th[i], rtemp, stemp, ftemp  = \
            spotAnalysis.threshold_avgperiod(threshold[i], spot_intensity[:,i],
                                             period, rp, sp, fp)
    Rdiff.append(rtemp)
    Sdiff.append(stemp)
    Fdiff.append(ftemp)

    avg[i,:], dFF[i] = spotAnalysis.avg_period(t, spot_intensity[:,i],
                        period, threshold[i])
    #diff[:,i]=spotAnalysis.difference(spot_intensity[:,i], period)

    diff_th1[:,i], diff_th2[:,i], dff_avg[i], Fstd[i] = \
        spotAnalysis.filted_diff(spot_intensity[:,i], period, threshold[i])

    rising[:,i], staying[:,i], falling[:,i] = \
        spotAnalysis.rf_selective_diff(spot_intensity[:,i], rp, sp, fp)

    temp1, temp2 = spotAnalysis.correlating_cycle(spot_intensity[:,i], nframes, period)
    #temp1, temp2 = spotAnalysis.norm_corr(spot_intensity[:,i], threshold[i], period)
    Maxcorr1.append(temp1)
    print("%d QD 1st > 2nd has %d consequetive cycle" % (i, len(temp1)/period))
    Maxcorr2.append(temp2)
    print("%d QD 1st < 2nd has %d consequetive cycle" % (i, len(temp2)/period))

    btemp=bsearch_py.burstarray(meanI[:,i])  # barray is (F[n+2]-2F[n+1]+F[n])/Favg for burst analysis
    barray.append(btemp)
    a=bsearch_py.burstsearch_py(np.array(btemp), 10,3)
    burstdata.append(a)
    #print("%d'th QD" % len(temp1)/period))
    #print(len(temp2)/period)
#figNum=multiplot.multiplot(avg_th)
barray = [np.array(m) for m in barray]
burstdata = [np.array(m) for m in burstdata]


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
burstmask=np.zeros((nframes, npoint), dtype=bool)
for i in range(npoint):
    burstmask[:,i] = bsearch_py.burstmask(nframes, period, burstdata[i])
    temp1, tavg1 = spotAnalysis.multiNdF(diff_th1[:,i])
    temp2, tavg2 = spotAnalysis.multiNdF(diff_th2[:,i])
    Ndiff1hist.append(temp1)
    Ndiff2hist.append(temp2)
    avg1.append(tavg1)
    avg2.append(tavg2)
    #Ndiff2.append(spotAnalysis.multiNdF(diff_th2[:,i]))
    oe_ratio[i] = multiplot.spotAnalysis(refimg, pts[i,:], meanI[:,i],
                                spot_intensity[:,i], threshold[i],
                                diff_th1[:,i], diff_th2[:,i], avg_th[i,:],
                                t, binNum, rising[:,i], staying[:,i],
                                falling[:,i], Rdiff[i], Sdiff[i], Fdiff[i])



NdF = []
Nstd = []
for i in range(npoint):
    temp_dF, temp_std = multiplot.Ndiffhisto(Ndiff1hist[i], Ndiff2hist[i], binNum)
    NdF.append(temp_dF)
    Nstd.append(temp_std)
even, odd, nhist, bins = spotAnalysis.evenodd(diff_th1)   #bins are dI/I

result=np.zeros((npoint, 21))  # summarize the data for excel input
for i in range(npoint):
    result[i, 0]=len(Maxcorr1[i])/period  # Maximum consequetive on > off cycle
    result[i, 1]=len(Maxcorr2[i])/period  # Maximum consequetive on < off cycle
    result[i, 4:12]=NdF[i]   # dF/F,  N cycle averaged
    result[i, 13:21]=Nstd[i]
result[:,3]=dff_avg  # dF/F
result[:,12]=Fstd
result[:,2]=oe_ratio  # odd/even gaussian peak

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



