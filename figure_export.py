# -*- coding: utf-8 -*-
"""
Created on Sun Nov 23 15:44:20 2014

@author: KyoungWon
"""
import matplotlib.gridspec as gridspec
from sub import plot
from scipy.optimize import leastsq


def gaussian(A, x, mu, sigma):
    return A*np.exp(-(x-mu)**2/(2*sigma**2))
def oddGaussian(A, x, mu, sigma):
    return A*np.exp(-(x-mu)**2/(2*sigma**2))*(x-mu)
def resEOG(p, y, x):
    A1, A2, m1, m2, sd1, sd2 = p
    m1=0  # 0 is fixed for zero mean gaussian
    m2=0
    y_fit = gaussian(A1, x, m1, sd1) + oddGaussian(A2, x, m2, sd2)
    err = y - y_fit
    return err
    
    
G = gridspec.GridSpec(3, 2)
plt.figure()
axes_1 = plt.subplot(G[0, 0:2])

plt.rc("font", size=16)
plt.plot(t, spot_intensity[:,2], color='b')
#thresh_line=np.ones(nframes)*threshold
#ax.plot(t, thresh_line, color='r')
plt.axhline(threshold[2], color='r')
#plt.xlabel("Time [s]", fontsize=16)
plt.ylabel("Intensity [a.u]", fontsize=16)

axes_4 = plt.subplot(G[1, 0])
abs_diff1=np.abs(diff_th1[:,2])
abs_diff2=np.abs(diff_th2[:,2])
max_diff=max(np.nanmax(abs_diff1), np.nanmax(abs_diff2))
hist_range=(-1*max_diff, max_diff)
n, bins, patches = axes_4.hist(diff_th1[:,2], binNum, range=hist_range, color='b', normed=False, alpha=0.4, label='1st-2nd')
n_shift, bins, patches = axes_4.hist(diff_th2[:,2], binNum, range=hist_range, color='r', normed=False, alpha=0.4, label='2nd-1st')
#axes_4.set_title('Histogram of $\Delta I/I_{Max}$ , thresholded')
axes_4.axis([max_diff*-1.2, max_diff*1.2, 0, n.max()*1.4])
plt.xticks([-0.6, -0.3, 0, 0.3, 0.6])

axes_6 = plt.subplot(G[1, 1])
axes_6.plot(avg_th[2,:], '-o')
plt.yticks([10.4, 10.6, 10.8])
plt.xticks([0, 1, 2, 3])


axes_5 = plt.subplot(G[2, 0])

onoff=Ndiff1hist[2][7]
offon=Ndiff2hist[2][7]
hmax=max(onoff.max(), offon.max())
hist_range=(-hmax, hmax)
nr, binss, patches0 = axes_5.hist(onoff, binNum, range=hist_range, color='b',  alpha=0.4)
nr, binss, patches01 = axes_5.hist(offon, binNum, range=hist_range, color='r',  alpha=0.4)
dF=(np.mean(onoff)-np.mean(offon))/2
#ax0.set_title("%d cycle, $\Delta I/I_{Max}$ = %.2f [%%]" % (number+2, dF*100) )
axes_5.set_xlim(-0.6, 0.6)
plt.xticks([-0.6, -0.3, 0, 0.3, 0.6])

#ax0.axvline(x=0, linewidth=2, color='k')

axes_7 = plt.subplot(G[2, 1])
xnew=np.zeros((len(bins)-1))
for i in range(len(bins)-1):
    xnew[i]= (bins[i]+bins[i+1])/2
mu = sum(xnew*n)/sum(n)
mu=0
sigma = np.sqrt(sum(n*(xnew-mu)**2)/sum(n))
A1 = n.max()
A2 = 0
p = [A1, A2, mu, mu, sigma, sigma]   # initial guess parameters
y_init=gaussian(A1, xnew, mu, sigma)  + oddGaussian(A2, xnew, mu, sigma)

plsq=leastsq(resEOG, p, args = (n, xnew))

y_est = gaussian(plsq[0][0], xnew, plsq[0][2], plsq[0][4]) + oddGaussian(plsq[0][1], xnew, plsq[0][2] + plsq[0][3], plsq[0][5])
y_even = gaussian(plsq[0][0], xnew, plsq[0][2], plsq[0][4])
y_odd = oddGaussian(plsq[0][1], xnew, plsq[0][2] + plsq[0][3], plsq[0][5])
if plsq[0][1] > 0:
    oe_ratio=y_odd.max()/y_even.max()
else:
    oe_ratio=-1*y_odd.max()/y_even.max()

axes_7.plot(xnew, n, label='Real Data')
#axes_7.plot(xnew, y_init, 'r.', label='Starting Guess')
axes_7.plot(xnew, y_est, 'g.', label='Fitted')
axes_7.plot(xnew, y_even, 'c', label='even Gaussian')
axes_7.plot(xnew, y_odd, 'm', label='odd Gaussian')
#axes_7.set_title('Fitting with Gaussians')
#axes_7.text(-0.1, -0.05, 'Odd /Even Gaussian ratio is %f ' % oe_ratio)


