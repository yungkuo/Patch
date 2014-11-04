# -*- coding: utf-8 -*-
"""
Created on Fri Oct 03 17:13:56 2014

@author: Philip
"""
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
import matplotlib.gridspec as gridspec



def spotAnalysis(refimg, pts, sortedI, intensity, threshold, diff1, diff2, filted_avg, t, binNum):
    frame=len(intensity)
    G = gridspec.GridSpec(3, 3)
    fig =figure(figsize=(13, 9))
    axes_0 = subplot(G[0, 0])    
    axes_0.imshow(refimg, cmap=cm.Greys_r, vmin=refimg.min(), vmax=refimg.max())
    axes_0.plot(pts[1], pts[0], 'ro')
    ncol=len(refimg[0,:])
    nrow=len(refimg[:,0])
    axes_0.set_xlim([0,ncol])
    axes_0.set_ylim([nrow,0])
    
    
    axes_1 = subplot(G[0, 1:])
    axes_1.plot(t, intensity, color='b')
    thresh_line=np.ones(frame)*threshold
    axes_1.plot(t, thresh_line, color='r')
    axes_1.set_xlabel("Time [s]")
    axes_1.xaxis.set_label_coords(1.05, -0.025)

    axes_1.set_title("Intensity trajectory")
    
    axes_2 = subplot(G[1, 0])
    Von=sortedI[:frame/2]
    Voff=sortedI[frame/2:]
    #f, ax=plt.subplots(2, sharex=True)
    maxbin=max(max(Von), max(Voff))
    minbin=min(min(Von), min(Voff))
    hist_range=(minbin, maxbin)
    n1, bins1, patches = axes_2.hist(Von, binNum, range=hist_range, color='b', normed=True, alpha=0.5, label='1st phase')  
    n2, bins2, patches = axes_2.hist(Voff, binNum, range=hist_range, color='r', normed=True, alpha=0.5, label='2nd phase')  
    axes_2.axis([minbin*0.5, maxbin*1.2, 0, n1.max()*1.2])
    axes_2.set_title("Histogram of $I_{Von}$ and $I_{Voff}$")            
        
    axes_3 = subplot(G[2, 0])
    axes_3.bar(bins1[:-1], n1-n2, width=bins1[1]-bins1[0])  
    axes_3.set_xlim(minbin*0.5, maxbin*1.2)
    axes_3.axvline(x=threshold, linewidth=4, color='r')
    axes_3.set_title("$I_{Von}$ - $I_{Voff}$")        
    #plt.xlabel("Intensity, [a.u]")
    
    axes_4 = subplot(G[1, 1])
    abs_diff1=abs(diff1)
    abs_diff2=abs(diff2)
    max_diff=max(np.nanmax(abs_diff1), np.nanmax(abs_diff2))
    hist_range=(-1*max_diff, max_diff)
    n, bins, patches = axes_4.hist(diff1, binNum, range=hist_range, color='b', normed=True, alpha=0.5, label='1st-2nd')
    n_shift, bins, patches = axes_4.hist(diff2, binNum, range=hist_range, color='r', normed=True, alpha=0.5, label='2nd-1st')
    axes_4.set_title('Histogram of $\Delta I/I_{Max}$ , thresholded')
    axes_4.axis([max_diff*-1.2, max_diff*1.2, 0, n.max()*1.2])

    axes_5 = subplot(G[2, 1])
    axes_5.bar(bins[:-1], n-n_shift, width=bins[1]-bins[0])  
    
    axes_5.set_xlim(-1.2*max_diff, max_diff*1.2)
    axes_5.set_title("Cycle - shifted cycle")        
    #plt.xlabel("Intensity, [a.u]")
    
    
    axes_6 = subplot(G[1, 2])
    axes_6.plot(filted_avg)
    axes_6.set_title('Averaged Sigal over multiple cycle')

    axes_7 = subplot(G[2, 2])
    xnew=np.zeros((len(bins)-1))
    for i in range(len(bins)-1):
        xnew[i]= (bins[i]+bins[i+1])/2
    mu = sum(xnew*n)/sum(n)           
    mu=0  
    sigma = sqrt(sum(n*(xnew-mu)**2)/sum(n))
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
    axes_7.plot(xnew, y_init, 'r.', label='Starting Guess')
    axes_7.plot(xnew, y_est, 'g.', label='Fitted')
    axes_7.plot(xnew, y_even, 'c', label='even Gaussian')
    axes_7.plot(xnew, y_odd, 'm', label='odd Gaussian')
    axes_7.set_title('Fitting with Gaussians')
    axes_7.text(-0.1, -0.05, 'Odd /Even Gaussian ratio is %f ' % oe_ratio)
    #axes_7.legend()
    #axes_7.show()        
    
    return oe_ratio 

def gaussian(A, x, mu, sigma):
    return A*exp(-(x-mu)**2/(2*sigma**2))
def oddGaussian(A, x, mu, sigma):
    return A*exp(-(x-mu)**2/(2*sigma**2))*(x-mu)
def resEOG(p, y, x):
    A1, A2, m1, m2, sd1, sd2 = p
    m1=0  # 0 is fixed for zero mean gaussian
    m2=0
    y_fit = gaussian(A1, x, m1, sd1) + oddGaussian(A2, x, m2, sd2)
    err = y - y_fit
    return err
    
    
def multiplot(result):
    npoint=len(result[:,0])
    spot_quotient=npoint/9   
    spot_reminder=mod(npoint,9)
      
    for i in range(spot_quotient+1):
        fig = plt.figure()
        fig.suptitle('Averaged Sigal over multiple cycle', fontsize=14, fontweight='bold')

        if i  == spot_quotient:
            if spot_reminder == 1:
                plt.plot(result[0+i*9,:])
            elif spot_reminder ==2:
                plt.subplot(211)
                plt.plot(result[0+i*9,:])
                plt.subplot(212)
                plt.plot(result[1+i*9,:])
            elif spot_reminder == 3:
                plt.subplot(221)
                plt.plot(result[0+i*9,:])
                plt.subplot(222)
                plt.plot(result[1+i*9,:])
                plt.subplot(223)
                plt.plot(result[2+i*9,:])        
            elif spot_reminder == 4:
                plt.subplot(221)
                plt.plot(result[0+i*9,:])
                plt.subplot(222)
                plt.plot(result[1+i*9,:])
                plt.subplot(223)
                plt.plot(result[2+i*9,:])
                plt.subplot(224)
                plt.plot(result[3+i*9,:])
            elif spot_reminder == 5:
                plt.subplot(321)
                plt.plot(result[0+i*9,:])
                plt.subplot(322)
                plt.plot(result[1+i*9,:])
                plt.subplot(323)
                plt.plot(result[2+i*9,:])
                plt.subplot(324)
                plt.plot(result[3+i*9,:])
                plt.subplot(325)
                plt.plot(result[4+i*9,:])
            elif spot_reminder == 6:
                plt.subplot(321)
                plt.plot(result[0+i*9,:])
                plt.subplot(322)
                plt.plot(result[1+i*9,:])
                plt.subplot(323)
                plt.plot(result[2+i*9,:])
                plt.subplot(324)
                plt.plot(result[3+i*9,:])
                plt.subplot(325)
                plt.plot(result[4+i*9,:])
                plt.subplot(326)
                plt.plot(result[5+i*9,:])
            elif spot_reminder == 7:
                plt.subplot(331)
                plt.plot(result[0+i*9,:])
                plt.subplot(332)
                plt.plot(result[1+i*9,:])
                plt.subplot(333)
                plt.plot(result[2+i*9,:])
                plt.subplot(334)
                plt.plot(result[3+i*9,:])
                plt.subplot(335)
                plt.plot(result[4+i*9,:])
                plt.subplot(336)
                plt.plot(result[5+i*9,:])
                plt.subplot(337)
                plt.plot(result[6+i*9,:])
            elif spot_reminder == 8:
                plt.subplot(331)
                plt.plot(result[0+i*9,:])
                plt.subplot(332)
                plt.plot(result[1+i*9,:])
                plt.subplot(333)
                plt.plot(result[2+i*9,:])
                plt.subplot(334)
                plt.plot(result[3+i*9,:])
                plt.subplot(335)
                plt.plot(result[4+i*9,:])
                plt.subplot(336)
                plt.plot(result[5+i*9,:])
                plt.subplot(337)
                plt.plot(result[6+i*9,:])
                plt.subplot(338)
                plt.plot(result[7+i*9,:])
            else:
                plt.close()
                i=i-1
        else:
            plt.subplot(331)
            plt.plot(result[0+i*9,:])
            plt.subplot(332)
            plt.plot(result[1+i*9,:])
            plt.subplot(333)
            plt.plot(result[2+i*9,:])
            plt.subplot(334)
            plt.plot(result[3+i*9,:])
            plt.subplot(335)
            plt.plot(result[4+i*9,:])
            plt.subplot(336)
            plt.plot(result[5+i*9,:])
            plt.subplot(337)
            plt.plot(result[6+i*9,:])
            plt.subplot(338)
            plt.plot(result[7+i*9,:])
            plt.subplot(339)
            plt.plot(result[8+i*9,:])
    figNum=4+i+1
    return figNum
        
def multihist(diff, threshold):       
    npoint=len(diff[0,:])
    spot_quotient=npoint/9   
    spot_reminder=mod(npoint,9)
    hist_range=[]  # list
    for i in range(npoint):
        abs_diff=abs(diff[:,i])
        hist_range.append((-1*np.nanmax(abs_diff), np.nanmax(abs_diff)))
        
      
    for i in range(spot_quotient+1):
        fig = plt.figure()
        if threshold==True:
            fig.suptitle('Histogram of $\Delta I/I_{Max}$, thresholded & Averaged over period', fontsize=14, fontweight='bold')
        else :
            fig.suptitle('Histogram of $\Delta I/I_{Max}$, unthresholded & Averaged over period', fontsize=14, fontweight='bold')
        if i  == spot_quotient:
            if spot_reminder == 1:
                n, bins, patches = plt.hist(diff[:,0+i*9], 30, range=hist_range[0+i*9])
            elif spot_reminder ==2:
                plt.subplot(211)
                n, bins, patches = plt.hist(diff[:,0+i*9], 30, range=hist_range[0+i*9])
                plt.subplot(212)
                n, bins, patches = plt.hist(diff[:,1+i*9], 30, range=hist_range[1+i*9])
            elif spot_reminder == 3:
                plt.subplot(221)
                n, bins, patches = plt.hist(diff[:,0+i*9], 30, range=hist_range[0+i*9])
                plt.subplot(222)
                n, bins, patches = plt.hist(diff[:,1+i*9], 30, range=hist_range[1+i*9])
                plt.subplot(223)
                n, bins, patches = plt.hist(diff[:,2+i*9], 30, range=hist_range[2+i*9])   
            elif spot_reminder == 4:
                plt.subplot(221)
                n, bins, patches = plt.hist(diff[:,0+i*9], 30, range=hist_range[0+i*9])
                plt.subplot(222)
                n, bins, patches = plt.hist(diff[:,1+i*9], 30, range=hist_range[1+i*9])
                plt.subplot(223)
                n, bins, patches = plt.hist(diff[:,2+i*9], 30, range=hist_range[2+i*9])   
                plt.subplot(224)
                n, bins, patches = plt.hist(diff[:,3+i*9], 30, range=hist_range[3+i*9])   
            elif spot_reminder == 5:
                plt.subplot(321)
                n, bins, patches = plt.hist(diff[:,0+i*9], 30, range=hist_range[0+i*9])
                plt.subplot(322)
                n, bins, patches = plt.hist(diff[:,1+i*9], 30, range=hist_range[1+i*9])
                plt.subplot(323)
                n, bins, patches = plt.hist(diff[:,2+i*9], 30, range=hist_range[2+i*9])   
                plt.subplot(324)
                n, bins, patches = plt.hist(diff[:,3+i*9], 30, range=hist_range[3+i*9])  
                plt.subplot(325)
                n, bins, patches = plt.hist(diff[:,4+i*9], 30, range=hist_range[4+i*9])  
            elif spot_reminder == 6:
                plt.subplot(321)
                n, bins, patches = plt.hist(diff[:,0+i*9], 30, range=hist_range[0+i*9])
                plt.subplot(322)
                n, bins, patches = plt.hist(diff[:,1+i*9], 30, range=hist_range[1+i*9])
                plt.subplot(323)
                n, bins, patches = plt.hist(diff[:,2+i*9], 30, range=hist_range[2+i*9])   
                plt.subplot(324)
                n, bins, patches = plt.hist(diff[:,3+i*9], 30, range=hist_range[3+i*9])  
                plt.subplot(325)
                n, bins, patches = plt.hist(diff[:,4+i*9], 30, range=hist_range[4+i*9])  
                plt.subplot(326)
                n, bins, patches = plt.hist(diff[:,5+i*9], 30, range=hist_range[5+i*9])  
            elif spot_reminder == 7:
                plt.subplot(331)
                n, bins, patches = plt.hist(diff[:,0+i*9], 30, range=hist_range[0+i*9])
                plt.subplot(332)
                n, bins, patches = plt.hist(diff[:,1+i*9], 30, range=hist_range[1+i*9])
                plt.subplot(333)
                n, bins, patches = plt.hist(diff[:,2+i*9], 30, range=hist_range[2+i*9])   
                plt.subplot(334)
                n, bins, patches = plt.hist(diff[:,3+i*9], 30, range=hist_range[3+i*9])  
                plt.subplot(335)
                n, bins, patches = plt.hist(diff[:,4+i*9], 30, range=hist_range[4+i*9])  
                plt.subplot(336)
                n, bins, patches = plt.hist(diff[:,5+i*9], 30, range=hist_range[5+i*9])  
                plt.subplot(337)
                n, bins, patches = plt.hist(diff[:,6+i*9], 30, range=hist_range[6+i*9])  
            elif spot_reminder == 8:
                plt.subplot(331)
                n, bins, patches = plt.hist(diff[:,0+i*9], 30, range=hist_range[0+i*9])
                plt.subplot(332)
                n, bins, patches = plt.hist(diff[:,1+i*9], 30, range=hist_range[1+i*9])
                plt.subplot(333)
                n, bins, patches = plt.hist(diff[:,2+i*9], 30, range=hist_range[2+i*9])   
                plt.subplot(334)
                n, bins, patches = plt.hist(diff[:,3+i*9], 30, range=hist_range[3+i*9])  
                plt.subplot(335)
                n, bins, patches = plt.hist(diff[:,4+i*9], 30, range=hist_range[4+i*9])  
                plt.subplot(336)
                n, bins, patches = plt.hist(diff[:,5+i*9], 30, range=hist_range[5+i*9])  
                plt.subplot(337)
                n, bins, patches = plt.hist(diff[:,6+i*9], 30, range=hist_range[6+i*9])  
                plt.subplot(338)
                n, bins, patches = plt.hist(diff[:,7+i*9], 30, range=hist_range[7+i*9])  
            else:
                plt.close()
                i=i-1
            #else:
            #    plt.subplot(331)
            #    plt.plot(result[0+i*9,:])
            #    plt.subplot(332)
            #    plt.plot(result[1+i*9,:])
            #    plt.subplot(333)
            #    plt.plot(result[2+i*9,:])
            #    plt.subplot(334)
            #    plt.plot(result[3+i*9,:])
            #    plt.subplot(335)
            #    plt.plot(result[4+i*9,:])
            #    plt.subplot(336)
            #    plt.plot(result[5+i*9,:])
            #    plt.subplot(337)
            #    plt.plot(result[6+i*9,:])
            #    plt.subplot(338)
            #    plt.plot(result[7+i*9,:])
            #    plt.subplot(339)
            #    plt.plot(result[8+i*9,:])
        else:
            plt.subplot(331)
            n, bins, patches = plt.hist(diff[:,0+i*9], 30, range=hist_range[0+i*9])
            plt.subplot(332)
            n, bins, patches = plt.hist(diff[:,1+i*9], 30, range=hist_range[1+i*9])
            plt.subplot(333)
            n, bins, patches = plt.hist(diff[:,2+i*9], 30, range=hist_range[2+i*9])   
            plt.subplot(334)
            n, bins, patches = plt.hist(diff[:,3+i*9], 30, range=hist_range[3+i*9])  
            plt.subplot(335)
            n, bins, patches = plt.hist(diff[:,4+i*9], 30, range=hist_range[4+i*9])  
            plt.subplot(336)
            n, bins, patches = plt.hist(diff[:,5+i*9], 30, range=hist_range[5+i*9])  
            plt.subplot(337)
            n, bins, patches = plt.hist(diff[:,6+i*9], 30, range=hist_range[6+i*9])  
            plt.subplot(338)
            n, bins, patches = plt.hist(diff[:,7+i*9], 30, range=hist_range[7+i*9])  
            plt.subplot(339)
            n, bins, patches = plt.hist(diff[:,8+i*9], 30, range=hist_range[8+i*9])  

def multihist2(Von, Voff):       
    npoint=len(Von[0,:])

    hist_range=[]  # list
    for i in range(npoint):
        maxbin=max(np.nanmax(Von[:,i]), np.nanmax(Voff[:,i]))
        minbin=min(np.nanmin(Von[:,i]), np.nanmin(Voff[:,i]))
        hist_range=(minbin, maxbin)
        plt.figure()
        plt.subplot(211)
        plt.title("Histogram of $I_{Von}$ and $I_{Voff}$, thresholded", fontsize=14, fontweight='bold')
        n1, bins1, patches = plt.hist(Von[:,i], 30, range=hist_range, color='b', normed=True, alpha=0.5, label='1st phase')  
        n2, bins2, patches = plt.hist(Voff[:,i], 30, range=hist_range, color='r', normed=True, alpha=0.5, label='2nd phase')  
        plt.subplot(212)
        binsize=bins1[1]-bins1[0]
        xnew=np.arange(minbin, maxbin-binsize, binsize)      
        plt.bar(bins1[:-1], n1-n2, width=bins1[1]-bins1[0], color='b', alpha=0.5)  
        plt.bar(bins1[:-1], n2-n1, width=bins1[1]-bins1[0], color='r', alpha=0.5)  
        #plt.plot(n1-n2)
        plt.title("$I_{Von}$ - $I_{Voff}$")  


        
        
     