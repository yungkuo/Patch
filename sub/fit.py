import numpy as np
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
from pylab import *  


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

def resEvenGauss(p, y, x):
    A, mu, sd = p
    y_fit = gaussian(A, x, mu, sd) 
    err = y - y_fit
    return err
    
def resOddGauss(p, y, x):
    A, mu, sd = p
    y_fit = oddGaussian(A, x, mu, sd)
    err = y - y_fit
    return err



def EOGdecomposition(y, x):  # decompose with even and odd gaussian y : nbins x: bins

    xnew=np.zeros((len(x)-1))
    for i in range(len(x)-1):
        xnew[i]= (x[i]+x[i+1])/2
    mu = sum(xnew*y)/sum(y)           
    mu=0  
    sigma = sqrt(sum(y*(xnew-mu)**2)/sum(y))
    A1 = y.max()
    A2 = 0
    p = [A1, A2, mu, mu, sigma, sigma]   # initial guess parameters
    y_init=gaussian(A1, xnew, mu, sigma)  + oddGaussian(A2, xnew, mu, sigma)
    
    plsq=leastsq(resEOG, p, args = (y, xnew))
    
    y_est = gaussian(plsq[0][0], xnew, plsq[0][2], plsq[0][4]) + oddGaussian(plsq[0][1], xnew, plsq[0][2] + plsq[0][3], plsq[0][5])
    y_even = gaussian(plsq[0][0], xnew, plsq[0][2], plsq[0][4]) 
    y_odd = oddGaussian(plsq[0][1], xnew, plsq[0][2] + plsq[0][3], plsq[0][5])
    if plsq[0][1] > 0:
        oe_ratio=y_odd.max()/y_even.max()
    else:
        oe_ratio=-1*y_odd.max()/y_even.max()
    
    f, ax = subplots()
    plt.plot(xnew, y, label='Real Data')
    plt.plot(xnew, y_init, 'r.', label='Starting Guess')
    plt.plot(xnew, y_est, 'g.', label='Fitted')
    plt.plot(xnew, y_even, 'c', label='even Gaussian')
    plt.plot(xnew, y_odd, 'm', label='odd Gaussian')
    plt.title('Histogram of $\Delta I/I_{Max}$, thresholded & Averaged over period', fontsize=14, fontweight='bold')
    plt.text(-0.1, -0.05, 'Odd /Even Gaussian ratio is %f ' % oe_ratio)
    plt.legend()
    plt.show()
    return  y_even, y_odd, oe_ratio
    
def fitGauss(y, x):  # decompose with even and odd gaussian

    xnew=np.zeros((len(x)-1))
    for i in range(len(x)-1):
        xnew[i]= (x[i]+x[i+1])/2
    mu = sum(xnew*y)/sum(y)           #note this correction
    sigma = sqrt(sum(y*(xnew-mu)**2)/sum(y))
    A = y.max()
    p = [A, mu, sigma]   # initial guess parameters
    y_init=gaussian(A, xnew, mu, sigma)
    
    plsq=leastsq(resEvenGauss, p, args = (y, xnew))
    
    y_est = gaussian(plsq[0][0], xnew, plsq[0][1], plsq[0][2]) 
    
    f, ax = subplots()
    plt.plot(xnew, y, label='Real Data')
    plt.plot(xnew, y_init, 'r.', label='Starting Guess')
    plt.plot(xnew, y_est, 'g.', label='Fitted')
    plt.legend()
    plt.show()
    return y_est
    
def fitOddGauss(y, x):  # decompose with even and odd gaussian

    xnew=np.zeros((len(x)-1))
    for i in range(len(x)-1):
        xnew[i]= (x[i]+x[i+1])/2
    mu = 0        #note this correction
    sigma = sqrt(sum(y*(xnew-mu)**2)/sum(abs(y)))
    A = y.max()
    p = [A, mu, sigma]   # initial guess parameters
    y_init=oddGaussian(A, xnew, mu, sigma)
    
    plsq=leastsq(resOddGauss, p, args = (y, xnew))
    
    y_est = oddGaussian(plsq[0][0], xnew, plsq[0][1], plsq[0][2]) 
    
    f, ax = subplots()
    plt.plot(xnew, y, label='Real Data')
    plt.plot(xnew, y_init, 'r.', label='Starting Guess')
    plt.plot(xnew, y_est, 'g.', label='Fitted')
    plt.legend()
    plt.show()
    return y_est
##################