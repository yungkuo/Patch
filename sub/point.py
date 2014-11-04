# -*- coding: utf-8 -*-
"""
Created on Sun Sep 21 20:47:36 2014

@author: Philip
"""
import numpy as np
import matplotlib.pyplot as plt
from pylab import * 

scan = 2

if __name__ == "__main__": 
    print("Your scan pixel is +-%d" %scan)    
    
def pIO(mov):
    print("Return to Figure 1, and click points of interests")
    print("If finished, press center button of your mouse")
    nrow=len(mov[1,:,1])
    ncol=len(mov[1,1,:])
    pts = ginput(0,0) 
    pts=np.array(pts)
    col_pts=np.around(pts[:,0])
    row_pts=np.around(pts[:,1])
    plot(col_pts,row_pts, 'ro')
    plt.xlim([0,ncol])
    plt.ylim([nrow,0])
    pts=pts.astype(int)
    pts_rc=zip(row_pts, col_pts)
    
    return pts_rc
    
def localmax(refimg, pts):
    
    drow_dcol= np.zeros((len(pts[:,0]),2))
    for i in range(len(pts[:,0])):
        local=mask(refimg, pts[i,:])
        #mask = refimg[pts[i,0]-scan:pts[i,0]+scan+1,pts[i,1]-scan:pts[i,1]+scan+1]
        drow_dcol[i,:]=np.array(zip(*np.where(local==local.max())))   # * unpack the tuple, return its value as an input element of zip)
    drow_dcol=drow_dcol.astype(int)-scan
    
    pts_new=pts+drow_dcol
    
    return pts_new
    
def mask(refimg, pts):
    local=refimg[pts[0]-scan:pts[0]+scan+1,pts[1]-scan:pts[1]+scan+1]
    return local