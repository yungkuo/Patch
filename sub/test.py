# -*- coding: utf-8 -*-
"""
Created on Tue Nov 11 09:55:24 2014

@author: Philip
"""
ncycle=frame/period
curve=spot_intensity[:,0]
curveRS=spot_intensity[:,0].reshape(ncycle, period)
cv_up=np.mean(curveRS[:, :period/2], axis=1)
cv_down=np.mean(curveRS[:, period/2:period], axis=1)


cv_norm=(cv_up+cv_down)/2
cv_diff=cv_up-cv_down

up=0
down=0
uplist=[]
downlist=[]
for i in range(len(cv_norm)):
    if cv_norm[i] > threshold[0]:
        if cv_diff[i] > 0:
            up=up+1
            if down>0:
                downlist.append((down, i))
            down=0
        else:
            down=down+1
            if up>0:
                uplist.append((up, i))
            up=0
    else:
        if up > 0:
            uplist.append((up, i))
        elif down > 0:
            downlist.append((down, i))
        up=0
        down=0
        
Maxcorr1=sorted(uplist, reverse=True)[0]
Maxcorr2=sorted(downlist, reverse=True)[0]
Phase1upArray= curve[(Maxcorr1[1]-Maxcorr1[0])*period : (Maxcorr1[1])*period]
Phase2upArray= curve[(Maxcorr2[1]-Maxcorr2[0])*period : (Maxcorr2[1])*period]
