# -*- coding: utf-8 -*-
"""
Created on Thu Nov 06 18:22:10 2014

@author: KyoungWon
"""
selected_rsf=sth[0]
selected_diff=rth[0]
r=[]
s=[]
f=[]

for i in range(len(selected_rsf)):
    for j in range(len(selected_rsf[i])):
        if selected_rsf[i][j]==1:
            r.append(selected_diff[i][j])
        elif selected_rsf[i][j]==0:
            s.append(selected_diff[i][j])
        elif selected_rsf[i][j]==-1:
            f.append(selected_diff[i][j])