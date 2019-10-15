#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 11:37:32 2019

@author: kiwitz
"""
import matplotlib as mpl
mpl.use('Agg')
import fenics as fcs
import numpy as np
import matplotlib.pyplot as plt
import json
from copy import deepcopy
import os
import pandas as pd
import seaborn as sns
import itertools
from scipy.constants import N_A

def sortDict(x,keys):
    res = x
    for i in keys:
        res = res[i]
    return res
def groupByKey(ls,keys):
    ls.sort(key=lambda x:sortDict(x,keys))
    l = []
    for i,e in itertools.groupby(ls,key=lambda x:sortDict(x,keys)):
        l.append(list(e))
    return l
def gradientPlot(path,xKey,**kwargs):
    
#    imgPath = path+"plots/"
#    if "imgPath" in kwargs:
#        imgPath = kwargs["imgPath"]+"/"+kwargs["l"]
#    os.makedirs(imgPath,exist_ok=True)
    
    with open(path+"gradient_dump.json","r") as file:
        dump = json.load(file)
        
#    sns.set_context("poster", font_scale=1.5,rc={
#            "lines.linewidth": 1,
#            "lines.markersize":20,
#            'xtick.labelsize': 'small',
#            'ytick.labelsize': 'small',
#            "xtick.major.width": 0.8,
#            "xtick.minor.width": 0.8,
#            "xtick.major.size": 3.5,
#            "xtick.minor.size": 3.5})
    sns.set_context("paper",font_scale=1.5,rc={
            "axes.labelsize":"large",
            "axes.titlesize":"large",
            "lines.markersize":25,
            'xtick.labelsize': 'large',
            'ytick.labelsize': 'large'}
            )
    dump = groupByKey(dump,["dict",xKey])
    for i,e in enumerate(dump):
        dump[i] = groupByKey(e,["field"])

    x = []
    y = []
    frames = []
    for singleScan in dump:
        
        
        for field_scan in singleScan:
            if field_scan[0]["field"] == "il2":
#            if True:
                gradients = [i["gradient"]/19.584 for i in field_scan]#volume quickfix!!
                con = [i["concentration"]/19.584 for i in field_scan]
                x_s = round(field_scan[0]["dict"][xKey]/(N_A**-1*10e9))
                x_tot = round((field_scan[0]["dict"]["R_il2_s"]+field_scan[0]["dict"]["R_il2_f"])/(N_A**-1*10e9))
                x = round(x_s/x_tot*100)
#                df = pd.DataFrame(
#                        {"x":[x]*len(gradients),"v":gradients,"type":"gradient","Cytokine":field_scan[0]["field"]},
#                        )
#                df = df.append(pd.DataFrame(
#                        {"x":[x]*len(gradients),"v":con,"type":"concentration","Cytokine":field_scan[0]["field"]}
#                        ))
                df = pd.DataFrame(
                        {"x":[x]*len(gradients),"v":gradients,"type":"gradient"},
                        )
                df = df.append(pd.DataFrame(
                        {"x":[x]*len(gradients),"v":con,"type":"concentration"}
                        ))
                
                frames.append(df)
    frame = frames[0]
    for i in frames[1:]:
        frame = frame.append(i)
        
    fig = plt.figure(figsize=(6.4,4))
    ax1 = plt.gca()
    color = "tab:red"
    sns.lineplot(x="x",y="v",ci="sd",ax=ax1,data=frame.where(frame["type"] == "gradient"),color=color,dashes={"gradient":(1,0)},markers={"gradient":"."},style="type",legend=False)
    ax1.set_ylabel(r'$\nabla$[IL-2] $\frac{nM}{dm}$',color=color)
    ax2 = ax1.twinx()
    color = "tab:blue"
    sns.lineplot(x="x",y="v",ci="sd",ax=ax2,data=frame.where(frame["type"] == "concentration"),color=color,dashes={"concentration":(1,0)},markers={"concentration":"."},style="type",legend=False)
    ax2.set_ylabel(r'[IL-2]  nM ',color=color)
    plt.tight_layout()
    ax1.set_xlabel(r'% of $R_{high}$ on secretors')
    
    ax1.set_ylim((0,ax1.get_ylim()[1]))
    ax2.set_ylim((0,ax2.get_ylim()[1]))
    os.makedirs(kwargs["plotPath"],exist_ok=True)
    if "plotPath" in kwargs:
        plt.savefig(kwargs["plotPath"]+"gradient.pdf")
    return fig
        
    
        
    
y = []
path = "/extra/kiwitz/results_parameter_scan/"
gradientPlot(path,"R_il2_s",plotPath="/home/kiwitz/parameter_scan_plots/")

#path = "/home/lukas/"
#gradientPlot(path,plotPath="/home/lukas/")

