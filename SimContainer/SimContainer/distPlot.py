#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 12:45:06 2019

@author: kiwitz
"""

import fenics as fcs
import numpy as np
import matplotlib.pyplot as plt
import json
from copy import deepcopy
import os
import pandas as pd
import seaborn as sns

def makeDataFrame_1(x,runs,**kwargs):
    frames = []
    l = kwargs["label"] if "label" in kwargs else "y"
    for run in runs:
        for i,row in enumerate(run):
                frames.append(pd.DataFrame({"x":[x[i]]*len(row),"u":row,"l":l}))
    frame = frames[0]
    for i in frames[:-1]:
        frame = frame.append(i)
    if "cutoff" in kwargs:
        frame = frame[frame.u  < kwargs["cutoff"]]
    return frame

def makeDataFrame_2(x,rows,**kwargs):
    frames = []
    l = kwargs["label"] if "label" in kwargs else "y"
    for i,e in enumerate(rows):
                frames.append(pd.DataFrame({"x":x[i]*len(e),"u":e,"l":l}))
    frame = frames[0]
    for i in frames[:-1]:
        frame = frame.append(i)
#    frame = frame[frame.u < 10]
    return frame

def myPlot(data,cellFunction,distFunction,runFunction,resultFunction):
    result = []
    for i,run in enumerate(data):
        runResult = []
        for o,dist in enumerate(run):
            distResult = []
            for n,cell in enumerate(dist):
                cellResult = cellFunction(n,cell)
                distResult.append(cellResult)
            distResult = distFunction(o,distResult)
            runResult.append(distResult)
        runResult = runFunction(i,runResult)
        result.append(runResult)
    return resultFunction(0,result)
    
def sanitize(dump):
    data = []
    x = []
    for dist in list(dump.values())[0][0]:
        cell = dist[0]
        x.append(cell["x"])
        del(cell)
        del(dist)
     
    
    for run in dump["il2"]:
        resultRun = []
        for dist in run:
            v = [{"il2":cell["v"]} for cell in dist]
            resultRun.append(v)
        data.append(resultRun)  
    
    il6Run = dump["il6"][0]
    for run in data:
        for i,dist in enumerate(run):
            for o,cell in enumerate(dist):
                cell["il6"] = il6Run[i][o]["v"]
                
    
    #sanitize list
    for run in data:
            for o,dist in enumerate(run):
                    indices = [dist.index(cell) for cell in dist if cell["il2"] < 10e-10 or cell["il6"] < 10e-10] 
                    run[o] = list(np.delete(dist,indices))
    return x,data
path = "./"
def distPlot(path,thresholds,**kwargs):
    imgPath = "plots/"
    os.makedirs(path+imgPath,exist_ok=True)
    
    with open(path+"distPlot_dump_il2.json","r") as file:
        dump_il2 = json.load(file)
    with open(path+"distPlot_dump_il6.json","r") as file:
        dump_il6 = json.load(file)  
    dump = {"il2":dump_il2,"il6":dump_il6}

    x,data = sanitize(dump)
#    x = np.arange(-1.5,1.5,0.2)
    il2Average = myPlot(
            data,
            lambda i,l : l["il2"], #cell
            lambda i,l: l, # for each distance value
            lambda i,l: l, # for single run
            lambda i,l: [row for row in np.transpose(l)] # for all runs
            )
    il2SD = myPlot(
            data,
            lambda i,l : l["il2"],
            lambda i,l: np.mean(l),
            lambda i,l: l,
            lambda i,l: [np.std(row) for row in np.transpose(l)]
            )
    il6Average = myPlot(
            data,
            lambda i,l : l["il6"],
            lambda i,l: l,
            lambda i,l: l,
            lambda i,l: [row for row in np.transpose(l)]
            )
    
    def treshhold(l,f):
        indices = [l.index(e) for e in l if not f(e)]
        return list(np.delete(l,indices))
    #t = [0.0025,0.015]
    t = thresholds
    il2p = myPlot(
            data,
            lambda i,l : l,
            lambda i,l: treshhold(l,lambda x: x["il2"] > t[0]),
            lambda i,l: [len(e) for e in l],
            lambda i,l: [e for e in np.transpose(l)]
            )
    il2pSD = myPlot(
            data,
            lambda i,l : l,
            lambda i,l: treshhold(l,lambda x: x["il2"] > t[0]),
            lambda i,l: [len(e) for e in l],
            lambda i,l: [np.std(e) for e in np.transpose(l)]
            )
    il6p = myPlot(
            data,
            lambda i,l : l,
            lambda i,l: treshhold(l,lambda x: x["il6"] > t[1]),
            lambda i,l: [len(e) for e in l],
            lambda i,l: [e for e in np.transpose(l)]
            )
    
    il2pil6p = myPlot(
            data,
            lambda i,l : l,
            lambda i,l: treshhold(l,lambda x: x["il2"] > t[0] and x["il6"] > t[1]),
            lambda i,l: [len(e) for e in l],
            lambda i,l: [e for e in np.transpose(l)]
            )
    il2pil6n = myPlot(
            data,
            lambda i,l : l,
            lambda i,l: treshhold(l,lambda x: x["il2"] > t[0] and x["il6"] < t[1]),
            lambda i,l: [len(e) for e in l],
            lambda i,l: [e for e in np.transpose(l)]
            )
    il2nil6p = myPlot(
            data,
            lambda i,l : l,
            lambda i,l: treshhold(l,lambda x: x["il2"] < t[0] and x["il6"] > t[1]),
            lambda i,l: [len(e) for e in l],
            lambda i,l: [e for e in np.transpose(l)]
            )
    il2nil6n = myPlot(
            data,
            lambda i,l : l,
            lambda i,l: treshhold(l,lambda x: x["il2"] < t[0] and x["il6"] < t[1]),
            lambda i,l: [len(e) for e in l],
            lambda i,l: [e for e in np.transpose(l)]
            )
    
    #conditionals
    il2_il6p = myPlot(
            data,
            lambda i,l : l,
            lambda i,l: [e["il2"] for e in treshhold(l,lambda x: x["il6"] > t[1])],
            lambda i,l: [e for e in l],
            lambda i,l: [e for e in np.transpose(l)]
            )
    il2_il6n = myPlot(
            data,
            lambda i,l : l,
            lambda i,l: [e["il2"] for e in treshhold(l,lambda x: x["il6"] < t[1])],
            lambda i,l: [e for e in l],
            lambda i,l: [e for e in np.transpose(l)]
            )
    il6_il2p = myPlot(
            data,
            lambda i,l : l,
            lambda i,l: [e["il6"] for e in treshhold(l,lambda x: x["il2"] > t[0])],
            lambda i,l: [e for e in l],
            lambda i,l: [e for e in np.transpose(l)]
            )
    il6_il2n = myPlot(
            data,
            lambda i,l : l,
            lambda i,l: [e["il6"] for e in treshhold(l,lambda x: x["il2"] < t[0])],
            lambda i,l: [e for e in l],
            lambda i,l: [e for e in np.transpose(l)]
            )
    
#    xl = r'distance from left boundary $[\mu m]$'
    
    il2Average = makeDataFrame_1(x,il2Average,label="IL-2",cutoff=kwargs["cutoff"][0])
    il6Average = makeDataFrame_1(x,il6Average,label="IL-6",cutoff=kwargs["cutoff"][1])
    
    il2p = makeDataFrame_2(x,il2p,label="IL-2")
    il6p = makeDataFrame_2(x,il6p,label="IL-6")
    simpleFractions = il2p.append(il6p)
#    
    il2pil6p = makeDataFrame_2(x,il2pil6p,label="IL2$^+$ IL6$^+$")
    il2pil6n = makeDataFrame_2(x,il2pil6n,label="IL2$^+$ IL6$^-$")
    il2nil6p = makeDataFrame_2(x,il2nil6p,label="IL2$^-$ IL6$^+$")
    il2nil6n = makeDataFrame_2(x,il2nil6n,label="IL2$^-$ IL6$^-$")
#    
    conditionalFractions = il2pil6p
    for i in [il2pil6n,il2nil6p,il2nil6n]:
        conditionalFractions = conditionalFractions.append(i)
        
#    il2_il6p = makeDataFrame_2(x,il2_il6p)
#    il2_il6n = makeDataFrame_2(x,il2_il6n)
#    il6_il2p = makeDataFrame_2(x,il6_il2p)
#    il6_il2n = makeDataFrame_2(x,il6_il2n)
#    
#    conditionalConcentrations= il2_il6p
#    for i in [il2_il6n,il6_il2p,il6_il2n]:
#        conditionalConcentrations = conditionalConcentrations.append(i)
    
    
#    sns.lineplot(x="x",y="u",ci="sd",data=il2p)
#    sns.lineplot(x="x",y="u",ci="sd",data=il6p)
    
    sns.set_context("paper", rc={"lines.linewidth": 0.1,"lines.markersize":10})
    
    xTicks = [i for i in np.arange(-1.5,2,0.5)]
    xScale = [round(i*100) for i in np.arange(0.15,3.5,0.5)]
    plt.figure(1)
    color = "tab:red"
    fig, ax1 = plt.subplots()
    sns.lineplot(x="x",y="u",ci="sd",data=il2Average,ax=ax1,color=color,markers={"IL-2":"."},style="l",legend=False)
    ax1.set_ylim(0,0.1)
    ax1.tick_params(axis="y",labelcolor=color)
    ax1.set_ylabel(r'[IL-2][nM]',color=color)
    ax2 = ax1.twinx()
    color = "tab:blue"
    sns.lineplot(x="x",y="u",ci="sd",data=il6Average,ax=ax2,color=color,markers=["."],style="l",legend=False)
    ax2.tick_params(axis="y",labelcolor=color)
    ax2.set_ylabel(r'[IL-6][nM]',color=color)
    ax2.set_ylim(0,1)
    ax1.set_xlabel(r'distance from left boundary $[\mu m]$')
    plt.xticks(xTicks,xScale)
    
    plt.savefig(path+imgPath+"averageConcentrations.png",dpi=600)
#    plt.close()
   
#    plt.figure(2)
#    sns.lineplot(x="x",y="u",ci="sd",data=simpleFractions,hue="l")
#
#    plt.legend()
#    plt.xticks(np.arange(-1,1.25,0.25),np.arange(15,225,25))
#    plt.xlabel(r'distance from left boundary $[\mu m]$')
#    plt.ylabel(r'fraction of positive cells')
#    plt.savefig(path+imgPath+"simpleFractions.png",dpi=600)
#   
#    plt.figure(3)
#    sns.lineplot(x="x",y="u",ci="sd",data=conditionalFractions,hue="l")
#    plt.xticks(np.arange(-1,1.25,0.25),np.arange(15,225,25))
#    plt.xlabel(r'distance from left boundary $[\mu m]$')
#    plt.ylabel(r'fraction of positive cells')
#   
#    plt.savefig(path+imgPath+"conditionalFractions.png",dpi=600)
##    
#    plt.figure(4)
#    color = "tab:red"
#    fig, ax1 = plt.subplots()
#    sns.lineplot(x="x",y="u",ci="sd",data=il2_il6p,ax=ax1,color=color)
#    sns.lineplot(x="x",y="u",ci="sd",data=il2_il6n,ax=ax1,color=color)
##    ax1.plot(x,il2_il6p,"-",label="<IL2>|IL6$^+$",color=color)
##    ax1.plot(x,il2_il6n,"--",label="<IL2>|IL6$^-$",color=color)
#    ax1.tick_params(axis="y",labelcolor=color)
#    ax1.set_ylabel(r'IL2 nM',color=color)
#    ax1.legend(loc="upper right")
#    ax2 = ax1.twinx()
#    
#    color = "tab:blue"
#    sns.lineplot(x="x",y="u",ci="sd",data=il6_il2p,ax=ax2,color=color)
#    sns.lineplot(x="x",y="u",ci="sd",data=il6_il2n,ax=ax2,color=color)
##    ax2.plot(x,il6_il2p,"-",label="<IL6>|IL2$^+$",color=color)
##    ax2.plot(x,il6_il2n,"--",label="<IL6>|IL2$^-$",color=color)
#    
#    
#    ax2.tick_params(axis="y",labelcolor=color)
#    ax2.set_ylabel(r'IL6 nM',color=color)
#    ax2.legend(loc="lower right")
#    ax1.set_xlabel(r'distance from left boundary $[\mu m]$')
#    plt.xticks(np.arange(-1,1.25,0.25),np.arange(15,225,25))
##    plt.savefig(path+imgPath+"conditionalConcentrations.png",dpi=600)
