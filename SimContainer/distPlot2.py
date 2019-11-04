#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 11:19:45 2019

@author: kiwitz
"""

import fenics as fcs
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import json
from copy import deepcopy
import os
import pandas as pd
import seaborn as sns
import lxml.etree as ET
import itertools
from myDictSorting import groupByKey

sns.set_context("paper", font_scale=1.5,rc={
            "lines.linewidth": 1,
            "lines.markersize":20,
            'xtick.labelsize': 'small',
            'ytick.labelsize': 'small',
            "xtick.major.width": 0.8,
            "xtick.minor.width": 0.8,
            "xtick.major.size": 3.5,
            "xtick.minor.size": 3.5})
    
def _prepData(inTree):
    frames = []
    for scan in inTree.findall("/scan"):
        scanIndex = float(scan.get("i"))
        timeSteps = np.unique([int(i.get("i")) for i in scan.findall("./timeStep")])
        
        for t in timeSteps:#scan.findall("./timeStep"):
            
#            timeIndex = float(timeStep.get("i"))
#            for fieldName in [i.get("field") for i in timeStep.findall("./file")]:
            files = scan.findall("./timeStep[@i='{t}']/file".format(t=t))
            offset= 3 #one file per field
            cellResults = np.empty((1500,len(files)+offset))
            xTicks = []
            for cellIndex,cell in enumerate(files[0].findall("./cellResults/cell")):
#                patch = int(cell.find("./patch").text)
                x = json.loads(cell.find("./center").text)[0]
                cellResults[cellIndex,0] = x
                cellResults[cellIndex,1] = t
                cellResults[cellIndex,2] = scanIndex
                nameList = []
                for o,file in enumerate(files):
                    cell = file.findall("./cellResults/cell")[cellIndex]
                    fieldName = file.get("field")
                    nameList.append(fieldName)
                    cellResults[cellIndex,o+offset] = float(cell.find("./surfaceConcentration").text)
            cellFrame = pd.DataFrame(cellResults,columns=["x","time","scanIndex"]+nameList)              
            frames.append(cellFrame)
    #join dataframes
    result = frames[0]
    for i in range(len(frames)-1):
        result = result.append(frames[i+1])
    return result
                
                
# prep
path = "/extra/kiwitz/results_parameter_scan_qStrength /"
tree = ET.parse(path+"postProcess.xml")  
res = _prepData(tree)
 
scanGroups = res.groupby(by=["scanIndex"])
th = [0.1,2]
#result = []
keys = list(scanGroups.groups.keys())
d = scanGroups.get_group(keys[0])

img = "/home/kiwitz/parameter_scan_qStrength_plots/"
xTicks = np.unique(d["x"].to_numpy())
xScale = [0,"",2,"",5,"",6,"",8,"",10]
for key in keys:
    scanIndex = d.get("scanIndex").values[0]
    imgPath = img+str(scanIndex)
    os.makedirs(imgPath,exist_ok=True)
    
    d = scanGroups.get_group(key)
    #averages
        
    average = series = d.groupby(["x","time"],as_index=False).mean()
    plt.figure(1)
    color = "tab:red"
    fig, ax1 = plt.subplots()
    sns.lineplot(x="x",y="il2",ci="sd",data=average,ax=ax1,color=color,markers=["."],legend=False,style=True,err_style="bars")
    ax1.set_ylim(0)
    ax1.tick_params(axis="y",labelcolor=color)
    ax1.set_ylabel(r'IL-2 on cell surface (nM)',color=color)
    ax2 = ax1.twinx()
    color = "tab:blue"
    sns.lineplot(x="x",y="il6",ci="sd",data=average,ax=ax2,color=color,markers=["."],legend=False,style=True,err_style="bars")
    ax2.tick_params(axis="y",labelcolor=color)
    ax2.set_ylabel(r'IL-21 on cell surface (nM)',color=color)
    ax2.set_ylim(0)
    ax1.set_xlabel(r'distance from left boundary $[\mu m]$')
    ax1.set_xlabel(r'cell distance from boundary')
    plt.xticks(xTicks,xScale)
    plt.tight_layout()
    plt.savefig(imgPath+"/"+"averageConcentrations.pdf",dpi=1200)
    plt.close()
    
#     simpleFractions
    
    
    def makeSimpleFractions(key,t,label,d):
        
        df = d.loc[(d[key] > t)]
        series = df.groupby(["x","time"]).size()
        df= pd.DataFrame([series.get_values()],["n"]).T
        mi = series.index
        mi = mi.to_frame(index=False,name=["x","time"])
        
        df= df.join(mi)
        df= df.assign(l=label)
        return df
    
    il2p = makeSimpleFractions("il2",th[0],"IL-2",d)
    il6p = makeSimpleFractions("il6",th[1],"IL-6",d)
    simpleFractions = il2p.append(il6p)
    
    plt.figure(2)
    color = {"IL-2":"tab:red","IL-6":"tab:blue"}
    fig, ax1 = plt.subplots()
    sns.lineplot(x="x",y="n",ci="sd",data=simpleFractions,markers={"IL-2":".","IL-6":"."},legend=False,style="l",hue="l",dashes={"IL-2":[1,0],"IL-6":[1,0]})
    ax1.set_ylim(0)
    ax1.tick_params(axis="y")
    ax1.set_ylabel(r'Fraction of IL-2 positive cells (%)')
    ax1.set_xlabel(r'cell distance from boundary')
    plt.xticks(xTicks,xScale)
    plt.tight_layout()
    plt.savefig(imgPath+"/"+"simpleFractions.pdf",dpi=1200)
    plt.close()
    
    # conditional Fractions
    def makeConditionalFractions(keys,th,label,p,d):
        if not p[0] and not p[1]:
            df = d.loc[(d[keys[0]] < th[0]) & (d[keys[1]] < th[1])]
        if  p[0] and  p[1]:
            df = d.loc[(d[keys[0]] > th[0]) & (d[keys[1]] > th[1])]
        if not p[0] and  p[1]:
            df = d.loc[(d[keys[0]] < th[0]) & (d[keys[1]] > th[1])]
        if p[0] and not p[1]:
            df = d.loc[(d[keys[0]] > th[0]) & (d[keys[1]] < th[1])]
        series = df.groupby(["x","time"]).size()
        df= pd.DataFrame([series.get_values()],["n"]).T
        mi = series.index
        mi = mi.to_frame(index=False,name=["x","time"])
        
        df= df.join(mi)
        df= df.assign(l=label)
        return df
        
    il2pil6n = makeConditionalFractions(["il2","il6"],th,"IL2$^+$ IL21$^-$",(True,False),d)
    il2nil6p = makeConditionalFractions(["il2","il6"],th,"IL2$^-$ IL21$^+$",(False,True),d)
    il2pil6p = makeConditionalFractions(["il2","il6"],th,"IL2$^+$ IL21$^+$",(True,True),d)
    il2nil6n = makeConditionalFractions(["il2","il6"],th,"IL2$^-$ IL21$^-$",(False,False),d)
    
    conditionalFractions = il2pil6n
    for i in [il2nil6p,il2pil6p,il2nil6n]:
        conditionalFractions = conditionalFractions.append(i)
    
    plt.figure(3)
    sns.lineplot(x="x",y="n",ci="sd",data=conditionalFractions,hue="l",err_style="bars",marker=".")
    plt.xlabel(r'distance from left boundary $[\mu m]$')
    #plt.legend(loc="upper right")
    plt.xlabel(r'cell distance from boundary')
    plt.ylabel(r'fraction of positive cells (%)')
    plt.ylim(-5,105)
    plt.yticks(np.arange(0,120,20),np.arange(0,120,20))
    plt.xticks(xTicks,xScale)
    plt.tight_layout()
    plt.savefig(imgPath+"/"+"conditionalFractions.pdf",dpi=600)
    plt.close()
        





