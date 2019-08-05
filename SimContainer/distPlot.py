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

with open("./distPlot_dump.json","r") as file:
    dump = json.load(file)    
x,data = sanitize(dump)

il2Average = myPlot(
        data,
        lambda i,l : l["il2"], #cell
        lambda i,l: np.mean(l), # for each distance value
        lambda i,l: l, # for single run
        lambda i,l: [np.mean(row) for row in np.transpose(l)] # for all runs
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
        lambda i,l: np.mean(l),
        lambda i,l: l,
        lambda i,l: [np.mean(row) for row in np.transpose(l)]
        )

def treshhold(l,f):
    indices = [l.index(e) for e in l if not f(e)]
    return list(np.delete(l,indices))
t = [0.64,0.54]

il2p = myPlot(
        data,
        lambda i,l : l,
        lambda i,l: treshhold(l,lambda x: x["il2"] > t[0]),
        lambda i,l: [len(e) for e in l],
        lambda i,l: [np.mean(e) for e in np.transpose(l)]
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
        lambda i,l: [np.mean(e) for e in np.transpose(l)]
        )

il2pil6p = myPlot(
        data,
        lambda i,l : l,
        lambda i,l: treshhold(l,lambda x: x["il2"] > t[0] and x["il6"] > t[1]),
        lambda i,l: [len(e) for e in l],
        lambda i,l: [np.mean(e) for e in np.transpose(l)]
        )
il2pil6n = myPlot(
        data,
        lambda i,l : l,
        lambda i,l: treshhold(l,lambda x: x["il2"] > t[0] and x["il6"] < t[1]),
        lambda i,l: [len(e) for e in l],
        lambda i,l: [np.mean(e) for e in np.transpose(l)]
        )
il2nil6p = myPlot(
        data,
        lambda i,l : l,
        lambda i,l: treshhold(l,lambda x: x["il2"] < t[0] and x["il6"] > t[1]),
        lambda i,l: [len(e) for e in l],
        lambda i,l: [np.mean(e) for e in np.transpose(l)]
        )
il2nil6n = myPlot(
        data,
        lambda i,l : l,
        lambda i,l: treshhold(l,lambda x: x["il2"] < t[0] and x["il6"] < t[1]),
        lambda i,l: [len(e) for e in l],
        lambda i,l: [np.mean(e) for e in np.transpose(l)]
        )

#conditionals
il2_il6p = myPlot(
        data,
        lambda i,l : l,
        lambda i,l: [e["il2"] for e in treshhold(l,lambda x: x["il6"] > t[1])],
        lambda i,l: [np.mean(e) for e in l],
        lambda i,l: [np.mean(e) for e in np.transpose(l)]
        )
il2_il6n = myPlot(
        data,
        lambda i,l : l,
        lambda i,l: [e["il2"] for e in treshhold(l,lambda x: x["il6"] < t[1])],
        lambda i,l: [np.mean(e) for e in l],
        lambda i,l: [np.mean(e) for e in np.transpose(l)]
        )
il6_il2p = myPlot(
        data,
        lambda i,l : l,
        lambda i,l: [e["il6"] for e in treshhold(l,lambda x: x["il2"] > t[0])],
        lambda i,l: [np.mean(e) for e in l],
        lambda i,l: [np.mean(e) for e in np.transpose(l)]
        )
il6_il2n = myPlot(
        data,
        lambda i,l : l,
        lambda i,l: [e["il6"] for e in treshhold(l,lambda x: x["il2"] < t[0])],
        lambda i,l: [np.mean(e) for e in l],
        lambda i,l: [np.mean(e) for e in np.transpose(l)]
        )


plt.figure(1)
color = "tab:red"
fig, ax1 = plt.subplots()
ax1.errorbar(x,il2Average,yerr=il2SD,capsize=5,color=color)
ax1.tick_params(axis="y",labelcolor=color)
ax1.set_ylabel(r'IL2 nM',color=color)
ax2 = ax1.twinx()
color = "tab:blue"
ax2.plot(x,il6Average,color=color)
ax2.tick_params(axis="y",labelcolor=color)
ax2.set_ylabel(r'IL6 nM',color=color)

ax1.set_xlabel(r'distance from left boundary $[\mu m]$')
plt.xticks(np.arange(-1,1.25,0.25),np.arange(0,225,25))
plt.savefig("averageConcentrations.png",dpi=600)

plt.close()
plt.figure(2)
plt.errorbar(x,il2p,yerr=il2pSD,label="il2",capsize=5)
plt.plot(x,il6p,label="il6")
plt.legend()
plt.xticks(np.arange(-1,1.25,0.25),np.arange(0,225,25))
plt.xlabel(r'distance from left boundary $[\mu m]$')
plt.ylabel(r'fraction of positive cells')
plt.savefig("simpleFractions.png",dpi=600)

plt.figure(3)
plt.plot(x,il2pil6p,label="IL2$^+$ IL6$^+$")
plt.plot(x,il2pil6n,label="IL2$^+$ IL6$^-$")
plt.plot(x,il2nil6p,label="IL2$^-$ IL6$^+$")
plt.plot(x,il2nil6n,label="IL2$^-$ IL6$^-$")
plt.legend(loc= "upper right")
plt.xticks(np.arange(-1,1.25,0.25),np.arange(0,225,25))
plt.xlabel(r'distance from left boundary $[\mu m]$')
plt.ylabel(r'fraction of positive cells')

plt.savefig("conditionalFractions.png",dpi=600)

plt.figure(4)
color = "tab:red"
fig, ax1 = plt.subplots()
ax1.plot(x,il2_il6p,"-",label="<IL2>|IL6$^+$",color=color)
ax1.plot(x,il2_il6n,"--",label="<IL2>|IL6$^-$",color=color)
ax1.tick_params(axis="y",labelcolor=color)
ax1.set_ylabel(r'IL2 nM',color=color)
ax1.legend(loc="upper right")
ax2 = ax1.twinx()

color = "tab:blue"
ax2.plot(x,il6_il2p,"-",label="<IL6>|IL2$^+$",color=color)
ax2.plot(x,il6_il2n,"--",label="<IL6>|IL2$^-$",color=color)


ax2.tick_params(axis="y",labelcolor=color)
ax2.set_ylabel(r'IL6 nM',color=color)
ax2.legend(loc="lower right")
ax1.set_xlabel(r'distance from left boundary $[\mu m]$')
plt.xticks(np.arange(-1,1.25,0.25),np.arange(0,225,25))
plt.savefig("conditionalConcentrations.png",dpi=600)
