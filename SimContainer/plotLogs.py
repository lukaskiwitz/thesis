#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 12:18:02 2019

@author: lukas
"""

import matplotlib.pyplot as plt
import numpy as np
import json
import os
import pandas as pd
import decimal



#for q in [10]:
#    data = []
#    loc = "./logs_q{q}/".format(q=q)
#    for i,o in enumerate(os.listdir(loc)):
#            with open(loc+str(o),"r") as file:
#                d = file.read()
#                data.append((float(o),json.loads(d)))
#    data.sort()
#    
#    series = {}
#    for d in data:
#        for key in d[1].keys():
#            try:
#                if not key in series:
#                    series[key] = [d[1][key]["solver_0"]["flux"]]
#                else:
#                    series[key].append(d[1][key]["solver_0"]["flux"])
#            except:
#                pass
#            
#    x = np.linspace(1,4,20)
#    for e in series.values():
4#        
#    plt.title("$q_{center}="+str(q)+"$")
#    plt.xticks([2,4],["$10^2$","$10^4$"])
#    plt.xlabel("$R$")
#    plt.ylabel("flux")
#    plt.savefig("q{q}.png".format(q=q),dpi=400)
        

data = []
loc = "./logs/"
for i,o in enumerate(os.listdir(loc)):
        with open(loc+str(o),"r") as file:
            d = file.read()
            data.append((float(o),json.loads(d)))
data.sort()

series = {}
keyStr = "level"
for d in data:
    for key in d[1].keys():
        #if key == "entity_2":
         #   continue
        try:
            if not key in series:
                series[key] = [d[1][key]["solver_0"][keyStr]]
            else:
                series[key].append(d[1][key]["solver_0"][keyStr])
        except:
            pass
        

for e in series.values():
#    if max(e) > 0e-13:
    plt.plot(e)
plt.ylim(bottom=0)
#plt.title("$q=10$")
#plt.ylabel("flux")
#plt.savefig("q{q}.png".format(q=q),dpi=400)
