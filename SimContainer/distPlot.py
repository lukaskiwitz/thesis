#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 12:45:06 2019

@author: kiwitz
"""

import fenics as fcs
import numpy as np
import matplotlib.pyplot as plt
import h5py
import json
import os
import re
from copy import deepcopy

mesh = fcs.Mesh()
with fcs.XDMFFile("./cache/meshCache_il2.xdmf") as f:
    f.read(mesh)

V = fcs.FunctionSpace(mesh,"P",1)
    
subPath = "./cache/boundary_markers.h5"
boundary_markers = fcs.MeshFunction("size_t",mesh,mesh.topology().dim() - 1)
with fcs.HDF5File(fcs.MPI.comm_world,subPath,"r") as f:
    f.read(boundary_markers,"/boundaries")


cell_dump_path = "./cell_dump.json"
with open(cell_dump_path,"r") as file:
    d = file.read()
    cell_data = json.loads(d)


results = {}
fields = ["il6"]
for field in fields:
    
    fileList = os.listdir("./sol/distplot")
    bList = [((re.compile(field)).search("{field}(?=.*\.h5)".format(field=i))) for i in fileList]
    fileList = np.delete(fileList,[i for i,e in enumerate(bList) if not e])
    fileList = fileList[-1:]
    tmpIndexList = []
    result = []
    for no,file in enumerate(fileList):
        
        u = fcs.Function(V)
        with fcs.HDF5File(fcs.MPI.comm_world,"./sol/distplot/"+file,"r") as f:
            f.read(u,field)
        print("reading file {file} ({no}/{tot})".format(file=file,no=no,tot=len(fileList)))
        
        ds = fcs.Measure("ds", domain=mesh, subdomain_data=boundary_markers)
        area = 4*np.pi*5**2
        for i in cell_data:   
#            print("patch no: "+str(i["patch"]))
            
            r = i["center"][0]
            i["v"] = fcs.assemble(u*ds(i["patch"]))/(4*np.pi*0.05**2)*10**9
            if r in tmpIndexList:
                result[tmpIndexList.index(r)].append(i)
            else:
                result.append([i])
                tmpIndexList.append(r)
    results[field] = deepcopy(result)
    
color = ["r","b"]
thresh = [10e-10,10e-10]

#with open("./distPlot_dump.json","w") as file:
#    file.write(json.dumps(results))
with open("./distPlot_dump.json","r") as file:
    results_old = json.load(file)
results.update(results_old)
    
for c,field in enumerate(["il2","il6"]):
    x = []
    d = []
    if field in results:
        for i in results[field]:
            v= [e["v"] for e in i]
            indecies = [i for (i,v) in enumerate(v) if v < thresh[c]]
            v = np.delete(v,indecies)
        #    print(indecies)
            x.append(i[0]["center"][0])
            d.append(v)
    
    #d = np.array()
    #d = d[2:]
    
    #d = np.transpose(d,axes=(1,0))
    y = [np.mean(v) for i,v in enumerate(d)]
    plt.plot(x,y,color[c]+"+",label=field)

#    plt.plot(x,[len(i) for i in d])
#    plt.ylim((np.min(d.flatten()),np.max(d.flatten())))
plt.savefig("./distPlot.png")
plt.xticks(np.arange(-1,1.25,0.25),np.arange(0,2.25,0.25)*100)
plt.xlabel(r'distance from left boundary $\left[ \mu \operatorname{m} \right]$')
plt.ylabel(r'nM')
plt.legend()
    


