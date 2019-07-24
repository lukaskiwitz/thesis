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


result = []
tmpIndexList = []
for file in os.listdir("./sol/distplot"):
    
    il2 = fcs.Function(V)
    with fcs.HDF5File(fcs.MPI.comm_world,"./sol/distplot/"+file,"r") as f:
        f.read(il2,"il2")
    
    
    ds = fcs.Measure("ds", domain=mesh, subdomain_data=boundary_markers)
    area = 4*np.pi*5**2
    for i in cell_data:   
        print("patch no: "+str(i["patch"]))
        r = i["center"][0]
        i["v"] = fcs.assemble(il2*ds(i["patch"]))/(4*np.pi*0.05**2)*10**9
        if r in tmpIndexList:
            result[tmpIndexList.index(r)].append(i)
        else:
            result.append([i])
            tmpIndexList.append(r)
    
x = []
d = []
for i in result:
    v= [e["v"] for e in i]
    indecies = [i for (i,v) in enumerate(v) if v < 0]
    v = np.delete(v,indecies)
#    print(indecies)
    x.append(i[0]["center"][0])
    d.append(v)

#d = np.array()
#d = d[2:]

#d = np.transpose(d,axes=(1,0))
for i,v in enumerate(d):
#    plt.plot(np.ones(len(v))*x[i],v,"+")
    pass
plt.plot(x,[len(i) for i in d])
#plt.ylim((np.min(d[0]),np.max(d[-1])))

    


