#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 12:09:23 2019

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
fields = ["il2","il6"]
for field in fields:
    
    fileList = os.listdir("./sol/distplot")
    bList = [((re.compile(field)).search("{field}(?=.*\.h5)".format(field=i))) for i in fileList]
    fileList = np.delete(fileList,[i for i,e in enumerate(bList) if not e])
#    fileList = fileList[-2:]
    for no,file in enumerate(fileList):
        tmpIndexList = []
        result = []
        u = fcs.Function(V)
        with fcs.HDF5File(fcs.MPI.comm_world,"./sol/distplot/"+file,"r") as f:
            f.read(u,field)
        print("reading file {file} ({no}/{tot})".format(file=file,no=no,tot=len(fileList)))
        
        ds = fcs.Measure("ds", domain=mesh, subdomain_data=boundary_markers)
        for i in cell_data:   
#            print("patch no: "+str(i["patch"]))
            r = i["center"][0]
            v = (fcs.assemble(u*ds(i["patch"]))/(4*np.pi*0.05**2)*10**9)
            if r in tmpIndexList:
                result[tmpIndexList.index(r)].append({"x":r,"v":v})
            else:
                result.append([{"x":r,"v":v}])
                tmpIndexList.append(r)
        if not field in results:
            results[field] = [deepcopy(result)]
        else:
            results[field].append(deepcopy(result))

with open("./distPlot_dump.json","w") as file:
    file.write(json.dumps(results))