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
from copy import copy,deepcopy
import mpi4py.MPI as MPI
import multiprocessing as mp
from math import ceil

def compute(no,h5File,field,extCache,**kwargs):
    mesh = kwargs["mesh"]
#    boundary_markers = kwargs["boundary_markers"]  
#    cell_data =  kwargs["cell_data"]
#    tmpIndexList = []
    
    u = kwargs["u"]
    V_vec = fcs.VectorFunctionSpace(mesh,"P",1)
    
#    ds = fcs.Measure("ds", domain=mesh, subdomain_data=boundary_markers)
    print("projecting gradient")
    grad = fcs.project(fcs.grad(u),V_vec,solver_type="gmres")
    result = fcs.assemble(fcs.sqrt(fcs.dot(grad,grad))*fcs.dX)
    print(result)
    return {"field":field,"result":result}

def job(files,extCache,cell_data,output):
    
#    output.put([i[1] if i else None for i in files])
#    return None
    comm = MPI.COMM_WORLD
#    rank = comm.Get_rank()
#    size = comm.Get_size()
    local = comm.Dup()

    mesh = fcs.Mesh()
    with fcs.XDMFFile(extCache+"cache/meshCache_il2.xdmf") as f:
        f.read(mesh)
    V = fcs.FunctionSpace(mesh,"P",1) 
    subPath = extCache+"cache/boundary_markers_il2.h5"
    boundary_markers = fcs.MeshFunction("size_t",mesh,mesh.topology().dim() - 1)
    with fcs.HDF5File(local,subPath,"r") as f:
        f.read(boundary_markers,"/boundaries")   
    resultList = []
    for no,i in enumerate(files):
        if not i:
            continue
        parameters = {
        "mesh":mesh,
        "V":V,
        "boundary_markers":boundary_markers,
        "cell_data":cell_data
        }
        u = fcs.Function(V)
        with fcs.HDF5File(local,i[1],"r") as f:
            f.read(u,i[2])
        parameters["u"] = copy(u)
        print("reading file {file} ({no}/{tot})".format(file=i[1],no=no,tot=len(files)))
        dataOUT = compute(i[0],i[1],i[2],extCache,**parameters)
        resultList.append(deepcopy(dataOUT))
    output.put(resultList)

def gradient_dump(path,extCache,fields,threads):
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    
    
    scatterList = []
    cell_dump_path = path+"cell_dump.json"
    with open(cell_dump_path,"r") as file:
        d = file.read()
        cell_data = json.loads(d)
    
    for field in fields:
        fileList = os.listdir(path+"sol/distplot")
        bList = [((re.compile(field)).search("{field}(?=.*\.h5)".format(field=i))) for i in fileList]
        fileList = np.delete(fileList,[i for i,e in enumerate(bList) if not e])
        for no,file in enumerate(fileList):
            scatterList.append([no,path+"sol/distplot/"+file,field])
            

    size = ceil(len(scatterList)/threads)
    partitionedList = [scatterList[x:x+size] for x in range(0,len(scatterList),size)]
    resultList = []
    output = mp.Queue(threads)
    
    jobs = [mp.Process(target=job,args=(i,extCache,cell_data,output)) for i in partitionedList]
    for j in jobs:
        j.start()
        
    for j in jobs:
        j.join(600)
    resultList = [output.get(True,10) for j in jobs]
    flattendList = []
    for i in resultList:
        for o in i:
            flattendList.append(o)
    outDict = {}
    for o in flattendList:
        f = o["field"]
        if f in outDict:    
            outDict[f].append(o["result"])
        else:
            outDict[f] = [o["result"]]
    for k,v in outDict.items():
        with open(path+"gradient_dump_"+k+".json","w") as file:
            file.write(json.dumps(v))
                    
#distPlot_dump("/extra/kiwitz/test/data/","/extra/kiwitz/extCache/",["il2","il6"],10)