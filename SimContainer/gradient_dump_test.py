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
def collectFilesFormLog(scanLog):
    pass
def compute(file,extCache,**kwargs):
    
    filePath = file["file"]
    field = file["field"]
    no = file["no"]
    mesh = kwargs["mesh"]    
    u = kwargs["u"]
    V_vec = fcs.VectorFunctionSpace(mesh,"P",1)
    grad = fcs.project(fcs.grad(u),V_vec,solver_type="gmres")
    gradient = fcs.assemble(fcs.sqrt(fcs.dot(grad,grad))*fcs.dX)*10**8
    conc = fcs.assemble(u*fcs.dX)*10**9
    return {"field":field,"gradient":gradient,"concentration":conc,"file":filePath,"dict":file["dict"]}

def job(files,extCache,output):
    

    comm = MPI.COMM_WORLD

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
    for n,file in enumerate(files):
        filePath = file["file"]
        field = file["field"]
        no = file["no"]
        if not file:
            continue
        parameters = {
        "mesh":mesh,
        "V":V,
        "boundary_markers":boundary_markers
        }
        u = fcs.Function(V)
        with fcs.HDF5File(local,filePath,"r") as f:
            f.read(u,"/"+field)
        parameters["u"] = copy(u)
        print("reading file {file} ({n}/{tot})".format(file=filePath,n=n,tot=len(files)))
        dataOUT = compute(file,extCache,**parameters)
        resultList.append(deepcopy(dataOUT))
    output.put(resultList)

def gradient_dump(path,extCache,fields,threads,**kwargs):
    
    for i in os.listdir(path):
        if i.endswith(".scan"):
            with open(path+"/"+i,"r") as file:
                d = file.read()
                scanLog = json.loads(d)
            break
    timesteps = []
    scatterList = []
    
    for singleScan in scanLog: ## loads timestep logs
        subFolder = singleScan["subfolder"]
        for subFile in os.listdir(subFolder):
            if subFile.endswith(".log"):
                with open(subFolder+subFile,"r") as f:
                    d = f.read()
                    d = "["+d[:-1]+"]"
                    timestep = json.loads(d)
                    timesteps.append(timestep)
            if subFile == "sol":
                for field in fields:
                    fileList = os.listdir(subFolder+subFile+"/distplot")
                    bList = [((re.compile(field)).search("{field}(?=.*\.h5)".format(field=i))) for i in fileList]
                    fileList = np.delete(fileList,[i for i,e in enumerate(bList) if not e])
                    for no,file in enumerate(fileList):
                        scatterList.append({"no":no,
                                            "file":subFolder+subFile+"/distplot/"+file,
                                            "field":field,
                                            "dict":singleScan["dict"]})
    print("scatter")    

    size = ceil(len(scatterList)/threads)
    partitionedList = [scatterList[x:x+size] for x in range(0,len(scatterList),size)]
    resultList = []
    output = mp.Queue(threads)
    
    jobs = [mp.Process(target=job,args=(i,extCache,output)) for i in partitionedList]
    for j in jobs:
        j.start()
        
    for j in jobs:
        j.join(600)
    resultList = [output.get(True,10) for j in jobs]
    flattendList = []
    for i in resultList:
        for o in i:
            flattendList.append(o)

    with open(path+"gradient_dump.json","w") as file:
            file.write(json.dumps(flattendList))