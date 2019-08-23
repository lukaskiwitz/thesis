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

def compute(no,h5File,field,path,extCache,**kwargs):
    mesh = kwargs["mesh"]
    boundary_markers = kwargs["boundary_markers"]  
    cell_data =  kwargs["cell_data"]
    tmpIndexList = []
    result = []
    u = kwargs["u"]
    
    ds = fcs.Measure("ds", domain=mesh, subdomain_data=boundary_markers)
    for i in cell_data:   
        
        r = i["center"][0]
        v = (fcs.assemble(u*ds(i["patch"]))/(4*np.pi*0.05**2)*10**9)
        
        if r in tmpIndexList:
            result[tmpIndexList.index(r)].append({"x":r,"v":v})
        else:
            result.append([{"x":r,"v":v}])
            tmpIndexList.append(r)
    return {"field":field,"result":result}
#    print(h5File)
#    print(u.vector()[5])
#    print(v)
#    print(result[10][0])
##    return result[10][0]





def distPlot_dump(path,extCache,fields):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = 1#comm.Get_size()
    local = comm.Dup()
    scatterList = []
    partitionedList = []
    scatterList = []
    if rank == 0:    
        mesh = fcs.Mesh()
        with fcs.XDMFFile(extCache+"cache/meshCache_il2.xdmf") as f:
            f.read(mesh)
        V = fcs.FunctionSpace(mesh,"P",1) 
        subPath = extCache+"cache/boundary_markers_il2.h5"
        boundary_markers = fcs.MeshFunction("size_t",mesh,mesh.topology().dim() - 1)
        with fcs.HDF5File(local,subPath,"r") as f:
            f.read(boundary_markers,"/boundaries")
        cell_dump_path = path+"cell_dump.json"
        with open(cell_dump_path,"r") as file:
            d = file.read()
            cell_data = json.loads(d)
        
        for field in fields:
            fileList = os.listdir(path+"sol/distplot")
            bList = [((re.compile(field)).search("{field}(?=.*\.h5)".format(field=i))) for i in fileList]
            fileList = np.delete(fileList,[i for i,e in enumerate(bList) if not e])
                      
            for no,file in enumerate(fileList):
                scatterList.append([no,file,field])
        partitionedList = [scatterList[x:x+size] for x in range(0,len(scatterList),size)]
        partitionedList = [e + [None] * (size - len(e)) for e in partitionedList]
    
        resultList = []
        for no,i in enumerate(partitionedList):
            if rank == 0:
                dataIN = i
            else:
                dataIN = None
            dataIN  = comm.scatter(dataIN,root = 0)
#            dataOUT = dataIN
            if dataIN:
                parameters = {
                "mesh":mesh,
                "V":V,
                "boundary_markers":boundary_markers,
                "cell_data":cell_data
                }
                u = fcs.Function(V)
                with fcs.HDF5File(local,path+"sol/distplot/"+dataIN[1],"r") as f:
                    f.read(u,dataIN[2])
                parameters["u"] = copy(u)
                print("reading file {file} ({no}/{tot})".format(file=dataIN[1],no=no,tot=len(scatterList)))
                dataOUT = compute(dataIN[0],dataIN[1],dataIN[2],path,extCache,**parameters)
    #            dataOUT = 
    #
            else:
                parameters = {
                "mesh":mesh,
                "V":V,
                "boundary_markers":boundary_markers,
                "cell_data":cell_data
                }
    #            u = fcs.Function(V)
    #            with fcs.HDF5File(local,path+"sol/distplot/"+partitionedList[0][0][1],"r") as f:
    #                f.read(u,partitionedList[0][0][2])
    #            parameters["u"] = copy(u)
    #            compute(partitionedList[0][0][0],partitionedList[0][0][1],partitionedList[0][0][2],path,extCache,**parameters)#dummy computation
    #            dataOUT = None
                            
#            data = dataOUT
            data = comm.gather(dataOUT,root = 0)
            if rank == 0:
                pass#print(data)
                resultList.append(deepcopy(data))
    
        if rank == 0:
            flattendList = np.array(resultList).flatten()
            flattendList = np.delete(flattendList,[i for i,e in enumerate(flattendList) if not e])
            flattendList = list(flattendList)
            outDict = {}
    #    if rank == 0 and resultList[0][0] == resultList[0][1]:
    #        print("TRUE!!!!!!!!!!!")
            
            for o in flattendList:
                f = o["field"]
                if f in outDict:    
                    outDict[f].append(o["result"])
                else:
                    outDict[f] = [o["result"]]
            for k,v in outDict.items():
                with open(path+"distPlot_dump_"+k+".json","w") as file:
                    file.write(json.dumps(v))