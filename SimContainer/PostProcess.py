#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 12:56:59 2019

@author: kiwitz
"""
import StateManager as ST

import fenics as fcs
import numpy as np
import matplotlib.pyplot as plt
import h5py
import json
import os
import re
import lxml.etree as ET
from copy import copy,deepcopy
import mpi4py.MPI as MPI
import multiprocessing as mp
from math import ceil

from myDictSorting import groupByKey

class PostProcessor:
    def __init__(self,path):
        self.stateManager = ST.StateManager(path)
        self.stateManager.loadXML()
        
        self.pDicts = []
        self.cellDump = []
        self.extCache = self.stateManager.elementTree.find("/cellDump/field/extCache").text
        self.subdomaincache = self.stateManager.elementTree.find("/cellDump/field/subdomains").text
        for s in self.stateManager.elementTree.findall("/scans/scan"):
            self.pDicts.append(self.stateManager.getParametersFromElement(s))
            
        for s in self.stateManager.elementTree.findall("/cellDump/field/cell"):
            patch = int(s.find("patch").text)
            center = json.loads(s.find("center").text)
            self.cellDump.append({"patch":patch,"center":center})
    def compute(self,file,extCache,**kwargs):
        
        
    
        filePath = file["file"]
        field = file["field"]
        no = file["no"]
        mesh = kwargs["mesh"]    
        u = kwargs["u"]
        cellData = kwargs["cellData"]
        boundary_markers = kwargs["boundary_markers"]
        if "toCompute" in kwargs:
            toCompute = kwargs["toCompute"]
        
        else :
            toCompute = []
        
        result = ET.Element("file")
        globalResults = ET.SubElement(result,"global")
        cellResults = ET.SubElement(result,"cellResults")
        
        result.set("field",str(field))
        result.set("path",str(filePath))
        result.set("dynamic",json.dumps(file["dynamic"]))
        result.set("scanIndex",str(kwargs["scanIndex"]))
        result.set("timeIndex",str(kwargs["timeIndex"]))
        
        
    
        if "gradient" in toCompute:
            V_vec = fcs.VectorFunctionSpace(mesh,"P",1)
            grad = fcs.project(fcs.grad(u),V_vec,solver_type="gmres")
            gradient = fcs.assemble(fcs.sqrt(fcs.dot(grad,grad))*fcs.dX)*10**8
            gradientResult = ET.SubElement(globalResults,"gradient")
            gradientResult.text = str(gradient)
        if "concentration" in toCompute:
            conc = fcs.assemble(u*fcs.dX)*10**9
            concResult = ET.SubElement(globalResults,"concentration")
            concResult.text = str(conc)
    
        if "surfaceConcentration" in toCompute:
            for i,cell in enumerate(cellData):
    
                cellElement = ET.SubElement(cellResults,"cell")
                patch = ET.SubElement(cellElement,"patch")
                center = ET.SubElement(cellElement,"center")
                patch.text = str(cell["patch"])
                center.text = json.dumps(list(cell["center"]))
                
                
                ds = fcs.Measure("ds", domain=mesh, subdomain_data=boundary_markers)
                v = (fcs.assemble(u*ds(cell["patch"]))/(4*np.pi*0.05**2)*10**9)
                
                surfaceConcentration= ET.SubElement(cellElement,"surfaceConcentration")
                surfaceConcentration.text = str(v)
                    
        return ET.tostring(result)
             
    def job(self,files,extCache,subdomaincache,cellData,output,toCompute,index):
        try:
            comm = MPI.COMM_WORLD
            
            local = comm.Dup()
            
            mesh = fcs.Mesh()
            with fcs.XDMFFile(extCache+".xdmf") as f:
                f.read(mesh)
            V = fcs.FunctionSpace(mesh,"P",1) 
            subPath = subdomaincache
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
                "boundary_markers":boundary_markers,
                "toCompute":toCompute,
                "cellData":cellData,
                "timeIndex":file["timeIndex"],
                "scanIndex":file["scanIndex"]
                }
                u = fcs.Function(V)
                with fcs.HDF5File(local,filePath,"r") as f:
                    f.read(u,"/"+field)
                parameters["u"] = copy(u)
                print("reading file {file} ({n}/{tot})".format(file=filePath,n=n,tot=len(files)))
                dataOUT = self.compute(file,extCache,**parameters)
                resultList.append(deepcopy(dataOUT))
            print("thread no {index}".format(index=index))
        
            output.put(resultList)   
        except Exception as e:
            print("error")
            
    
    def dump(self,path,threads,**kwargs):
#        timesteps = []
        scatterList = []
        
#        for field in fields:
        cellData = self.stateManager.elementTree.findall("/cellDump/field[@name='il2']/cell")
        cellData = [{"patch":int(i.find("patch").text),"center":json.loads(i.find("center").text)} for i in cellData]
        for scan in self.stateManager.elementTree.findall("scans/scan"): ## loads timestep logs
            dynamic = scan.findall("parameters/dynamic/parameter")
            dynamic= [{"name":i.get("name"),"value":i.text} for i in dynamic]
            
            for step in scan.findall("timeSeries/field/step"):
                field = step.getparent().get("name")
                scatterList.append({"no":0,
                                "file":step.find("distPlotPath").text,
                                "field":field,
                                "dynamic":dynamic,
                                "scanIndex":scan.get("i"),
                                "timeIndex":step.get("n")
                                })
        print("scatter")    
#        scatterList = scatterList[0:16]
        size = ceil(len(scatterList)/threads)
        partitionedList = [scatterList[x:x+size] for x in range(0,len(scatterList),size)]
        resultList = []
        output = mp.Queue(threads)
        toCompute = ["concentration","gradient","surfaceConcentration"]
        jobs = [mp.Process(target=self.job,args=(i,self.extCache,self.subdomaincache,deepcopy(cellData),output,toCompute,index)) for index,i in enumerate(partitionedList)]
        for j in jobs:
            j.start()
        
        for j in jobs:
            try:
                j.join(60)
            except(Exception):
                print("Join Timeout")
        print("joined jobs")
        resultList = [output.get(True,60) for j in jobs]
        print("collected output from threads")
        flattendList = []
        
        for i in resultList:
            for o in i:
                flattendList.append(ET.XML(o))
        indexedList = [
                {"scanIndex":i.get("scanIndex"),
                 "timeIndex":i.get("timeIndex"),
                 "entry":i} 
                for i in flattendList]
        indexedList = groupByKey(indexedList,["scanIndex"])
        for i,e in enumerate(indexedList):
            indexedList[i] = groupByKey(e,["timeIndex"])
        
        
        postProcessResult = ET.Element("postProcess")
        tree = ET.ElementTree(element=postProcessResult)
        for s in indexedList:
            scan =ET.SubElement(postProcessResult,"scan")
#            print(s[0])
            scan.set("i",str(s[0][0]["scanIndex"]))
            for t in s:
                for i in t:
                    time = ET.SubElement(scan,"timeStep")
                    time.set("i",i["timeIndex"])
                    time.append(i["entry"])
            
     
        tree.write(path+"postProcess.xml",pretty_print=True)
                
path = "/extra/kiwitz/results_parameter_scan_qStrength /"
pp = PostProcessor(path)
pp.dump(path,64)
