#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 14:58:01 2019

@author: kiwitz
"""

import FieldProblem as fp
import Entity
import MySolver
import fenics as fcs
import numpy as np
import matplotlib.pyplot as plt
import BC as bc
import SimContainer as SC
from bcFunctions import cellBC_il2,cellBC_il6,outerBC_il2,outerBC_il6,absorbing_il2
from copy import copy,deepcopy
import BooleanInternalSolver as intSolver
import json
import os
import time
import random
from scipy.constants import N_A
import xml.etree.ElementTree as ET
import StateManager

import mpi4py.MPI as MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def updateState(p,sc,t):
    ran = random.Random()
    ran.seed(t)
    for i,e in enumerate(sc.entityList):
        p_temp = deepcopy(p)
        draw = ran.random()
        if draw > 1-p["fraction"]:
            p_temp["q_il2"] = p["q_il2_s"]
            p_temp["R_il2"] = p["R_il2_s"]
            
            p_temp["q_il6"] = p["q_il6_s"]
            p_temp["R_il6"] = p["R_il6_s"]
            
            p_temp["type"] = 1
        elif draw > 1-2*p["fraction"]:
            p_temp["q_il2"] = p["q_il2_f"]
            p_temp["R_il2"] = p["R_il2_f"]
            
            p_temp["q_il6"] = p["q_il6_f"]
            p_temp["R_il6"] = p["R_il6_f"]
            
            p_temp["type"] = 2
        else:
            p_temp["q_il2"] = p["q_il2_n"]
            p_temp["R_il2"] = p["R_il2_n"]

            p_temp["q_il6"] = p["q_il6_n"]            
            p_temp["R_il6"] = p["R_il6_n"]
            
            p_temp["type"] = 3
                
        e.p = p_temp
        
def makeCellListGrid(p,xLine,yLine,zLine):
    
    cellList = []
    ran = random.Random()
    ran.seed(1)
    for x in xLine:
        for y in yLine:
            for z in zLine:                    
                p_temp = deepcopy(p)
                cell = Entity.Cell([x,y,z],p_temp["rho"],
               [
                       bc.Integral(cellBC_il2,fieldQuantity="il2"),
                       bc.Integral(cellBC_il6,fieldQuantity="il6")
                ])
                cell.name = "cell_{xCoord}_{yCoord}_{zCoord}".format(xCoord=x,yCoord=y,zCoord=z)
                cell.p = p_temp
                cellList.append(cell)
    return cellList
def scan_log(topLevel,output):
    root = ET.Element("run")
    for o in output:
        scan = ET.SubElement(root,"scan")
        
        path = ET.SubElement(scan,"path")
        path.text = str(o["subfolder"])
        
        parameters = ET.SubElement(scan,"parameters")
        
        constant = ET.SubElement(parameters,"constant")
        for k,v in o["constant"].items():
            par = ET.SubElement(constant,k)
            par.text = str(v)
        dynamic = ET.SubElement(parameters,"dynamic")
        for k,v in o["dynamic"].items():
            par = ET.SubElement(dynamic,k)
            par.text = str(v)
        number = ET.SubElement(scan,"number")
        number.text = str(o["number"])
        
        timeSeries = ET.SubElement(scan,"timeSeries")
    
    tree = ET.ElementTree(element=root)
    tree.write(topLevel+"log.scan")
#    with open(topLevel+"log.scan","w") as f:
#        f.write(json.dumps(output))
        
def run(p,T,domainBC,path,**kwargs):
    
#    p_c = p_global["p_c"]
#    p_d = p_global["p_d"]
#    p_il2 = p_global["p_il2"]
#    p_il6 = p_global["p_il6"]

    p_domain = deepcopy(p)
    p_domain.update({
            "R_il2":p["R_il2_b"],
            "q_il2":p["q_il2_b"],
            "R_il6":p["R_il6_b"],
            "q_il6":p["q_il6_b"]
            })
    domain = Entity.DomainCube([-p["d_x"],-p["d_y"],-p["d_z"]],[p["d_x"],p["d_y"],p["d_z"]],p_domain,domainBC)
    """IL-2"""
    solver_il2 = MySolver.PoissonSolver()
    
    fieldProblem_il2 = fp.FieldProblem()
    fieldProblem_il2.fieldName = "il2"
    fieldProblem_il2.fieldQuantity = "il2"
    
    
    fieldProblem_il2.setSolver(solver_il2)
    fieldProblem_il2.p = deepcopy(p)
    fieldProblem_il2.extCache = kwargs["extCache"]+"cache/meshCache_il2"
#    fieldProblem_il2.meshCached = "cache/meshCache_il2"
    fieldProblem_il2.setOuterDomain(domain)
    
    """IL-6"""
    solver_il6 = MySolver.PoissonSolver()
    
    fieldProblem_il6 = fp.FieldProblem()
    fieldProblem_il6.fieldName = "il6"
    fieldProblem_il6.fieldQuantity = "il6"
    
    
    fieldProblem_il6.setSolver(solver_il6)
    fieldProblem_il6.p = deepcopy(p)
    fieldProblem_il6.extCache = kwargs["extCache"]+"cache/meshCache_il2"
#    fieldProblem_il6.meshCached = "cache/meshCache_il2"
    fieldProblem_il6.setOuterDomain(domain)
#Setup
    sc = SC.SimContainer()
    sc.path = path
    
    for i in makeCellListGrid(p,p["x"],p["y"],p["z"]):
        sc.addEntity(i)

    sc.addField(fieldProblem_il2)
    sc.addField(fieldProblem_il6)
    
    sc.initialize(load_subdomain=kwargs["extCache"]+"cache/boundary_markers_il2.h5")
    print("init complete")
    
    updateState(p,sc,0)
    sc.saveSubdomains()
    print("subdomains saved")
    sc.saveDomain()
    print("domain saved")
    stMan = StateManager.StateManager(path)
    
    if "scan" in kwargs:
        scan = kwargs["scan"]
        stMan.scan_log(scan,p)
        stMan.addCellDump(sc,0)
        stMan.writeElementTree()
        
        for number,s in enumerate(scan):
            p = stMan.updateSimContainer(sc,number)
            sc.path = stMan.getScanFolder(number)            
            updateState(p,sc,0)
            sc.initXdmfFiles()
#            sc.saveSubdomains()
                
            for n,t in enumerate(T):
                updateState(p,sc,t)
                start = time.process_time()
                
                sc.step(1)
                end = time.process_time()
                print("time: "+str(end-start)+"s for step number "+str(n))
                resultPaths = sc.saveFields(n)
                for k,v in resultPaths.items():
                    (distplot,sol) = v
                    stMan.addTimeStep(number,n,t,displot=distplot,sol=sol,fieldName=k)
                    stMan.writeElementTree()
                
    else:
        for n,t in enumerate(T):
            updateState(p,sc,t)
            start = time.process_time()
            sc.step(1)
            end = time.process_time()
            print("time: "+str(end-start)+"s for step number "+str(n))
            sc.saveFields(n)