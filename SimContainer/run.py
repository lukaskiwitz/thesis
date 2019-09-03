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

def updateState(p,sc,t):
    ran = random.Random()
    ran.seed(t)
    for i,e in enumerate(sc.entityList):
        p_temp = deepcopy(p)
        draw = ran.random()
        if draw > 3/4:
            p_temp["q_il2"] = p["q_il2_s"]#two molecules per second
            p_temp["R_il2"] = p["R_il2_s"]
        elif draw > 2/4:
            p_temp["R_il2"] = p["R_il2_f"]
        else:
            p_temp["R_il2"] = p["R_il2_n"]
            
        draw = ran.random()
        if draw > 3/4:
            p_temp["q_il6"] = p["q_il6_s"]#two molecules per second
            p_temp["R_il6"] = p["R_il6_s"]
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
    
    for n,t in enumerate(T):
        updateState(p,sc,t)
        start = time.process_time()
        sc.step(1)
        end = time.process_time()
        print("time: "+str(end-start)+"s for step number "+str(n))
        sc.saveFields(n)