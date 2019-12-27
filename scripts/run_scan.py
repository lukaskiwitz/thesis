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
from bcFunctions import cellBC_il2,cellBC_il6,outerBC_il2,outerBC_il6,absorbing_il2, cellBC_infg
from copy import copy,deepcopy
import BooleanInternalSolver as intSolver
import json
import os
import time
import random
from scipy.constants import N_A
import xml.etree.ElementTree as ET
import StateManager
from CellType import CellType
from InternalSolver import InternalSolver,ODE_Solver, RuleBasedSolver
import cell_types

import mpi4py.MPI as MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def make_cell_type(p,update,solver,name):
    p_temp = deepcopy(p)
    p_temp.update(update)
    p_temp["type_name"]=name
    return CellType(p_temp,solver,name)


def updateState(p,sc,t):



    ran = random.Random()
    ran.seed(t)
    for i,e in enumerate(sc.entity_list):

        draw = ran.random()
        e.set_cell_type(cell_types.Tn)
        if draw > 1-p["fraction"]:
            #e.p["q_il2"] = p["q_h"]
            #e.p["R_il2"] = p["R_h"]
            e.set_cell_type(cell_types.Tfh)
        elif draw > 1-2*p["fraction"]:
            #e.p["q_infg"] = p["q_h"]
            #e.p["R_infg"] = p["R_l"]
            e.set_cell_type(cell_types.Th1)

        e.p["id"] = e.id
        
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
                bc.Integral(cellBC_il6,fieldQuantity="il6"),
                bc.Integral(cellBC_infg, fieldQuantity="infg")
                ])
                cell.name = "cell_{xCoord}_{yCoord}_{zCoord}".format(xCoord=x,yCoord=y,zCoord=z)
                # cell.p = p_temp
                cellList.append(cell)
    return cellList

        
def run(p, T, dt, domainBC, path, **kwargs):

    p_domain = deepcopy(p)
    p_domain.update({
            "R_il2":p["R_il2_b"],
            "q_il2":p["q_il2_b"],
            "R_il6":p["R_il6_b"],
            "q_il6":p["q_il6_b"],
            "R_infg": p["R_infg_b"],
            "q_infg": p["q_infg_b"]
            })

    domain = Entity.DomainCube(
        [
        p["x"][0]-p["margin"],
        p["y"][0]-p["margin"],
        p["z"][0]-p["margin"]
        ],
        [
        p["x"][-1]+p["margin"],
        p["y"][-1]+p["margin"],
        p["z"][-1]+p["margin"]
        ],
        p_domain,domainBC)
    """IL-2"""
    solver_il2 = MySolver.PoissonSolver()
    
    fieldProblem_il2 = fp.FieldProblem()
    fieldProblem_il2.fieldName = "il2"
    fieldProblem_il2.fieldQuantity = "il2"
    
    
    fieldProblem_il2.setSolver(solver_il2)
    fieldProblem_il2.p = deepcopy(p)

    if "extCache" in kwargs:
        fieldProblem_il2.extCache = kwargs["extCache"]+"cache/meshCache_il2"

    fieldProblem_il2.setOuterDomain(domain)
    
    """IL-6"""
    solver_il6 = MySolver.PoissonSolver()
    
    fieldProblem_il6 = fp.FieldProblem()
    fieldProblem_il6.fieldName = "il6"
    fieldProblem_il6.fieldQuantity = "il6"
    
    
    fieldProblem_il6.set_solver(solver_il6)
    fieldProblem_il6.p = deepcopy(p)
    if "extCache" in kwargs:
        fieldProblem_il6.extCache = kwargs["extCache"]+"cache/meshCache_il2"

    fieldProblem_il6.set_outer_domain(domain)

    """INFg"""
    solver_infg = MySolver.PoissonSolver()

    fieldProblem_infg = fp.FieldProblem()
    fieldProblem_infg.fieldName = "infg"
    fieldProblem_infg.fieldQuantity = "infg"

    fieldProblem_infg.set_solver(solver_infg)
    fieldProblem_infg.p = deepcopy(p)

    if "extCache" in kwargs:
        fieldProblem_infg.extCache = kwargs["extCache"] + "cache/meshCache_il2"

    fieldProblem_infg.set_outer_domain(domain)
#Setup
    sc = SC.SimContainer()
    sc.path = path
    
    for i in makeCellListGrid(p,p["x"],p["y"],p["z"]):
        sc.add_entity(i)

    sc.add_field(fieldProblem_il2)
    sc.add_field(fieldProblem_il6)
    sc.add_field(fieldProblem_infg)

    updateState(p, sc, 0)
    if "extCache" in kwargs:
        sc.initialize(load_subdomain=kwargs["extCache"]+"cache/boundary_markers_il2.h5")
    else:
        sc.initialize()

    print("init complete")
    


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
            sc.init_xdmf_files()
#            sc.saveSubdomains()
                
            for n in T:
                # updateState(p,sc,t)
                start = time.process_time()
                
                sc.step(dt)
                end = time.process_time()
                print("time: "+str(end-start)+"s for step number "+str(n))
                resultPaths = sc.save_fields(n)
                for k,v in resultPaths.items():
                    (distplot, sol,cells) = v
                    stMan.addTimeStep(0, n, sc.T, displot=distplot, sol=sol, fieldName=k, cell_list=cells)
                    stMan.writeElementTree()

    sc.save_subdomains()
    print("subdomains saved")
    sc.save_domain()
    print("domain saved")