#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 14:58:01 2019

@author: kiwitz
"""

import FieldProblem as fp
import Entity
import MySolver
import numpy as np
import BC as bc
import SimContainer as SC
from bcFunctions import cellBC_il2,cellBC_il6, cellBC_infg
from copy import deepcopy
import time
import random
import StateManager
from InternalSolver import InternalSolver
import mpi4py.MPI as MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

class RuleBasedSolver(InternalSolver):

    name = "RuleBasedSolver"

    def __init__(self):
        self.transition = (False, 0, "Tn")
    def step(self,t,dt,p,entity=None):

        if entity.type_name == "Tn":
            if self.transition[0]:
                if self.transition[1] <= 0:
                    entity.change_type = self.transition[2]
                else:
                    self.transition = (True,self.transition[1]-dt,self.transition[2])
            elif np.random.rand(1) > 0.5:  # chance for type change; uniform distribution
                draw = np.random.normal(1,0.2)  # time to type change; normal distribution

                if p["surf_c_il2"]*10e9 < 0.7 and p["surf_c_il6"]*10e9 > 0.35:  # il2neg and il6pos
                    self.transition = (True, draw, "Tfh")
                elif p["surf_c_infg"]*10e9 > 0.15:  #infg pos
                        self.transition = (True, draw, "Th1")

        return p



def updateState(p,sc,t):

    ran = random.Random()
    ran.seed(t)

    for i,e in enumerate(sc.entity_list):

        draw = ran.random()
        e.set_cell_type(sc.get_entity_type_by_name("Tn"))
        if draw > 1-p["fraction"]:
            e.set_cell_type(sc.get_entity_type_by_name("Tfh"))
        elif draw > 1-2*p["fraction"]:
            e.set_cell_type(sc.get_entity_type_by_name("Th1"))
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
                bc.Integral(cellBC_il2,field_quantity="il2"),
                bc.Integral(cellBC_il6,field_quantity="il6"),
                bc.Integral(cellBC_infg, field_quantity="infg")
                ])
                cell.name = "cell_{xCoord}_{yCoord}_{zCoord}".format(xCoord=x,yCoord=y,zCoord=z)
                cell.p = p_temp
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
    fieldProblem_il2.field_name = "il2"
    fieldProblem_il2.field_quantity = "il2"
    
    
    fieldProblem_il2.set_solver(solver_il2)
    fieldProblem_il2.p = deepcopy(p)

    if "extCache" in kwargs:
        fieldProblem_il2.ext_cache = kwargs["extCache"]+"cache/meshCache_il2"

    fieldProblem_il2.set_outer_domain(domain)
    
    """IL-6"""
    solver_il6 = MySolver.PoissonSolver()
    
    fieldProblem_il6 = fp.FieldProblem()
    fieldProblem_il6.field_name = "il6"
    fieldProblem_il6.field_quantity = "il6"
    
    
    fieldProblem_il6.set_solver(solver_il6)
    fieldProblem_il6.p = deepcopy(p)
    if "extCache" in kwargs:
        fieldProblem_il6.ext_cache = kwargs["extCache"]+"cache/meshCache_il2"

    fieldProblem_il6.set_outer_domain(domain)

    """INFg"""
    solver_infg = MySolver.PoissonSolver()

    fieldProblem_infg = fp.FieldProblem()
    fieldProblem_infg.field_name = "infg"
    fieldProblem_infg.field_quantity = "infg"

    fieldProblem_infg.set_solver(solver_infg)
    fieldProblem_infg.p = deepcopy(p)

    if "extCache" in kwargs:
        fieldProblem_infg.ext_cache = kwargs["extCache"] + "cache/meshCache_il2"

    fieldProblem_infg.set_outer_domain(domain)
#Setup
    sc = SC.SimContainer()
    sc.path = path
    
    for i in makeCellListGrid(p,p["x"],p["y"],p["z"]):
        sc.add_entity(i)

    sc.add_field(fieldProblem_il2)
    sc.add_field(fieldProblem_il6)
    sc.add_field(fieldProblem_infg)

    sc.add_internal_solver(RuleBasedSolver)


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
            sc.T = 0
            for n in T:

                # updateState(p,sc,t)
                start = time.process_time()
                
                sc.step(dt)
                end = time.process_time()
                print("time: "+str(end-start)+"s for step number "+str(n))
                resultPaths = sc.save_fields(n)
                for k,v in resultPaths.items():
                    (distplot, sol,cells) = v
                    print(len(cells))
                    stMan.addTimeStep(number, n, sc.T, displot=distplot, sol=sol, field_name=k, cell_list=cells)
                    stMan.writeElementTree()

    # sc.save_subdomains()
    # print("subdomains saved")
    # sc.save_domain()
    # print("domain saved")