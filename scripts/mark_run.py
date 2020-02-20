#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 14:58:01 2019

@author: kiwitz
"""

import random
from copy import deepcopy

import BC as bc
import Entity
import FieldProblem as fp
import MySolver
import SimContainer as SC
from bcFunctions import cellBC_il2, cellBC_il6
from my_debug import message


def updateState(p,sc,t):
    ran = random.Random()
    ran.seed(t)
    for i,e in enumerate(sc.entityList):
        p_temp = deepcopy(p)
        draw = ran.random()
        if draw > 3/4:
            p_temp["q_il2"] = p["q_il2_s"]#two molecules per second
            p_temp["R_il2"] = p["R_il2_s"]
            p_temp["type"] = 1
        elif draw > 2/4:
            p_temp["R_il2"] = p["R_il2_f"]
            p_temp["type"] = 2
        else:
            p_temp["R_il2"] = p["R_il2_n"]
            p_temp["type"] = 3
            
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
                       bc.Integral(cellBC_il2,field_quantity="il2"),
                       bc.Integral(cellBC_il6,field_quantity="il6")
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
    fieldProblem_il2.field_name = "il2"
    fieldProblem_il2.field_quantity = "il2"
    
    
    fieldProblem_il2.setSolver(solver_il2)
    fieldProblem_il2.p = deepcopy(p)
    fieldProblem_il2.extCache = kwargs["extCache"]+"cache/meshCache_il2"
#    fieldProblem_il2.meshCached = "cache/meshCache_il2"
    fieldProblem_il2.setOuterDomain(domain)
    
#Setup
    sc = SC.SimContainer()
    sc.path = path
    
    for i in makeCellListGrid(p,p["x"],p["y"],p["z"]):
        sc.addEntity(i)

    sc.addField(fieldProblem_il2)
    
    sc.initialize(load_subdomain=kwargs["extCache"]+"cache/boundary_markers_il2.h5")
    message("init complete")
    
    for t in T:
        updateState(p,sc,t)
        sc.saveSubdomains()
        message("subdomains saved for timestep {t}".format(t=t))
#        sc.saveDomain()
#        message("domain saved for timestep {t}".format(t=t))