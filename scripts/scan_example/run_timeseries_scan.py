#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 14:58:01 2019

@author: kiwitz
"""

import random
import time
from copy import deepcopy

import BC as bc
import Entity
import FieldProblem as fp
import MySolver
import SimContainer as SC
import StateManager
import fenics as fcs
import mpi4py.MPI as MPI
import numpy as np
from InternalSolver import InternalSolver
from bcFunctions import cellBC_il2, cellBC_il6, cellBC_infg
from my_debug import message, total_time
from scipy.constants import N_A

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


class RuleBasedSolver(InternalSolver):

    name = "RuleBasedSolver"

    def __init__(self):
        self.transition = (False, 0, "Tn")
        self.il2_threshold = 0.71
        self.il6_threshold = 0.059
        self.infg_threshold = 0.06

    def step(self,t,dt,p,entity=None):

        if entity.type_name == "Tn":
            if self.transition[0]:
                if self.transition[1] <= 0:
                    entity.change_type = self.transition[2]
                else:
                    self.transition = (True,self.transition[1]-dt,self.transition[2])
            elif np.random.rand(1) > 0.5:  # chance for type change; uniform distribution
                draw = np.random.normal(1,0.2)  # time to type change; normal distribution

                if p["surf_c_il2"]*1e9 < self.il2_threshold and p["surf_c_il6"]*1e9 > self.il6_threshold:  # il2neg and il6pos
                    self.transition = (True, draw, "Tfh")
                elif p["surf_c_infg"]*1e9 > self.infg_threshold:  #infg pos
                        self.transition = (True, draw, "Th1")

        return p

def updateState(p, sc, t):
    ran = random.Random()
    ran.seed(t)

    for i, e in enumerate(sc.entity_list):
        draw = ran.random()
        e.set_cell_type(sc.get_entity_type_by_name("Tn"))
        if draw > 1 - p["fraction"]:
            e.set_cell_type(sc.get_entity_type_by_name("Tfh"))
        elif draw > 1 - 2 * p["fraction"]:
            e.set_cell_type(sc.get_entity_type_by_name("Th1"))
        e.p["id"] = e.id

def makeCellListGrid(p, xLine, yLine, zLine):
    cellList = []
    ran = random.Random()
    ran.seed(1)
    for x in xLine:
        for y in yLine:
            for z in zLine:
                p_temp = deepcopy(p)
                cell = Entity.Cell([x, y, z], p_temp["rho"],
               [
                        bc.Integral(cellBC_il2,field_quantity="il2"),
                        bc.Integral(cellBC_il6,field_quantity="il6"),
                        bc.Integral(cellBC_infg, field_quantity="infg")
                ])
                cell.name = "cell_{xCoord}_{yCoord}_{zCoord}".format(xCoord=x, yCoord=y, zCoord=z)
                cell.p = p_temp
                cellList.append(cell)
    return cellList

def setup(p, path, ext_cache=""):

    message("setup")

    """Domain setup"""

    # x = [0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8]
    # y = [-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8]
    # z = [0]

    x = np.arange(-1,1,0.2)
    y = x
    z = [0]

    margin = 0.2

    # domainBC = [
    #     bc.OuterIntegral(lambda u, p: fcs.Constant(0), "!near(x[0],{d_x})".format(d_x=x[0] - margin),
    #                      field_quantity="il2"),
    #     bc.OuterIntegral(outerBC_il2, "near(x[0],{d_x})".format(d_x=x[0] - margin), field_quantity="il2"),
    #     bc.OuterIntegral(outerBC_il6, "near(x[0],{d_x})".format(d_x=x[0] - margin), field_quantity="il6"),
    #     bc.OuterIntegral(lambda u, p: fcs.Constant(0), "!near(x[0],{d_x})".format(d_x=x[0] - margin),
    #                      field_quantity="il6"),
    #     bc.OuterIntegral(lambda u, p: fcs.Constant(0), "near(x[0],{d_x})".format(d_x=x[0] - margin),
    #                      field_quantity="infg"),
    #     bc.OuterIntegral(lambda u, p: fcs.Constant(0), "!near(x[0],{d_x})".format(d_x=x[0] - margin),
    #                      field_quantity="infg")
    # ]
    domainBC = [
        bc.OuterIntegral(lambda u, p: fcs.Constant(0), "true".format(d_x=x[0] - margin),
                         field_quantity="il2"),
        bc.OuterIntegral(lambda u, p: fcs.Constant(0), "true".format(d_x=x[0] - margin),
                         field_quantity="il6"),
        bc.OuterIntegral(lambda u, p: fcs.Constant(0), "true".format(d_x=x[0] - margin),
                         field_quantity="infg"),
    ]


    domain = Entity.DomainCube(
        [
            x[0] - margin,
            y[0] - margin,
            z[0] - margin
        ],
        [
            x[-1] + margin,
            y[-1] + margin,
            z[-1] + margin
        ],
        deepcopy(p), domainBC)


    """IL-2"""
    solver_il2 = MySolver.PoissonSolver()

    fieldProblem_il2 = fp.FieldProblem()
    fieldProblem_il2.field_name = "il2"
    fieldProblem_il2.field_quantity = "il2"

    fieldProblem_il2.set_solver(solver_il2)
    fieldProblem_il2.p = deepcopy(p)

    if not ext_cache == "":
        fieldProblem_il2.ext_cache = ext_cache + "cache/meshCache_il2"

    fieldProblem_il2.set_outer_domain(domain)

    """IL-6"""
    solver_il6 = MySolver.PoissonSolver()

    fieldProblem_il6 = fp.FieldProblem()
    fieldProblem_il6.field_name = "il6"
    fieldProblem_il6.field_quantity = "il6"

    fieldProblem_il6.set_solver(solver_il6)
    fieldProblem_il6.p = deepcopy(p)
    if not ext_cache == "":
        fieldProblem_il6.ext_cache = ext_cache + "cache/meshCache_il2"

    fieldProblem_il6.set_outer_domain(domain)

    """INFg"""
    solver_infg = MySolver.PoissonSolver()

    fieldProblem_infg = fp.FieldProblem()
    fieldProblem_infg.field_name = "infg"
    fieldProblem_infg.field_quantity = "infg"

    fieldProblem_infg.set_solver(solver_infg)
    fieldProblem_infg.p = deepcopy(p)

    if not ext_cache == "":
        fieldProblem_infg.ext_cache = ext_cache + "cache/meshCache_il2"

    fieldProblem_infg.set_outer_domain(domain)

    # Setup
    sc = SC.SimContainer()
    """top level path"""
    sc.path = path

    """adds cell to simulation"""
    for i in makeCellListGrid(p, x, y, z):
        sc.add_entity(i)

    sc.add_field(fieldProblem_il2)
    sc.add_field(fieldProblem_il6)
    sc.add_field(fieldProblem_infg)

    """adds internal solver"""
    sc.add_internal_solver(RuleBasedSolver)

    """sets external path for subdomain markers"""
    if not ext_cache == "":
        sc.initialize(load_subdomain=ext_cache + "cache/boundary_markers_il2.h5")
    else:
        sc.initialize()

    message("initialization complete")
    return sc

def run(sc, p, T, dt, path, **kwargs):

    start_run = time.time()

    if "scan" in kwargs:
        scan = kwargs["scan"]
        stMan = StateManager.StateManager(path)
        stMan.scan_log(scan, p)

        for number, s in enumerate(scan):
            stMan.addCellDump(sc, 0)
            stMan.writeElementTree()

            p = stMan.updateSimContainer(sc, number)
            sc.path = stMan.getScanFolder(number)
            updateState(p, sc, 0)
            sc.init_xdmf_files()
            sc.T = 0
            for n in T:

                # updateState(p,sc,t)
                start = time.process_time()

                sc.step(dt)
                end = time.process_time()
                total_time(end - start, pre = "total time ", post=" for step number ")
                resultPaths = sc.save_fields(n)
                for k, v in resultPaths.items():
                    (distplot, sol, cells) = v
                    stMan.addTimeStep(number, n, sc.T, displot=distplot, sol=sol, field_name=k, cell_list=cells)
                    stMan.writeElementTree()



    total_time(time.time() - start_run, pre = "total time ")

def get_update_dict(dict, update):
    dict_copy = deepcopy(dict)
    dict_copy.update(update)
    return dict_copy

p = {
    "k_on": 1e9 * 111.6 / 60 ** 2,  # 111.6 per hour
    "rho": 0.05,  # mu
    "D": (10 ** 0.5 * 0.01) ** 2,  # mu² per s
    "R_h": 4000 * N_A ** -1 * 1e9,
    "R_l": 100 * N_A ** -1 * 1e9,
    "kd": 0,#0.1/(60*2),
    "q_h": 10 * N_A ** -1 * 1e9,
    "q_l": 1 * N_A ** -1 * 1e9,
    "fraction": 0.05
}

p_bc_defaults = {  # default values for boundary functions
    "R_il2": 0,
    "q_il2": 0,
    "R_il6":0,
    "q_il6":0,
    "R_infg":0,
    "q_infg":0
}

p_boundary = {# outer boundary
        "R_il2_b":p["R_h"],
        "q_il2_b":0,
        "R_il6_b":0,
        "q_il6_b":p["q_h"],
        "R_infg_b":0,
        "q_infg_b":0,
    }

""" cell types"""
Tn = {
        "q_il2": p["q_h"],
        "R_il2": p["R_h"],
        "q_il6": 0,#p["q_l"],
        "R_il6": p["R_l"],
        "R_infg":p["R_l"],
        "q_infg":0,
        "type_int":1
}

Tfh = {
        "q_il2": p["q_l"],
        "R_il2": p["R_h"],
        "q_il6": p["q_h"],
        "R_il6": p["R_h"],
        "R_infg":0,
        "q_infg":0,
        "type_int":2
}

Th1 = {
        "q_il2": 0,#p["q_h"],
        "R_il2": 0,#p["R_l"],
        "q_il6": 0,#p["q_l"],
        "R_il6": 0,#p["R_l"],
        "q_infg": p["q_h"],
        "R_infg": p["R_h"],
        "type_int":3
}

"""scan list"""
scan_default= [
        {
            "name":"scan",

        "parameters":{
            "D":p["D"]
        },
        "entity_types":[
            {
            "parameters":Tn,
             "name":"Tn",
             "internal_solver":"RuleBasedSolver"
             },
            {
            "parameters":Th1,
             "name":"Th1",
             "internal_solver":"RuleBasedSolver"
             },
            {
            "parameters":Tfh,
             "name":"Tfh",
             "internal_solver":"RuleBasedSolver"
             }
        ]

        } for v1 in [1]
]

"""reshapes into 1-d list"""
scan = list(np.ravel(scan_default))



path = "/extra/kiwitz/scan_example_small/"
ext_cache="/extra/kiwitz/scan_example_small_ext_cache/"


p = {**p,**p_bc_defaults,**p_boundary}

T = range(100)
dt = 1

sc = setup(p, path, ext_cache)
run(sc,p,T,dt,path,scan=scan)




