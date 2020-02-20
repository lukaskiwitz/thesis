#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 14:58:01 2019

@author: kiwitz
"""

import getpass
import random
import time
from copy import deepcopy

import fenics as fcs
import matplotlib.pyplot as plt
import mpi4py.MPI as MPI
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.constants import N_A

import BC as bc
import Entity
import FieldProblem as fp
import MySolver
import SimContainer as SC
import StateManager
from PostProcess import PostProcessor
from bcFunctions import cellBC_il2
from my_debug import message

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()



def updateState(p, sc, t):
    ran = random.Random()
    ran.seed(t)

    for i, e in enumerate(sc.entity_list):

        draw = ran.random()
        e.set_cell_type(sc.get_entity_type_by_name("Tn"))
        if draw > 1 - p["fraction"]:
            e.set_cell_type(sc.get_entity_type_by_name("sec"))
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
                                       bc.Integral(cellBC_il2, field_quantity="il2")
                                   ])
                cell.name = "cell_{xCoord}_{yCoord}_{zCoord}".format(xCoord=x, yCoord=y, zCoord=z)
                cell.p = p_temp
                cellList.append(cell)
    return cellList

def setup(p, path, ext_cache=""):

    message("setup")

    """Domain setup"""
    x = np.round(np.arange(-0.8,0.8,0.2),2)
    y = x
    z = [0]

    margin = 0.2

    domainBC = [
        bc.OuterIntegral(lambda u, p: fcs.Constant(0), "true", field_quantity="il2")
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

    # Setup
    sc = SC.SimContainer()
    """top level path"""
    sc.path = path

    """adds cell to simulation"""
    for i in makeCellListGrid(p, x, y, z):
        sc.add_entity(i)

    sc.add_field(fieldProblem_il2)

    """sets external path for subdomain markers"""
    if not ext_cache == "":
        sc.initialize(load_subdomain=ext_cache + "cache/boundary_markers_il2.h5")
    else:
        sc.initialize()

    message("initialization complete")
    return sc

def run(sc, p, T, dt, path, **kwargs):

    start = time.process_time()

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
                message("time: " + str(end - start) + "s for step number " + str(n))
                resultPaths = sc.save_fields(n)
                for k, v in resultPaths.items():
                    (distplot, sol, cells) = v
                    stMan.addTimeStep(number, n, sc.T, displot=distplot, sol=sol, field_name=k, cell_list=cells)
                    stMan.writeElementTree()

    end = time.process_time()
    message("--------------------- total Time: {t} m ---------------------".format(t=(end - start) / 60))

def get_update_dict(dict, update):
    dict_copy = deepcopy(dict)
    dict_copy.update(update)
    return dict_copy

p_bc_defaults = {  # default values for boundary functions
    "R_il2": 0,
    "q_il2": 0
}

p = {
    "k_on": 1e9 * 111.6 / 60 ** 2,  # 111.6 per hour
    "rho": 0.05,  # mu
    "D": (10 ** 0.5 * 0.01) ** 2,  # mu² per s
    "R_h": 400 * N_A ** -1 * 1e9,
    "R_l": 10 * N_A ** -1 * 1e9,
    "kd": 0,#0.1/(60*2),
    "q_h": 10 * N_A ** -1 * 1e9,
    "q_l": 1 * N_A ** -1 * 1e9,
    "fraction": 0.1
}

p_boundary = {  # outer boundary
    "R_il2_b": p["R_l"],
    "q_il2_b": p["q_l"]
}

"""cell types"""
Tn = {
        "q_il2": p["q_l"],
        "R_il2": p["R_l"],
        "type_int":1
}

sec = {
        "q_il2": p["q_h"],
        "R_il2": p["R_h"],
        "type_int":2,
}

"""scan list"""
scan_default= [
    [
        {
            "name":"scan",

        "parameters":{
            "fraction":v2,
            "D":p["D"]*v1
        },
        "entity_types":[
            {
            "parameters":Tn,
             "name":"Tn",
             "internal_solver":""
             },
            {
            "parameters":sec,
             "name":"sec",
             "internal_solver":""
             }
        ]

        } for v1 in np.logspace(-1,1,3)
    ] for v2 in np.linspace(0,0.99,10)
]

"""reshapes into 1-d list"""
scan = list(np.ravel(scan_default))


user = getpass.getuser()
path = "/extra/{u}/scan_example/".format(u=user)
ext_cache="/extra/{u}/scan_example_ext_cache/".format(u=user)


# p = {**p,**p_bc_defaults,**p_boundary}
#
# T = range(1)
# dt = 1
#
# sc = setup(p, path, ext_cache)
# run(sc,p,T,dt,path,scan=scan)
#

pp = PostProcessor(path)
pp.write_post_process_xml(4)

pp.make_dataframes(kde=False)
pp.global_dataframe.to_hdf(path + 'global_df.h5', key="data", mode="w")
pp.cell_dataframe.to_hdf(path+"cell_df.h5", key="df", mode="w")
pp.cell_stats.to_hdf(path+"cell_stats_df.h5", key="df", mode="w")

global_df = pd.read_hdf(path+"global_df.h5", mode="r")
cell_df= pd.read_hdf(path+"cell_df.h5", mode="r")
cell_stats = pd.read_hdf(path+"cell_stats_df.h5", mode="r")

plt.figure()
sns.lineplot(x="fraction", y="concentration", data=global_df,hue="D")
plt.show()
plt.figure()
sns.lineplot(x="fraction", y="gradient", data=global_df,hue="D")
plt.show()
plt.figure()
sns.lineplot(x="fraction", y="sd", data=global_df,hue="D")
plt.show()

