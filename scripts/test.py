# -*- coding: utf-8 -*-


import time
from copy import deepcopy

import fenics as fcs
import mpi4py.MPI as MPI
import numpy as np

import BC as bc
from bcFunctions import outerBC_il2, outerBC_il6
from run_scan import run
from sim_parameters import Tn, Th1, Tfh, p, p_d, p_sim
from my_debug import message

x = [-0.4, -0.2, 0, 0.2, 0.4]
y = x
z = y

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


def get_update_dict(dict, update):
    dict_copy = deepcopy(dict)
    dict_copy.update(update)
    return dict_copy

scan_default= [
        {
            "name":"Diffusion_test",

        "parameters":{
            "D":p["D"]*v
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

        } for v in np.logspace(-1,1,4)
    ]


scan_list = [
    scan_default
]



for scan in scan_list:
    def bD(n):
        return n/10

    p_boundary = {
        "R_il2_b":p["R_h"],
        "q_il2_b":0,
        "R_il6_b":0,
        "q_il6_b":p["q_h"],
        "R_infg_b":0,
        "q_infg_b":0,
    }

    domainBC = [
            bc.OuterIntegral(lambda u,p: fcs.Constant(0),"!near(x[0],{d_x})".format(d_x=p_d["x"][0]-p_d["margin"]),field_quantity="il2"),
            bc.OuterIntegral(outerBC_il2,"near(x[0],{d_x})".format(d_x=p_d["x"][0]-p_d["margin"]),field_quantity="il2"),
            bc.OuterIntegral(outerBC_il6,"near(x[0],{d_x})".format(d_x=p_d["x"][0]-p_d["margin"]),field_quantity="il6"),
            bc.OuterIntegral(lambda u,p: fcs.Constant(0),"!near(x[0],{d_x})".format(d_x=p_d["x"][0]-p_d["margin"]),field_quantity="il6"),
            bc.OuterIntegral(lambda u,p: fcs.Constant(0),"near(x[0],{d_x})".format(d_x=p_d["x"][0]-p_d["margin"]),field_quantity="infg"),
            bc.OuterIntegral(lambda u,p: fcs.Constant(0),"!near(x[0],{d_x})".format(d_x=p_d["x"][0]-p_d["margin"]),field_quantity="infg")
            ]


    T = range(1)
    start = time.process_time()
    run({**p,**p_d,**p_sim,**p_boundary},T,1,domainBC,"/extra/kiwitz/parameter_scan_{n}/".format(n=scan[0]["name"]),extCache="/extra/kiwitz/Diffusion_test_ext",scan=scan)
    end = time.process_time()
    # if rank == 0:
    message("--------------------- total Time: {t} m ---------------------".format(t=(end-start)/60))







