# -*- coding: utf-8 -*-

import BC as bc
import fenics as fcs
from bcFunctions import outerBC_il2, outerBC_il6
import time
from run_scan import run
from copy import deepcopy

import mpi4py.MPI as MPI
from sim_parameters import Tn,Th1,Tfh,p, p_d, p_sim

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


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

def get_update_dict(dict, update):
    dict_copy = deepcopy(dict)
    dict_copy.update(update)
    return dict_copy
domainBC = [
        bc.OuterIntegral(lambda u,p: fcs.Constant(0),"!near(x[0],{d_x})".format(d_x=p_d["x"][0]-p_d["margin"]),field_quantity="il2"),
        bc.OuterIntegral(outerBC_il2,"near(x[0],{d_x})".format(d_x=p_d["x"][0]-p_d["margin"]),field_quantity="il2"),
        bc.OuterIntegral(outerBC_il6,"near(x[0],{d_x})".format(d_x=p_d["x"][0]-p_d["margin"]),field_quantity="il6"),
        bc.OuterIntegral(lambda u,p: fcs.Constant(0),"!near(x[0],{d_x})".format(d_x=p_d["x"][0]-p_d["margin"]),field_quantity="il6"),
        bc.OuterIntegral(lambda u,p: fcs.Constant(0),"near(x[0],{d_x})".format(d_x=p_d["x"][0]-p_d["margin"]),field_quantity="infg"),
        bc.OuterIntegral(lambda u,p: fcs.Constant(0),"!near(x[0],{d_x})".format(d_x=p_d["x"][0]-p_d["margin"]),field_quantity="infg")
        ]
scan = [
    {

    "parameters":{
        "dummy":0
    },
    "entity_types":[
        {
        "parameters":get_update_dict(Tn,{"R_infg":v}),
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

    } for v in [1*p["R_l"]]
]

T = range(5)
start = time.process_time()
run({**p,**p_d,**p_sim,**p_boundary},T,1,domainBC,"/extra/kiwitz/parameter_scan_large/",extCache="/extra/kiwitz/parameter_scan_large/",scan=scan)
end = time.process_time()
# if rank == 0:
print("--------------------- total Time: {t} m ---------------------".format(t=(end-start)/60))







