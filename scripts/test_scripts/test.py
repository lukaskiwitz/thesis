# -*- coding: utf-8 -*-



import BC as bc
import fenics as fcs
import numpy as np
from bcFunctions import outerBC_il2, outerBC_il6
from cell_types import p, p_d, p_sim
import time

from run_scan import run

import mpi4py.MPI as MPI

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


domainBC = [
        bc.outerIntegral(lambda u,p: fcs.Constant(0),"!near(x[0],{d_x})".format(d_x=p_d["x"][0]-p_d["margin"]),fieldQuantity="il2"),
        bc.outerIntegral(outerBC_il2,"near(x[0],{d_x})".format(d_x=p_d["x"][0]-p_d["margin"]),fieldQuantity="il2"),
        bc.outerIntegral(outerBC_il6,"near(x[0],{d_x})".format(d_x=p_d["x"][0]-p_d["margin"]),fieldQuantity="il6"),
        bc.outerIntegral(lambda u,p: fcs.Constant(0),"!near(x[0],{d_x})".format(d_x=p_d["x"][0]-p_d["margin"]),fieldQuantity="il6"),
        bc.outerIntegral(lambda u,p: fcs.Constant(0),"near(x[0],{d_x})".format(d_x=p_d["x"][0]-p_d["margin"]),fieldQuantity="infg"),
        bc.outerIntegral(lambda u,p: fcs.Constant(0),"!near(x[0],{d_x})".format(d_x=p_d["x"][0]-p_d["margin"]),fieldQuantity="infg")
        ]
scan = [{"dummy":0}]
T = range(30)
start = time.process_time()
run({**p,**p_d,**p_sim,**p_boundary},T,1,domainBC,"/extra/kiwitz/test_results/",extCache="/extra/kiwitz/extCache_test/",scan=scan)
end = time.process_time()
# if rank == 0:
print("--------------------- total Time: {t} m ---------------------".format(t=(end-start)/60))







