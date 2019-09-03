# -*- coding: utf-8 -*-

import matplotlib as mpl
mpl.use('Agg')

import FieldProblem as fp
import Entity
#import MySolver
import fenics as fcs
import numpy as np
#import matplotlib.pyplot as plt
import BC as bc
import SimContainer as sc
from bcFunctions import cellBC_il2,cellBC_il6,outerBC_il2,outerBC_il6,absorbing_il2
from copy import copy,deepcopy
import BooleanInternalSolver as intSolver
import json
import os
import time
import random
from scipy.constants import N_A

from scipy.constants import N_A
from run import run
from distPlot_dump import distPlot_dump
from distPlot import distPlot

#path = "/extra/results/sim_1/"
#factor = 1e0
def bD(n):
    return n/10

x = np.linspace(-bD(15),bD(15),int(15))
y = np.linspace(-bD(10),bD(10),int(10))
z = np.linspace(-bD(10),bD(10),int(10))


d_x = x[-1]+0.2
d_y = y[-1]+0.2
d_z = z[-1]+0.2

p_d = {"x":x,
       "y":y,
       "z":z,
       "d_x":d_x,
       "d_y":d_y,
       "d_z":d_z
       }
p_sim = {#default values
         "R_il2":0,
         "q_il2":0,
         "R_il6":0,
         "q_il6":0
         }
p_c = {
         "k_on":10e9*111.6/60**2,#111.6 per hour
         "rho": 0.05,#mu
         "D":(10**0.5*0.01)**2,#mu² per s
         "R_h":400*N_A**-1*10e9,
         "R_l":10*N_A**-1*10e9,
         "kd":0,#0.1/(60*2),
         "q_h":10*N_A**-1*10e9,
         "q_l":0,
         "radius":2,
         "N":6
         }

#distPlot_dump("/extra/results/sim_test/","/extra/extCache/")
#distPlot("/extra/results/sim_test/")

#T = range(25)
#
##
#"""secHigh_fLux_fractionLow"""
#p_il2 = {
#     "R_il2_s":p_c["R_h"],#secretors
#     "R_il2_f":p_c["R_l"],#fraction 
#     "R_il2_n":p_c["R_l"],#normal
#     "R_il2_b":p_c["R_h"],#boundary
#     "q_il2_s":p_c["q_h"],
#     "q_il2_f":p_c["q_l"],
#     "q_il2_n":p_c["q_l"],
#     "q_il2_b":p_c["q_l"]
#     }
#p_il6 = {
#     "R_il6_s":p_c["R_l"],
#     "R_il6_f":p_c["R_l"],
#     "R_il6_n":0,
#     "R_il6_b":0,
#     "q_il6_s":p_c["q_h"]*0.5,
#     "q_il6_f":p_c["q_l"],
#     "q_il6_n":p_c["q_l"],
#     "q_il6_b":p_c["q_h"],
#     }
#
#domainBC = [
#       bc.outerIntegral(lambda u,p: fcs.Constant(0),"!near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il2"),
##       bc.outerDirichletBC(0,"near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il2"),
#       bc.outerIntegral(outerBC_il2,"near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il2"),
#       bc.outerIntegral(outerBC_il6,"near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il6"),
#       bc.outerIntegral(lambda u,p: fcs.Constant(0),"!near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il6")
#        ]
#
#run({**p_d,**p_c,**p_il2,**p_il6,**p_sim},T,domainBC,"/extra/kiwitz/results_fluxRun/secHigh/",extCache="/extra/kiwitz/extCache/")
#
#"""secLow_flux_fractionHigh"""
#p_il2 = {
#     "R_il2_s":p_c["R_l"],#secretors
#     "R_il2_f":p_c["R_h"],#fraction 
#     "R_il2_n":p_c["R_l"],#normal
#     "R_il2_b":p_c["R_h"],#boundary
#     "q_il2_s":p_c["q_h"],
#     "q_il2_f":p_c["q_l"],
#     "q_il2_n":p_c["q_l"],
#     "q_il2_b":p_c["q_l"]
#     }
#p_il6 = {
#     "R_il6_s":p_c["R_l"],
#     "R_il6_f":p_c["R_l"],
#     "R_il6_n":0,
#     "R_il6_b":0,
#     "q_il6_s":p_c["q_h"]*0.5,
#     "q_il6_f":p_c["q_l"],
#     "q_il6_n":p_c["q_l"],
#     "q_il6_b":p_c["q_h"],
#     }
#
#domainBC = [
#       bc.outerIntegral(lambda u,p: fcs.Constant(0),"!near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il2"),
##       bc.outerDirichletBC(0,"near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il2"),
#       bc.outerIntegral(outerBC_il2,"near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il2"),
#       bc.outerIntegral(outerBC_il6,"near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il6"),
#       bc.outerIntegral(lambda u,p: fcs.Constant(0),"!near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il6")
#        ]
#
#run({**p_d,**p_c,**p_il2,**p_il6,**p_sim},T,domainBC,"/extra/kiwitz/results_fluxRun/secLow/",extCache="/extra/kiwitz/extCache/")


from distPlot_dump_test import distPlot_dump
from distPlot import distPlot
import mpi4py.MPI as MPI
import seaborn as sns

path = "/extra/kiwitz/results_fluxRun/"
extCache = "/extra/kiwitz/extCache/"
fields = ["il2","il6"]

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

l = [
"secHigh",
"secLow"]
#l = l[2:-1]
#l = l[2:3]
threshholds = [
        [0.04,3],
        [0.12,3]
        ]
#for i in l:
#    pass
#    distPlot_dump(path+i+"/",extCache,fields,12)
for i,e in enumerate(l):  
    pass
    distPlot(path+e+"/",threshholds[i],cutoff=[10,10],imgPath="/home/kiwitz/plots_fluxRun/",l=e,il2_limits=(0,0.3),il6_limits=(0,6))






