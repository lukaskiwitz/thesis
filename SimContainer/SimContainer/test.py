# -*- coding: utf-8 -*-

import FieldProblem as fp
import Entity
#import MySolver
import fenics as fcs
import numpy as np
#import matplotlib.pyplot as plt
import BC as bc
#import SimContainer as sc
from bcFunctions import cellBC_il2,cellBC_il6,outerBC_il2,outerBC_il6,absorbing_il2
from copy import deepcopy

from scipy.constants import N_A
from run import run
#from distPlot_dump import distPlot_dump
#from distPlot import distPlot

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
         "k_on": 10e9*111.6/60**2,#111.6 per hour
         "rho": 0.05,#mu
         "D":(10**0.5*0.01)**2,#mu² per s
         "R_h":400*N_A**-1*10e9,
         "R_l":10*N_A**-1*10e9,
         "kd":0.1/(60*2),
         "q_h":10*N_A**-1*10e9,
         "q_l":0,
         "radius":2,
         "N":6
         }

#distPlot_dump("/extra/results/sim_test/","/extra/extCache/")
#distPlot("/extra/results/sim_test/")

T = range(1)

"""secHigh_dirichlet_fractionLow"""
p_il2 = {
     "R_il2_s":p_c["R_h"],#secretors
     "R_il2_f":p_c["R_l"],#fraction 
     "R_il2_n":p_c["R_l"],#normal
     "R_il2_b":p_c["R_h"],#boundary
     "q_il2_s":p_c["q_h"],
     "q_il2_f":p_c["q_l"],
     "q_il2_n":p_c["q_l"],
     "q_il2_b":p_c["q_l"]
     }
p_il6 = {
     "R_il6_s":p_c["R_h"],
     "R_il6_f":p_c["R_l"],
     "R_il6_n":0,
     "R_il6_b":0,
     "q_il6_s":p_c["q_h"]*1,
     "q_il6_f":p_c["q_l"]*1,
     "q_il6_n":p_c["q_l"]*1,
     "q_il6_b":p_c["q_h"],
     }

domainBC = [
       bc.outerIntegral(lambda u,p: fcs.Constant(0),"!near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il2"),
       bc.outerDirichletBC(0,"near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il2"),
#       bc.outerIntegral(outerBC_il2,"near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il2"),
       bc.outerIntegral(outerBC_il6,"near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il6"),
       bc.outerIntegral(lambda u,p: fcs.Constant(0),"!near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il6")
        ]

run({**p_d,**p_c,**p_il2,**p_il6,**p_sim},T,domainBC,"/extra/kiwitz/results/bc_Test/",extCache="/extra/kiwitz/extCache/")


#T = range(100)
#
#"""secHigh_dirichlet_fractionLow"""
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
#     "R_il6_s":p_c["R_h"],
#     "R_il6_f":p_c["R_l"],
#     "R_il6_n":0,
#     "R_il6_b":0,
#     "q_il6_s":p_c["q_h"]*1,
#     "q_il6_f":p_c["q_l"]*1,
#     "q_il6_n":p_c["q_l"]*1,
#     "q_il6_b":p_c["q_h"],
#     }
#
#domainBC = [
#       bc.outerIntegral(lambda u,p: fcs.Constant(0),"!near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il2"),
#       bc.outerDirichletBC(0,"near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il2"),
##       bc.outerIntegral(outerBC_il2,"near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il2"),
#       bc.outerIntegral(outerBC_il6,"near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il6"),
#       bc.outerIntegral(lambda u,p: fcs.Constant(0),"!near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il6")
#        ]
#
#run({**p_d,**p_c,**p_il2,**p_il6,**p_sim},T,domainBC,"/extra/kiwitz/results/secHigh_dirichlet_fractionLow/",extCache="/extra/kiwitz/extCache/")

#"""secLow_dirichlet_fractionHigh"""
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
#     "q_il6_s":p_c["q_h"],
#     "q_il6_f":p_c["q_l"],
#     "q_il6_n":p_c["q_l"],
#     "q_il6_b":p_c["q_h"],
#     }
#
#domainBC = [
#       bc.outerIntegral(lambda u,p: fcs.Constant(0),"!near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il2"),
#       bc.outerDirichletBC(0,"near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il2"),
##       bc.outerIntegral(outerBC_il2,"near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il2"),
#       bc.outerIntegral(outerBC_il6,"near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il6"),
#       bc.outerIntegral(lambda u,p: fcs.Constant(0),"!near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il6")
#        ]
#
#run({**p_d,**p_c,**p_il2,**p_il6,**p_sim},T,domainBC,"/extra/kiwitz/secLow_dirichlet_fractionHigh/",extCache="/extra/kiwitz/extCache/")
#
#"""secHigh_noFLux_fractionLow"""
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
#     "q_il6_s":p_c["q_h"],
#     "q_il6_f":p_c["q_l"],
#     "q_il6_n":p_c["q_l"],
#     "q_il6_b":p_c["q_h"],
#     }
#
#domainBC = [
#       bc.outerIntegral(lambda u,p: fcs.Constant(0),"!near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il2"),
##       bc.outerDirichletBC(0,"near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il2"),
#       bc.outerIntegral(lambda u,p: fcs.Constant(0),"near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il2"),
#       bc.outerIntegral(outerBC_il6,"near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il6"),
#       bc.outerIntegral(lambda u,p: fcs.Constant(0),"!near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il6")
#        ]
#
#run({**p_d,**p_c,**p_il2,**p_il6,**p_sim},T,domainBC,"/extra/kiwitz/secHigh_noFLux_fractionLow/",extCache="/extra/kiwitz/extCache/")
#
#"""secLow_noFlux_fractionHigh"""
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
#     "q_il6_s":p_c["q_h"],
#     "q_il6_f":p_c["q_l"],
#     "q_il6_n":p_c["q_l"],
#     "q_il6_b":p_c["q_h"],
#     }
#
#domainBC = [
#       bc.outerIntegral(lambda u,p: fcs.Constant(0),"!near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il2"),
##       bc.outerDirichletBC(0,"near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il2"),
#       bc.outerIntegral(lambda u,p: fcs.Constant(0),"near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il2"),
#       bc.outerIntegral(outerBC_il6,"near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il6"),
#       bc.outerIntegral(lambda u,p: fcs.Constant(0),"!near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il6")
#        ]
#
#run({**p_d,**p_c,**p_il2,**p_il6,**p_sim},T,domainBC,"/extra/kiwitz/secLow_noFlux_fractionHigh/",extCache="/extra/kiwitz/extCache/")
#
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
#     "q_il6_s":p_c["q_h"],
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
#run({**p_d,**p_c,**p_il2,**p_il6,**p_sim},T,domainBC,"/extra/kiwitz/secHigh_fLux_fractionLow/",extCache="/extra/kiwitz/extCache/")
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
#     "q_il6_s":p_c["q_h"],
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
#run({**p_d,**p_c,**p_il2,**p_il6,**p_sim},T,domainBC,"/extra/kiwitz/secLow_flux_fractionHigh/",extCache="/extra/kiwitz/extCache/")











