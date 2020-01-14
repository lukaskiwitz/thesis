# -*- coding: utf-8 -*-

import matplotlib as mpl
mpl.use('Agg')



#import MySolver
import fenics as fcs
import numpy as np
#import matplotlib.pyplot as plt
import BC as bc
from scipy.constants import N_A
from mark_run import run
from bcFunctions import outerBC_il2,outerBC_il6
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



T = [0]

#
"""secHigh_fLux_fractionLow"""
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
     "R_il6_s":p_c["R_l"],
     "R_il6_f":p_c["R_l"],
     "R_il6_n":0,
     "R_il6_b":0,
     "q_il6_s":p_c["q_h"]*0.5,
     "q_il6_f":p_c["q_l"],
     "q_il6_n":p_c["q_l"],
     "q_il6_b":p_c["q_h"],
     }

domainBC = [
       bc.outerIntegral(lambda u,p: fcs.Constant(0),"!near(x[0],-{d_x})".format(d_x=d_x),field_quantity="il2"),
#       bc.outerDirichletBC(0,"near(x[0],-{d_x})".format(d_x=d_x),field_quantity="il2"),
       bc.outerIntegral(outerBC_il2,"near(x[0],-{d_x})".format(d_x=d_x),field_quantity="il2"),
       bc.outerIntegral(outerBC_il6,"near(x[0],-{d_x})".format(d_x=d_x),field_quantity="il6"),
       bc.outerIntegral(lambda u,p: fcs.Constant(0),"!near(x[0],-{d_x})".format(d_x=d_x),field_quantity="il6")
        ]

run({**p_d,**p_c,**p_il2,**p_il6,**p_sim},T,domainBC,"/home/kiwitz/markings/",extCache="/extra/kiwitz/extCache/")