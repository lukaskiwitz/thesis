# -*- coding: utf-8 -*-

import FieldProblem as fp
import Entity
#import MySolver
import fenics as fcs
import numpy as np
#import matplotlib.pyplot as plt
import BC as bc
<<<<<<< HEAD
#import SimContainer as sc
from bcFunctions import cellBC_il2,cellBC_il6,outerBC_il2,outerBC_il6,absorbing_il2
from copy import deepcopy
=======
import SimContainer as sc
from bcFunctions import cellBC_il2,cellBC_il6,outerBC_il2,outerBC_il6,absorbing_il2
from copy import copy,deepcopy
import BooleanInternalSolver as intSolver
import json
import os
import time
import random
from scipy.constants import N_A
>>>>>>> 346cb68085eb45681e048f422269d86c3476166d

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
<<<<<<< HEAD
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

T = range(1)

"""secHigh_dirichlet_fractionLow"""
p_il2 = {
     "R_il2_s":p_c["R_l"],#secretors
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
     "q_il6_s":p_c["q_h"]*0.5,
     "q_il6_f":p_c["q_l"]*1,
     "q_il6_n":p_c["q_l"]*1,
     "q_il6_b":p_c["q_h"],
     }

domainBC = [
       bc.outerIntegral(lambda u,p: fcs.Constant(0),"!near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il2"),
#       bc.outerDirichletBC(0,"near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il2"),
       bc.outerIntegral(outerBC_il2,"near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il2"),
       bc.outerIntegral(outerBC_il6,"near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il6"),
       bc.outerIntegral(lambda u,p: fcs.Constant(0),"!near(x[0],-{d_x})".format(d_x=d_x),fieldQuantity="il6")
        ]

run({**p_d,**p_c,**p_il2,**p_il6,**p_sim},T,domainBC,"/extra/kiwitz/results/test_kd/",extCache="/extra/kiwitz/extCache/")

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
#     "q_il6_s":p_c["q_h"]*0.5,
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
#     "q_il6_s":p_c["q_h"]*0.5,
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
#run({**p_d,**p_c,**p_il2,**p_il6,**p_sim},T,domainBC,"/extra/kiwitz/results/secLow_dirichlet_fractionHigh/",extCache="/extra/kiwitz/extCache/")
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
#     "q_il6_s":p_c["q_h"]*0.5,
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
#run({**p_d,**p_c,**p_il2,**p_il6,**p_sim},T,domainBC,"/extra/kiwitz/results/secHigh_noFLux_fractionLow/",extCache="/extra/kiwitz/extCache/")
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
#    "q_il6_s":p_c["q_h"]*0.5,
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
#run({**p_d,**p_c,**p_il2,**p_il6,**p_sim},T,domainBC,"/extra/kiwitz/results/secLow_noFlux_fractionHigh/",extCache="/extra/kiwitz/extCache/")
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
#run({**p_d,**p_c,**p_il2,**p_il6,**p_sim},T,domainBC,"/extra/kiwitz/results/secHigh_fLux_fractionLow/",extCache="/extra/kiwitz/extCache/")
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
#run({**p_d,**p_c,**p_il2,**p_il6,**p_sim},T,domainBC,"/extra/kiwitz/results/secLow_flux_fractionHigh/",extCache="/extra/kiwitz/extCache/")
#
#
#from distPlot_dump_test import distPlot_dump
#from distPlot import distPlot
#import mpi4py.MPI as MPI
#import seaborn as sns
#
#path = "/extra/kiwitz/results/"
#extCache = "/extra/kiwitz/extCache/"
#fields = ["il2","il6"]
#
#comm = MPI.COMM_WORLD
#rank = comm.Get_rank()
#size = comm.Get_size()
#
#l = ["secHigh_dirichlet_fractionLow",
#"secLow_dirichlet_fractionHigh",
#"secHigh_noFLux_fractionLow",
#"secLow_noFlux_fractionHigh",
#"secHigh_fLux_fractionLow",
#"secLow_flux_fractionHigh"]
##l = l[2:-1]
##l = l[2:3]
#threshholds = [
#        [0.2,0.5]
#        ]
#for i in l:
#    pass
#    distPlot_dump(path+i+"/",extCache,fields,12)
##for i in l:  
##    pass
##    distPlot(path+i+"/",threshholds[0],cutoff=[10,10],imgPath="/home/kiwitz/plots/",l=i)
##distPlot("/extra/kiwitz/test/data/",threshholds[0])

=======
         "R_il6":10*N_A**-1*10e9*factor,
         "q_il6":0,
         "q_il2":0,
         "k_on": 10e9*111.6/60**2,#111.6 per hour
         "rho": 0.05,#mu
         "D":0.1**2,#mu² per s
         "high":400*N_A**-1*10e9*factor,
         "low":10*N_A**-1*10e9*factor,
         "kd":0.1/(60*2),
         "q_high":2*N_A**-1*10e9*factor,
         "dd":dd
        }

    
def updateState(p,sc,t):
    ran = random.Random()
    ran.seed(t)
    for i,e in enumerate(sc.entityList):
        p_temp = deepcopy(p)
        draw = ran.random()
        if draw > 3/4:
            p_temp["q_il2"] = p_global["q_high"]#two molecules per second
            p_temp["R_il2"] = p_temp["low"]
        elif draw > 2/4:
            p_temp["R_il2"] = p_temp["high"]
        else:
            p_temp["R_il2"] = p_temp["low"]
            
        draw = ran.random()
        if draw > 3/4:
            pass#p_temp["q_il6"] = p_global["q_high"]#two molecules per second
        e.p = p_temp
        
def makeCellListGrid(p,xLine,yLine,zLine):
    
    interSolver = intSolver.BooleanInternalSolver()
    cellList = []
    ran = random.Random()
    ran.seed(1)
    for x in xLine:
        for y in yLine:
            for z in zLine:                    
                p_temp = deepcopy(p)
                cell = Entity.Cell([x,y,z],p_temp["rho"],
               [
                       bc.Integral(cellBC_il2,fieldQuantity="il2"),
                       bc.Integral(cellBC_il6,fieldQuantity="il6")
                ])
                cell.name = "cell_{xCoord}_{yCoord}_{zCoord}".format(xCoord=x,yCoord=y,zCoord=z)
                cell.p = p_temp
#                cell.addSolver(deepcopy(interSolver))
                cellList.append(cell)
   
    return cellList
p_domain = deepcopy(p_global)
p_domain.update({
         "R_il2":p_global["high"],
         "R_il6":0,
         "q_il6":p_global["q_high"],
         "q_il2":0})

domain = Entity.DomainCube([-dd,-dd,-dd],[dd,dd,dd],p_domain,[
       bc.outerIntegral(lambda u,p: fcs.Constant(0),"!near(x[0],-{dd}) && on_boundary".format(dd=dd),fieldQuantity="il2"),
#       bc.outerDirichletBC(0,"near(x[0],-{dd}) && on_boundary".format(dd=dd),fieldQuantity="il2"),
       bc.outerIntegral(outerBC_il2,"near(x[0],-{dd}) && on_boundary".format(dd=dd),fieldQuantity="il2"),
       bc.outerIntegral(outerBC_il6,"near(x[0],-{dd}) && on_boundary".format(dd=dd),fieldQuantity="il6"),
       bc.outerIntegral(lambda u,p: fcs.Constant(0),"!near(x[0],-{dd}) && on_boundary".format(dd=dd),fieldQuantity="il6")
        ])


"""IL-2"""
solver_il2 = MySolver.PoissonSolver()

fieldProblem_il2 = fp.FieldProblem()
fieldProblem_il2.fieldName = "il2"
fieldProblem_il2.fieldQuantity = "il2"


fieldProblem_il2.setSolver(solver_il2)
fieldProblem_il2.p = deepcopy(p_global)
fieldProblem_il2.meshCached = "./cache/meshCache_il2"
fieldProblem_il2.setOuterDomain(domain)

"""IL-6"""
solver_il6 = MySolver.PoissonSolver()

fieldProblem_il6 = fp.FieldProblem()
fieldProblem_il6.fieldName = "il6"
fieldProblem_il6.fieldQuantity = "il6"


fieldProblem_il6.setSolver(solver_il6)
fieldProblem_il6.p = deepcopy(p_global)
fieldProblem_il6.meshCached = "./cache/meshCache_il2"
fieldProblem_il6.setOuterDomain(domain)

sc = sc.SimContainer()
boxDim = 1

n = 2*boxDim/0.2
#n = 3

x = np.linspace(-boxDim,boxDim,int(n))
y = np.linspace(-boxDim,boxDim,int(n))
z = np.linspace(-boxDim,boxDim,int(n))

#x = [0]
#y = x
#z = x

for i in makeCellListGrid(p_global,x,y,z):
    sc.addEntity(i)

sc.addField(fieldProblem_il2)
#sc.addField(fieldProblem_il6)


sc.initialize(load_subdomain="./cache/boundary_markers_il2.h5")

#sc.initialize()
print("init complete")
if not os.path.isdir("./logs"):
    try:
        os.mkdir("./logs")
    except:
        pass
else:
    for file in os.listdir("./logs"):
        try:
            os.remove("./logs/"+str(file))
        except:
            pass
times = []

#sc.saveSubdomains()
#print("subdomains saved")
#sc.saveDomain()
#print("domain saved")
>>>>>>> 346cb68085eb45681e048f422269d86c3476166d



<<<<<<< HEAD
=======
for n,i in enumerate(range(10)):
>>>>>>> 346cb68085eb45681e048f422269d86c3476166d






