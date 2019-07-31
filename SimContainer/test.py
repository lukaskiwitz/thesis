# -*- coding: utf-8 -*-

import FieldProblem as fp
import Entity
import MySolver
import fenics as fcs
import numpy as np
import matplotlib.pyplot as plt
import BC as bc
import SimContainer as sc
from bcFunctions import cellBC_il2,cellBC_il6,outerBC_il2,outerBC_il6
from copy import copy,deepcopy
import BooleanInternalSolver as intSolver
import json
import os
import time
import random
from scipy.constants import N_A

factor = 1e0
dd = 1.2
p_global= {
         "R_il2":0,
         "R_il6":100*N_A**-1*10e9*factor,
         "q_il6":0,
         "q_il2":0,
         "k_on": 10e-9*111.6/60**2,#111.6 per hour
         "rho": 0.05,#mu
         "D":0.01**2,#mu² per s
         "high":4000*N_A**-1*10e9*factor,
         "low":100*N_A**-1*10e9*factor,
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
#                draw = ran.random()
#                if draw > 3/4:
#                    p_temp["q_il2"] = p_global["q_high"]#two molecules per second
#                    p_temp["R_il2"] = p_temp["low"]
#                elif draw > 2/4:
#                    p_temp["R_il2"] = p_temp["high"]
#                else:
#                    p_temp["R_il2"] = p_temp["low"]
                cell.name = "cell_{xCoord}_{yCoord}_{zCoord}".format(xCoord=x,yCoord=y,zCoord=z)
                cell.p = p_temp
#                cell.addSolver(deepcopy(interSolver))
                cellList.append(cell)
   
    return cellList


domain = Entity.DomainCube([-dd,-dd,-dd],[dd,dd,dd],[
       bc.Integral(outerBC_il2,fieldQuantity="il2"),
       bc.Integral(outerBC_il6,fieldQuantity="il6")
        ])
p_domain = deepcopy(p_global)
p_domain.update({
         "R_il2":p_global["high"],
         "R_il6":0,
         "q_il6":p_global["q_high"],
         "q_il2":0})
domain.p = p_domain

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

sc.initialize(load_subdomain="./cache/boundary_markers.h5")
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


xScale = []
for n,i in enumerate(range(250)):

    xScale.append([n,i])
    updateState(p_global,sc,n)
    start = time.process_time()
    sc.step(1)
    end = time.process_time()
    print("time: "+str(end-start)+"s for step number "+str(n))
    times.append(end-start)
    dump = json.dumps(sc.log())
    with open("./logs/"+str(i),"w") as file:
        file.write(dump)
    sc.saveFields(n)


#with open("./timing","w") as file:
#        file.write(times)
    


