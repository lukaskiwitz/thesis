# -*- coding: utf-8 -*-

import FieldProblem as fp
import Entity
import MySolver
import fenics as fcs
import numpy as np
import matplotlib.pyplot as plt
import BC as bc
import SimContainer as sc
from bcFunctions import cellBC,outerBC
from copy import copy,deepcopy
import BooleanInternalSolver as intSolver
import json
import os
import time

def makeCellListSphere(p):
    p_center = deepcopy(p)
    p_center["R"] = 10**2
#    p_center["q"] = 10
    cellList = []
    cell_center = Entity.Cell([0,0,0],p_center["rho"],[
                bc.Integral(cellBC,fieldQuantity="C")
                ])
    cell_center.p = p_center
    cellList.append(cell_center)
#    interSolver = intSolver.BooleanInternalSolver()
    
#    for l,r in enumerate([0.15,0.3]):
#        n = 6*(l+1)
#        print(np.linspace(0,2*np.pi*(1-1/n),n))
#        for i,phi in enumerate(np.linspace(0,2*np.pi*(1-1/n),n)):
#            p_temp = deepcopy(p)
#            p_temp["R"] = 10**8#10**4 if (i+l) % 2 == 0 else 10**2
#            x = r*np.cos(phi)
#            y = r*np.sin(phi)
#            cell = Entity.Cell([x,y,0],p_temp["rho"],[
#                bc.Integral(cellBC,fieldQuantity="C")
#                ])
#            cell.p = p_temp
#            cell.name = "cell_{r}_{phi}".format(r=r,phi=phi)
#            cell.addSolver(deepcopy(interSolver))
#            cellList.append(cell)
#    print(len(cellList))
    return cellList

def makeCellListGrid(p,xLine,yLine,zLine):
    
    interSolver = intSolver.BooleanInternalSolver()
    cellList = []
    for x in xLine:
        for y in yLine:
            for z in zLine:                    
                p_temp = deepcopy(p)
                cell = Entity.Cell([x,y,z],p_temp["rho"],[
                bc.Integral(cellBC,fieldQuantity="C")
                ])
                if x == 0 and y == 0 and z == 0:
                    p_temp["R"] = p_global["low"]
                    cell.name = "center"
                cell.p = p_temp
                cell.addSolver(deepcopy(interSolver))
                cellList.append(cell)
   
    return cellList
factor = 1e+12
p_global= {
         "R":6.642156e-12*factor,
         "q": 1.6605391e-14*factor,
         "k_on": (0.031**(-1))*factor,
         "rho": 0.05,
         "D":0.01**2,
         "high":6.642156e-12*factor,
         "low":1.6605391e-13*factor,
         "threshold":0.0001,
         "decay":0.01,
         "L":5,
         "R_resp":pow(10,2),
         "N":100
        }



#domain = Entity.DomainSphere([0,0,0],0.5,[
#        bc.Integral(outerBC,fieldQuantity="C")
#        ])
domain = Entity.DomainCube([-0.5,-0.5,-0.3],[0.5,0.5,0.3],[
        bc.Integral(outerBC,fieldQuantity="C")
        ])
domain.p = p_global


solver = MySolver.PoissonSolver()

fieldProblem = fp.FieldProblem()
fieldProblem.fieldName = "cytokine"
fieldProblem.fieldQuantity = "C"
vtkfile = fcs.File("./sol/solution.pvd")

fieldProblem.setSolver(solver)
fieldProblem.p = p_global
fieldProblem.meshCached = "./cache/meshCache_64_cytokine.xml"
fieldProblem.res = 64
fieldProblem.setOuterDomain(domain)

sc = sc.SimContainer()

x = np.linspace(-0.2,0.2,3)
y = np.linspace(-0.2,0.2,3)
z = np.array([0])#np.linspace(-0.2,0.2,2)

for i in makeCellListGrid(p_global,x,y,z):
    sc.addEntity(i)

sc.addField(fieldProblem)

sc.initialize()

if not os.path.isdir("./logs"):    
    os.mkdir("./logs")

times = []
for n,i in enumerate(np.linspace(p_global["low"],p_global["high"],10)):
    start = time.process_time()
    sc.step(0)
    end = time.process_time()
    
    print("time: "+str(end-start)+"s for step number "+str(n))
    times.append(end-start)
    sc.getEntityByName("center").p["R"] = i
    print("changing entity: "+str(sc.getEntityByName("center").name)+" R="+str(i))
    
    
    dump = json.dumps(sc.log())
    with open("./logs/"+str(i),"w") as file:
        file.write(dump)
    u = sc.fields[0].getFields()
    
#    mesh = solver.mesh
#    V_vec = fcs.FunctionSpace(mesh,"CG",1)
#    n = fcs.FacetNormal(mesh)
#    
#    ds = fcs.Measure("dS", domain=mesh)
#    flux = fcs.project(fcs.div(fcs.grad(u)),V_vec)
#    
##    flux_n = -1*fcs.assemble(flux)
    
    vtkfile << u
    
    


