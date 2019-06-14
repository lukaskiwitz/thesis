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

def makeCellList(p):
    p_center = deepcopy(p)
    p_center["R"] = 10**2
#    p_center["q"] = 10
    cellList = []
    cell_center = Entity.Cell([0,0,0],p_center["rho"],[
                bc.Integral(cellBC,fieldQuantity="C")
                ])
    cell_center.p = p_center
    cellList.append(cell_center)
    interSolver = intSolver.BooleanInternalSolver()
    
    for l,r in enumerate([0.15,0.3]):
        n = 6*(l+1)
        print(np.linspace(0,2*np.pi*(1-1/n),n))
        for i,phi in enumerate(np.linspace(0,2*np.pi*(1-1/n),n)):
            p_temp = deepcopy(p)
            p_temp["R"] = 10**8#10**4 if (i+l) % 2 == 0 else 10**2
            x = r*np.cos(phi)
            y = r*np.sin(phi)
            cell = Entity.Cell([x,y,0],p_temp["rho"],[
                bc.Integral(cellBC,fieldQuantity="C")
                ])
            cell.p = p_temp
            cell.name = "cell_{r}_{phi}".format(r=r,phi=phi)
            cell.addSolver(deepcopy(interSolver))
            cellList.append(cell)
    print(len(cellList))
    return cellList


    
p_global= {
         "R":pow(10,2),
         "q": 10,
         "k_on": 111.6*6.02214076*10**(23-9)/60**2,
         "rho": 0.05,
         "D":0.01**2,
         "L":1,
         "R_resp":pow(10,2),
         "N":100
        }



domain = Entity.DomainSphere([0,0,0],0.5,[
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

for i in makeCellList(p_global):
    sc.addEntity(i)

sc.addField(fieldProblem)

sc.initialize()

if not os.path.isdir("./logs"):    
    os.mkdir("./logs")
    
for i in range(100):
    #sc.entityList[0].p["R"] = 10**(i)
    sc.step(20)
    
    
    
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
    
    


