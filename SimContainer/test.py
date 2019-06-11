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

def makeCellList(p):
    p_center = deepcopy(p)
    p_center["R"] = 10**2 
    cellList = []
    cell_center = Entity.Cell([0,0,0],p_center["rho"],[
                bc.Integral(cellBC,fieldQuantity="C")
                ])
    cell_center.p = p_center
    cellList.append(cell_center)
    interSolver = intSolver.BooleanInternalSolver()
    
    for l,r in enumerate([0.2,0.4,0.6]):
        n = 6*(l+1)
        print(np.linspace(0,2*np.pi*(1-1/n),n))
        for i,phi in enumerate(np.linspace(0,2*np.pi*(1-1/n),n)):
            p_temp = deepcopy(p)
            p_temp["R"] = 10**4#10**4 if (i+l) % 2 == 0 else 10**2
            x = r*np.cos(phi)
            y = r*np.sin(phi)
            cell = Entity.Cell([x,y,0],p_temp["rho"],[
                bc.Integral(cellBC,fieldQuantity="C")
                ])
            cell.p = p_temp
            cell.addSolver(deepcopy(interSolver))
            cellList.append(cell)
    print(len(cellList))
    return cellList

p_global= {
         "R":pow(10,2),
         "q": 10,
         "k_on": 111.6,
         "rho": 0.05,
         "D":10,
         "L":1,
         "R_resp":pow(10,2),
         "N":100
        }

#cell = Entity.Cell([0,0,0],p_global["rho"],[
#        bc.Integral(cellBC,fieldQuantity="C")
#        ])
#cell.p = p_global


domain = Entity.DomainSphere([0,0,0],1,[
        bc.DirichletBC("0",fieldQuantity="C")
        ])
domain.p = p_global


solver = MySolver.PoissonSolver()

fieldProblem = fp.FieldProblem()
fieldProblem.fieldName = "cytokine"
fieldProblem.fieldQuantity = "C"
vtkfile = fcs.File("./sol/solution.pvd")

fieldProblem.setSolver(solver)
fieldProblem.p = p_global
fieldProblem.meshCached = "./cache/meshCache_90_cytokine.xml"
fieldProblem.res = 90
fieldProblem.setOuterDomain(domain)

sc = sc.SimContainer()

for i in makeCellList(p_global):
    sc.addEntity(i)

sc.addField(fieldProblem)

sc.initialize()

for i in range(50):
    sc.step(25)
    u = sc.fields[0].getFields()
    vtkfile << u
    
    


