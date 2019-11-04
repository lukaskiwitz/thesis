#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 14:30:17 2019

@author: lukas kiwitz
"""

import MyEntities as entity
import myMeshing as myMesh
import PoissonSolver as mySolver
import fenics as fcs
import dolfin as dlf
from copy import copy,deepcopy

import time
import numpy as np
import matplotlib.pyplot as plt



def cellBC(u,p):
    """
    Defines the flux boundary conditions for the cell.
    Can be passed to the solver in the "Rec" field of the "bcDict" dictionary
    """
    R = p["R"]
    q = p["q"]
    k_on = p["k_on"]
    print(R)
    return (q-u*k_on*R)
def outerBC(u,p):
    """
    Defines the flux boundary condition on the outer boundary.
    Can be passed to the solver in the "Rec" field of the "bcDict" dictionary
    """
    k_on = p["k_on"]
    rho = p["rho"]
    L = p["L"]
    N = p["N"]
    R_resp = p["R_resp"]
    
    return k_on*u*(L+rho)*N*R_resp
    
#output path
vtkfile = fcs.File("./sol/solution.pvd")
#parameter vector
p_global= {
         "R":pow(10,2),
         "q": 10,
         "k_on": 111.6,
         "rho": 0.1,
         "D":10,
         "L":1,
         "R_resp":pow(10,2),
         "N":100
        }

def makeCellList(p):
    cellList = []
    for l,r in enumerate([0.4,0.8]):
        for i,phi in enumerate(np.linspace(0,2*np.pi,6)):
            p_temp = deepcopy(p)
            p_temp["R"] = 10**4 if (i+l) % 2 == 0 else 10**2
            x = r*np.cos(phi)
            y = r*np.sin(phi)
            entry = {"center":[x,y,0],"radius":p_temp["rho"],"bcDict":{"Rec":cellBC,"p":p_temp}}
            cellList.append(entry)
    return cellList
def run(res,p):
    """Mesh Generation"""
    
    boxBC = {"Rec":outerBC}
    domRad = p["L"]+p["rho"]#radius of the domain
    box = entity.OuterSphere([0,0,0],domRad,bcDict=boxBC)#instantiates SubDomain at origin;
    
    meshGen = myMesh.MeshGenerator(outerBoundary=box)#instantiates meshgenerator object
    cellList=makeCellList(p)#defines list of cells(as SubDomain objects) in simulation
    meshGen.cellList= cellList
    meshGen.dim = 3#explicitly sets mesh dimension
#    
    meshGen.compileSubdomains()
#    #print(meshGen.domain)
    mesh,subdomains, boundary_markers = meshGen.meshGen(res)#returns mesh, subdomains and boundary_markers
#    
    file=dlf.File("mesh_64.xml")
    file<< mesh
#    
#    """Solver Setup"""
    solver = mySolver.PoissonSolver(mesh,subdomains, boundary_markers,3)#instantiates my solver class
    solver.p = p # set modell parameters
    solver.compileSolver()
    
    #alters parameters of fenics solver
    solver.solver.parameters["newton_solver"]["linear_solver"] = "gmres"
    solver.solver.parameters["newton_solver"]["preconditioner"] = "ilu"
        
    
#    u_int = fcs.interpolate(fcs.Expression(u_exact_str(p),degree=1),solver.V)
#    y_sol = []
#    x = np.linspace(p["rho"]*(1+0.01),domRad*(1-0.01),100)
#    for i in x:
#        y_sol.append(u_int([i,0,0]))
#    return 1,[x,y_sol]

    """runs simulation"""
    start = time.process_time_ns()
    u = solver.solve()#computes solution
    end = time.process_time_ns()
    vtkfile << u#writes solution to file
#    y_sol = []
#    x = np.linspace(p["rho"]*(1+0.001),domRad*(1-0.01),100)
#    for i in x:
#        y_sol.append(u([i,0,0]))
#
#    #plt.plot(x,y_sol)
#    return end-start,[x,y_sol]



run(64,p_global)







