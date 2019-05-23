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

import time
import numpy as np
import matplotlib.pyplot as plt
import gc


    
def u_exact(x,p):
    R = p["R"]
    q = p["q"]
    k_on = p["k_on"]
    rho = p["rho"]
    D= p["D"]
    
    return q*rho/(np.linalg.norm(x)*(k_on*R+D*4*np.pi*rho))

def u_exact_str(p):
    R = str(p["R"])
    q = str(p["q"])
    k_on = str(p["k_on"])
    rho = str(p["rho"])
    D= str(p["D"])
    
    exp = "({q}*{rho}/(sqrt(pow(x[0],2)+pow(x[1],2))*({k_on}*{R}+{D}*4*"+str(np.pi)+"*{rho})))"
    exp = exp.format(q=q,rho=rho,k_on=k_on,R=R,D=D)
    return exp

#output path
vtkfile = fcs.File("./sol/solution.pvd")
#parameter vector
p = {
         "R":pow(10,2),
         "q": 10,
         "k_on": 111.6,
         "rho": 0.5,
         "D":10
        }


def run(res,options):
   
    """Mesh Generation"""
    boxBC = {"D":u_exact_str(p)}#bcDict; set dirichlet on the outer boundary to exact value
    domRad = 1#radius of the domain
    box = entity.OuterSphere([0,0,0],domRad,bcDict=boxBC)#instantiates SubDomain at origin;
    
    meshGen = myMesh.MeshGenerator(outerBoundary=box)#instantiates meshgenerator object
    cellList=[{"center":[0,0,0],"radius":p["rho"],"bcDict":{"Rec":""}}]#defines list of cells(as SubDomain objects) in simulation
    meshGen.cellList= cellList
    meshGen.dim = 3#explicitly sets mesh dimension
    
    meshGen.compileSubdomains()
    #print(meshGen.domain)
    mesh,subdomains, boundary_markers = meshGen.meshGen(t)#returns mesh, subdomains and boundary_markers

    
    """Solver Setup"""
    solver = mySolver.PoissonSolver(mesh,subdomains, boundary_markers,3)#instantiates my solver class
    solver.p = p # set modell parameters
    solver.compileSolver()
    
    #alters parameters of fenics solver
    solver.solver.parameters["newton_solver"]["linear_solver"] = "lu"
    #solver.solver.parameters["newton_solver"]["preconditioner"] = "ilu"
        
    """runs simulation"""
    start = time.process_time_ns()
    u = solver.solve()#computes solution
    end = time.process_time_ns()
    
    y_sol = []
    x = np.linspace(p["rho"]*(1+0.01),domRad*(1-0.01),100)
    for i in x:
        y_sol.append(u([i,0,0]))

    vtkfile << u#writes solution to file
    #plt.plot(x,y_sol)
    return end-start,[x,y_sol]



timing = []
for t in [16,32,64,128]:
    solution = []
    try:
        dT,u = run(t,"")
        timing.append([t,dT])
    except:
        print("Failed lu")
    
np.save("lu",np.array(timing))
   
        
        







