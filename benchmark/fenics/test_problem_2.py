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

def cellBC(u,p):
        R = p["R"]
        q = p["q"]
        k_on = p["k_on"]

        return (q-u*k_on*R)
def outerBC(u,p):
    k_on = p["k_on"]
    rho = p["rho"]
    L = p["L"]
    N = p["N"]
    R_resp = p["R_resp"]
    
    return k_on*u*(L+rho)*N*R_resp
    
def u_exact(x,p):
    R = p["R"]
    q = p["q"]
    k_on = p["k_on"]
    rho = p["rho"]
    D= p["D"]
    r = np.linalg.norm(x)
    L = p["L"]
    N = p["N"]
    R_resp = p["R_resp"]
    
    
    A = q*rho/(k_on*r)
    num = 4*np.pi*D*r*(L + rho) + k_on*(L-r+rho)*N*R_resp
    denom = k_on*L*R*N*R_resp + 4*np.pi*D*rho*(L+rho)*(R+N*R_resp)
    return A * (num/denom)

def u_exact_str(p):
    R = p["R"]
    q = p["q"]
    k_on = p["k_on"]
    rho = p["rho"]
    D= p["D"]
    L = p["L"]
    N = p["N"]
    R_resp = p["R_resp"]
    
    exp = "{q}*{rho}/({k_on}*(sqrt(pow(x[0],2) + pow(x[1],2))))*((4*{pi}*{D}*(sqrt(pow(x[0],2) + pow(x[1],2)))*({L} + {rho}) + {k_on}*({L}-(sqrt(pow(x[0],2) + pow(x[1],2)))+{rho})*{N}*{R_resp})/({k_on}*{L}*{R}*{N}*{R_resp} + 4*{pi}*{D}*{rho}*({L}+{rho})*({R}+{N}*{R_resp})))"
    exp = exp.format(q=q,rho=rho,k_on=k_on,R=R,D=D,L=L,N=N,R_resp=R_resp,pi=np.pi)
    return exp

#output path
vtkfile = fcs.File("./sol/solution.pvd")
#parameter vector
p = {
         "R":pow(10,1),
         "q": 10,
         "k_on": 111.6,
         "rho": 0.1,
         "D":10,
         "L":1,
         "R_resp":pow(10,1),
         "N":100
        }


def run(res,options):
    """Mesh Generation"""
    #boxBC = {"D":u_exact_str(p)}#bcDict; set dirichlet on the outer boundary to exact value
    boxBC = {"Rec":outerBC}
    domRad = p["L"]+p["rho"]#radius of the domain
    box = entity.OuterSphere([0,0,0],domRad,bcDict=boxBC)#instantiates SubDomain at origin;
    
    meshGen = myMesh.MeshGenerator(outerBoundary=box)#instantiates meshgenerator object
    cellList=[{"center":[0,0,0],"radius":p["rho"],"bcDict":{"Rec":cellBC}}]#defines list of cells(as SubDomain objects) in simulation
    meshGen.cellList= cellList
    meshGen.dim = 3#explicitly sets mesh dimension
    
    meshGen.compileSubdomains()
    #print(meshGen.domain)
    mesh,subdomains, boundary_markers = meshGen.meshGen(res)#returns mesh, subdomains and boundary_markers

    
    """Solver Setup"""
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
    
    y_sol = []
    x = np.linspace(p["rho"]*(1+0.001),domRad*(1-0.01),100)
    for i in x:
        y_sol.append(u([i,0,0]))

    #vtkfile << u#writes solution to file
    #plt.plot(x,y_sol)
    return end-start,[x,y_sol]



dT,sol = run(64,"")

x = np.linspace(p["rho"],p["L"]+p["rho"],100)

y = [u_exact(i,p) for i in x]

plt.plot(x,y)
plt.plot(sol[0],sol[1])
plt.legend(["exact","FEM"])
plt.xlabel("r")
plt.ylabel("u")
plt.title("problem 2")
#plt.savefig("problem_2.png",dpi=200)   
        







