#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 20 15:27:17 2019

@author: lukas
"""


from __future__ import print_function
from fenics import *
import numpy as np
from mshr import *
from math import sin, cos, pi,sqrt
from dolfin import *
from time import process_time_ns
import matplotlib.pyplot as plt


from util import *

def u_exact(x,p):
    R = p["R"]
    q = p["q"]
    k_on = p["k_on"]
    rho = p["rho"]
    D= p["D"]
    
    return q*rho/(np.linalg.norm(x)*(k_on*R+D*4*np.pi*rho))

def u_exact_exp(p):
    R = str(p["R"])
    q = str(p["q"])
    k_on = str(p["k_on"])
    rho = str(p["rho"])
    D= str(p["D"])
    
    exp = "({q}*{rho}/(sqrt(pow(x[0],2)+pow(x[1],2))*({k_on}*{R}+{D}*4*"+str(np.pi)+"*{rho})))"
    exp = exp.format(q=q,rho=rho,k_on=k_on,R=R,D=D)
    return Expression(exp,degree=1)
            
def mySolve(mesh, boundary_markers,subdomains,p):
    neumann = []
    dirichlet = []
    nonLinBC = []
    
    ds = Measure("ds", domain=mesh, subdomain_data=boundary_markers)
    
    V = FunctionSpace(mesh, 'P', 2)
    u = Function(V)
    v = TestFunction(V)
    f = Constant(0)
    D = Constant(p["D"])
    
    def q(u,p):
        R = p["R"]
        q = p["q"]
        k_on = p["k_on"]
        
        return q - k_on*u*R
    
    for i in subdomains:
        if "D" in i.bcDict:
            value = i.bcDict["D"]
            #value = Expression(str(value),degree=2)
            bc = DirichletBC(V, value ,boundary_markers,i.patch)
            dirichlet.append(bc)
        if "N" in i.bcDict:
            value = i.bcDict["N"]
            value = Expression(str(value),degree=2)
            neumann.append(value*v*ds(i.patch))
        if "Rec" in i.bcDict:
            #p = Expression(str(i.bcDict["Rec"]),degree=2)
            nonLinBC.append(q(u,p)*v*ds(i.patch))
    
    F = -D*(dot(grad(u), grad(v))*dx) + f*v*dx + D*(sum(neumann) + sum(nonLinBC))
    
    
    problem = NonlinearVariationalProblem(F,u, dirichlet,J=derivative(F, u))
    solver = NonlinearVariationalSolver(problem)
    solver.parameters["newton_solver"]["linear_solver"] = "cg"
    solver.parameters["newton_solver"]["preconditioner"] = "ilu"
    
    solver.solve()
    
    u.rename('u','u')
    vtkfile << u
    u_e = interpolate(u_exact_exp(p),V)
    return u,u_e
    #plot(u)
#    errorNorm.append(errornorm(u_exact(p), u))
#    
#    x = [0.2,0.2]
#    errorPoint.append(u(x)-f(x))

vtkfile = File("./sol/solution.pvd")

errorNorm = []
errorPoint = []
timing = []

for t in [32]:#np.arange(16,32,16):
    
    p = {
         "R":pow(10,2),
         "q": 10,
         "k_on": 111.6,
         "rho": 0.2,
         "D":10
         }
    
    cellList = [
            {"center":[0,0,0],"radius":p["rho"],"bcDict":{"Rec":True}} #positive values are secretion
            ]
    mesh, boundary_markers,subdomains = meshGen(Sphere(Point(0,0,0),1),t,cellList)
    subdomains[0].bcDict["D"] = u_exact_exp(p)
    start = process_time_ns()
    u,u_e = mySolve(mesh, boundary_markers,subdomains,p)
    end = process_time_ns()
    elapsed = end-start
    timing.append([t,elapsed])
    x = np.arange(p["rho"],1,0.05)
    y_exact = []
    y_exact_project = []
    y_sol = []
    
    for i in x:
        y_exact.append(u_exact([i,0,0],p))
#    for i in x:
#        y_exact_project.append(u_e([i,0,0]))
    for i in x:
        y_sol.append(u([i,0,0]))
    
    plt.plot(x,y_exact)
    #plt.plot(x,y_exact_project)
    plt.plot(x,y_sol)
    











