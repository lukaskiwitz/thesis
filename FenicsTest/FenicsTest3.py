#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 18 13:24:43 2019

@author: lukas
"""

from __future__ import print_function
from fenics import *
import numpy as np
from mshr import *
from math import sin, cos, pi,sqrt
from dolfin import *
import matplotlib.pyplot as plt

class Box(SubDomain):
    bcDict = {"D":0}
    patch = 0
    def inside(self,x,on_boundary):
        tol = 10e-10
        return on_boundary and (near(x[1],1,tol) or near(x[0],1,tol) or near(x[0],-1,tol))
class Bottom(SubDomain):
    bcDict = {"D":"exp(-pow(5*x[0],2))"}
    patch = 0
    def inside(self,x,on_boundary):
        tol = 10e-10
        return on_boundary and (near(x[1],-1,tol))
    
class Cell(SubDomain):
    center = [0,0]
    r = 0.1
    bcDict = {}
    patch = 0    
    def init(self,c,r, bcDict):
        self.center = c
        self.r = r
        self.bcDict = bcDict

    def inside(self,x,on_boundary):
        tol = 10e-2
        return on_boundary and near(sqrt((x[0]-self.center[0])**2+(x[1]-self.center[1])**2) - self.r,0,tol)
    
def meshGen(domain,cells):
    subdomains = []
    box = Box()
    bottom = Bottom()
    subdomains.append(box)
    subdomains.append(bottom)
    
    for cell in cells:
        c = cell["center"]
        r = cell["radius"]
        bcDict = cell["bcDict"] if "bcDict" in cell else {}
        domain -= Circle(Point(c),r)
        cellSubDomain = Cell()
        cellSubDomain.init(c,r,bcDict)
        
        subdomains.append(cellSubDomain)
        
    mesh = generate_mesh(domain,16*2)
    
    boundary_markers = MeshFunction("size_t" , mesh,mesh.topology().dim() - 1)
    boundary_markers.set_all(0)  
    
    box.mark(boundary_markers,1)
    box.patch = 1
    
    for i,o in enumerate(subdomains):
        o.mark(boundary_markers,i+1)
        o.patch = i+1
    return mesh,boundary_markers,subdomains
#START
vtkfile = File("./sol/solution.pvd")
def cellGrid(x,y):
    cellList = []
    for i in x:
        for o in y:
            cellList.append( {"center":[i,o],"radius":0.05,"bcDict":{"Rec":1}})
    return cellList

#for t in (np.arange(1,2,1)):
t = 1
print(t)

x = np.linspace(-0.5,0.5,5)
cellList = cellGrid(x,x)
mesh, boundary_markers,subdomains = meshGen(Rectangle(Point(-1,1),Point(1,-1)),cellList)
subdomains[1].bcDict["D"] = "%s*exp(-pow(5*x[0],2))"%(t)
#plot(mesh)
neumann = []
dirichlet = []
nonLinBC = []

ds = Measure("ds", domain=mesh, subdomain_data=boundary_markers)

V = FunctionSpace(mesh, 'P', 3)
u = Function(V)
v = TestFunction(V)
f = Constant(0)
D = Constant(100)

def q(u):
    return 10*u/(10+u)

for i in subdomains:
    if "D" in i.bcDict:
        D = i.bcDict["D"]
        value = Expression(str(D),degree=2)
        bc = DirichletBC(V, value ,boundary_markers,i.patch)
        dirichlet.append(bc)
    if "N" in i.bcDict:
        N = i.bcDict["N"]
        value = Expression(str(N),degree=2)
        neumann.append(value*v*ds(i.patch))
    if "Rec" in i.bcDict:
        nonLinBC.append(q(u)*v*ds(i.patch))
        
F = dot(grad(u), grad(v))*dx - f*v*dx + sum(neumann)- sum(nonLinBC)
solve(F == 0, u, dirichlet)
u.rename('u','u')
vtkfile << u
plot(u)
# 