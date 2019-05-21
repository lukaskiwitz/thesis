#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 13:56:50 2019

@author: chrx
"""
from __future__ import print_function
from fenics import *
from fenics import SubDomain,MeshFunction
import numpy as np
from mshr import *
from mshr import Circle,Sphere, Rectangle, Box
import mshr as mshr
from math import sin, cos, pi,sqrt
from dolfin import *
from dolfin import FunctionSpace,Point

class DummyError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
        
class MySubDomain(SubDomain):
    def __init__(self):
        self.patch = 0
    def inside(self,x,on_boundary):
        return on_boundary
class OuterBoundary(MySubDomain):
    pass
class OuterSphere(OuterBoundary):
    bcDict = {"D":0}
    center = [0,0]
    radius = 0.5
    def __init__(self,c,r):
        self.center = c
        self.radius = r
class OuterCube(OuterBoundary):
    bcDict = {"D":0}
    def __init__(self,p1,p2):
        self.p1 = p1
        self.p2 = p2
    def getGeometry(self,dim):
        if dim==2:
            return Rectangle(Point(self.p1[0],self.p1[1]),Point(self.p2[0],self.p2[1]))
        elif dim==3:
            return Box(Point(self.p1[0],self.p1[1],self.p1[2]),Point(self.p2[0],self.p2[1],self.p2[2]))
    def inside(self,x,on_boundary):
        pass
class Cell(MySubDomain):
    bcDict = {"D":1}
    def __init__(self,c,r):
        self.radius = r
        self.center = c
    def getGeometry(self,dim):
        if dim == 2:
            return Circle(Point(self.center[0],self.center[1]),self.radius)
        elif dim == 3:
            return Sphere(Point(self.center[0],self.center[1],self.center[2]),self.radius)
    
class MeshGenerator:
    outerBoundary = OuterCube([-1,-1],[1,1])
    cellList = []
    subdomains = []
    dim = 2
    def __init__(self):
            pass
    def meshGen(self,resolution):
        domain = self.outerBoundary.getGeometry(self.dim)
        for i in self.subdomains:
            domain -= i.getGeometry(self.dim)
        mesh = generate_mesh(domain,resolution)
        
        #boundary_markers = MeshFunction("size_t" , mesh,mesh.topology().dim() - 1)
        #boundary_markers.set_all(0)
        #self.outerBoundary.mark(boundary_markers,1)
        
        domain.set_subdomain(1,domain)
        for i,o in enumerate(self.subdomains):
            domain.set_subdomain(i+2,o.getGeometry(self.dim))
            o.patch = i+2
        boundary_markers = MeshFunction("size_t" , mesh,self.dim-1,mesh.domains())
        return mesh,self.subdomains, boundary_markers

    def compileSubdomains(self):
        for i in self.cellList:
            r = i["radius"]
            c = i["center"]
            if "bcDict" in i:
                bcDict = i["bcDict"]
            else:
                pass
                #raise DummyError("non bcDict!")
            subdomains.append(Cell(c,r))
        self.subdomains = subdomains
        
class PoissonSolver:
    neumann = []
    dirichlet = []
    nonLinBC = []
    subdomains = []
    dim = 2
    def __init__(self,mesh,subdomains,boundary_markers,dim):
        self.mesh = mesh
        self.dim = dim
        self.subdomains = subdomains
        self.boundary_markers = boundary_markers
        self.V = FunctionSpace(self.mesh,"P",1)
    def solve(self):
        
        u = Function(self.V)
        v = TestFunction(self.V)
        f = Constant(0)
        D = Constant(1)
        ds = Measure("ds", domain=self.mesh, subdomain_data=self.boundary_markers)
        
        
        for i in self.subdomains:
            if "D" in i.bcDict:
                value = i.bcDict["D"]
                value = Expression(str(value),degree=2)
                bc = DirichletBC(self.V, value ,self.boundary_markers,i.patch)
                self.dirichlet.append(bc)
            if "N" in i.bcDict:
                value = i.bcDict["N"]
                value = Expression(str(value),degree=2)
                self.neumann.append(value*v*ds(i.patch))
            if "Rec" in i.bcDict:
                #p = Expression(str(i.bcDict["Rec"]),degree=2)
                pass#self.nonLinBC.append(q(u,p)*v*ds(i.patch))
        
        F = -D*(dot(grad(u), grad(v))*dx) + f*v*dx + D*(sum(self.neumann) + sum(self.nonLinBC))
    
        problem = NonlinearVariationalProblem(F,u, self.dirichlet,J=derivative(F, u))
        solver = NonlinearVariationalSolver(problem)
        #solver.parameters["newton_solver"]["linear_solver"] = "cg"
        #solver.parameters["newton_solver"]["preconditioner"] = "ilu"
        
        solver.solve()
        u.rename('u','u')
        return u