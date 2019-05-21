#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 14:30:17 2019

@author: chrx
"""
import PoissonSolver as MySol
from fenics import *
from mshr import generate_mesh

#vtkfile = File("./sol/solution.pvd")

meshGen = MySol.MeshGenerator()
cellList=[{"center":[0,0],"radius":0.5,"bcDict":{"N":-1}}]
meshGen.cellList= cellList

meshGen.compileSubdomains()
#meshGen.compileSubdomains()
mesh,subdomains, boundary_markers = meshGen.meshGen(8)
solver = MySol.PoissonSolver(mesh,subdomains, boundary_markers,2)
plot(solver.solve())

