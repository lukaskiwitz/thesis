#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 20 15:28:26 2019

@author: lukas
"""


from __future__ import print_function
from fenics import *
import numpy as np
from mshr import *
from math import sin, cos, pi,sqrt
from dolfin import *

class Box(SubDomain):
    bcDict = {"D":0}
    patch = 0
    def inside(self,x,on_boundary):
        tol = 10e-10
        return on_boundary and (near(x[1],1,tol) or near(x[0],1,tol) or near(x[0],-1,tol))
class Bottom(SubDomain):
    bcDict = {"D":0}
    patch = 0
    def inside(self,x,on_boundary):
        tol = 10e-10
        return on_boundary and (near(x[1],-1,tol))
class Disk(SubDomain):
    bcDict = {"D":0}
    patch = 0
    def inside(self,x,on_boundary):
        tol = 10e-2
        return on_boundary and near(np.linalg.norm(x),1,tol)
    
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
        return on_boundary and near(np.linalg.norm(x),self.r,tol)
    
def meshGen(domain,res,cells):
    subdomains = []
    box = Disk()
    #bottom = Bottom()
    subdomains.append(box)
    #subdomains.append(bottom)
    
    for cell in cells:
        c = cell["center"]
        r = cell["radius"]
        bcDict = cell["bcDict"] if "bcDict" in cell else {}
        domain -= Sphere(Point(c),r)
        cellSubDomain = Cell()
        cellSubDomain.init(c,r,bcDict)
        
        subdomains.append(cellSubDomain)
        
    mesh = generate_mesh(domain,res)
    
    boundary_markers = MeshFunction("size_t" , mesh,mesh.topology().dim() - 1)
    boundary_markers.set_all(0)  
    
    box.mark(boundary_markers,1)
    box.patch = 1
    
    for i,o in enumerate(subdomains):
        o.mark(boundary_markers,i+1)
        o.patch = i+1
    return mesh,boundary_markers,subdomains

def cellGrid(x,y):
    cellList = []
    for i in x:
        for o in y:
            cellList.append( {"center":[i,o],"radius":0.05,"bcDict":{"Rec":1}})
    return cellList