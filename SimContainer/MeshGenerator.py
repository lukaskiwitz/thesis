#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 12:06:41 2019

@author: lukas
"""

import dolfin as dlf
import mshr
import fenics as fcs


class MeshGenerator:
    outerBoundary = None
    cellList = []
    subdomains = []
    dim = 2
    def __init__(self,**kwargs):
            if "outerBoundary" in kwargs:
                self.outerBoundary = kwargs["outerBoundary"]
    def meshGen(self,resolution,**kwargs):
        domain = self.outerBoundary.getGeometry(self.dim)
        for i in self.subdomains:
            domain -= i.getGeometry(self.dim)
        if "load" in kwargs:
            path = kwargs["load"]
            mesh = dlf.Mesh(path)
            print("loading Mesh from:"+path)
        else:
            mesh = mshr.generate_mesh(domain,resolution)
        
        
        boundary_markers = fcs.MeshFunction("size_t",mesh, mesh.topology().dim() - 1)
        boundary_markers.set_all(0)
        
        
        
        for i,o in enumerate(self.subdomains):
            o.mark(boundary_markers,i+2)
            self.subdomains[i].patch = i+2
            
        self.outerBoundary.mark(boundary_markers,1)
        self.outerBoundary.patch = 1
        self.subdomains.insert(0,self.outerBoundary)
        return mesh,boundary_markers