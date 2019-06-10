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
    outerDomain = None
    entityList = []
    
    dim = 2
    def __init__(self,**kwargs):
            if "outerDomain" in kwargs:
                self.outerDomain = kwargs["outerDomain"]
    def meshGen(self,resolution,**kwargs):
        domain = self.outerDomain["entity"].getSubDomain().getGeometry(self.dim)
        for i in self.entityList:
            domain -= i["entity"].getSubDomain().getGeometry(self.dim)
        if "load" in kwargs:
            path = kwargs["load"]
            mesh = dlf.Mesh(path)
            print("loading Mesh from:"+path)
        else:
            mesh = mshr.generate_mesh(domain,resolution)
        
        
        boundary_markers = fcs.MeshFunction("size_t",mesh, mesh.topology().dim() - 1)
        boundary_markers.set_all(0)
        
        
        
        self.outerDomain["entity"] .getSubDomain().mark(boundary_markers,1)
        self.outerDomain["patch"] = 1
        
        for i,o in enumerate(self.entityList):
            a = self.outerDomain["patch"]+1
            o["entity"].getSubDomain().mark(boundary_markers,i+a)
            o["patch"] = i+a
            
    
        return mesh,boundary_markers