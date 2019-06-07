#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 12:06:41 2019

@author: lukas
"""
import MyEntities as entity
import dolfin as dlf
import mshr
import fenics as fcs


class MeshGenerator:
    outerBoundary = entity.OuterCube([-1,-1],[1,1])
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
        return mesh,self.subdomains, boundary_markers

    def compileSubdomains(self):
        self.subdomains = []
        for i in self.cellList:
            r = i["radius"]
            c = i["center"]
            if "bcDict" in i:
                bcDict = i["bcDict"]
            else:
                pass
                #raise DummyError("non bcDict!")
            newCell = entity.Cell(c,r)
            newCell.bcDict = bcDict
            self.subdomains.append(newCell)