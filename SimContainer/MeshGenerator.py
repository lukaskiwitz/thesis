#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 12:06:41 2019

@author: lukas
"""

import dolfin as dlf
#from dolfin import *
import mshr
import fenics as fcs
import numpy as np
import random
#import h5py
import pygmsh 
import meshio
import os
import json


class MeshGenerator:
    outerDomain = None
    entityList = []
    
    dim = 2
    def __init__(self,**kwargs):
            if "outerDomain" in kwargs:
                self.outerDomain = kwargs["outerDomain"]
    def meshGen(self,resolution,**kwargs):
        geom = pygmsh.opencascade.Geometry(
          characteristic_length_min=0.05,
          characteristic_length_max=0.1
          )
        
        p1 = self.outerDomain["entity"].getSubDomain().p1
        p2 = self.outerDomain["entity"].getSubDomain().p2
        
        domain =  geom.add_box(p1,np.array(p2)-np.array(p1))
        entities = []
        path = kwargs["path"] if "path" in kwargs else ""
        for i in self.entityList[1:]:
            r = i["entity"].radius
            c = i["entity"].center
            ball = geom.add_ball(c,r,char_length=0.05)
            geom.add_physical(ball,label=i["entity"].name)
            entities.append(ball)
        if not kwargs["load"]:
            if len(entities) > 0:
                geom.boolean_difference([domain],entities)
            mesh = pygmsh.generate_mesh(geom)
            print(mesh)
            meshio.write(path+".xdmf", meshio.Mesh(points=mesh.points, cells={"tetra": mesh.cells["tetra"]}))
        else:
            print("loading Mesh from:"+path+".xdmf")
        mesh = dlf.Mesh()
        with dlf.XDMFFile(dlf.MPI.comm_world, path+".xdmf") as f:
            f.read(mesh)
        print("mesh loaded")
        
        
        if "load_subdomain" in kwargs and os.path.isfile(kwargs["load_subdomain"]):
            subPath = kwargs["load_subdomain"]
            print("loading subdomain from "+kwargs["load_subdomain"])
            boundary_markers = fcs.MeshFunction("size_t",mesh,mesh.topology().dim() - 1)
            with fcs.HDF5File(fcs.MPI.comm_world,subPath,"r") as f:
                f.read(boundary_markers,"/boundaries")
            self.outerDomain["patch"] = 1
            for i,o in enumerate(self.entityList):
                a = self.outerDomain["patch"]+1
                o["patch"] = i+a
                
        else:
            boundary_markers = fcs.MeshFunction("size_t",mesh, mesh.topology().dim() - 1)
            boundary_markers.set_all(0)
            print("boundaries marked")
            self.outerDomain["entity"] .getSubDomain().mark(boundary_markers,1)
            self.outerDomain["patch"] = 1
            print("outer domain set")
            for i,o in enumerate(self.entityList):
                a = self.outerDomain["patch"]+1
                o["entity"].getCompiledSubDomain().mark(boundary_markers,i+a)
                o["patch"] = i+a
            print("loop complete")
            if "load_subdomain" in kwargs:
                subPath = kwargs["load_subdomain"]
                with fcs.HDF5File(fcs.MPI.comm_world,subPath,"w") as f:
                    f.write(boundary_markers,"/boundaries")
        self.jsonDump()
        return mesh,boundary_markers
    def jsonDump(self):
        print("dumping json cell data")
        dumpList = []
        for i in self.entityList:
            dumpList.append({"patch":i["patch"],"center":i["entity"].center})
        dump = json.dumps(dumpList)
        with open("./cell_dump.json","w") as file:
            file.write(dump)