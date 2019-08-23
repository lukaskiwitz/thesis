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
import Entity

class DomainTypeError(Exception):
    def __init__(self,text):
        self.text = text
    def __str__(self):
        return "Uknown Domain Type {t}".format(t=self.text)
class MeshGenerator:
    outerDomain = None
    entityList = []
    
    dim = 2
    def __init__(self,**kwargs):
            if "outerDomain" in kwargs:
                self.outerDomain = kwargs["outerDomain"]
    def meshGen(self,resolution,**kwargs):
#        print(kwargs)
        geom = pygmsh.opencascade.Geometry(
          characteristic_length_min=kwargs["min_char_length"] if "min_char_length" in kwargs else 0.001,
          characteristic_length_max=kwargs["max_char_length"]if "min_char_length" in kwargs else 0.05,
          )

        if isinstance(self.outerDomain,Entity.DomainCube):
            p1 = self.outerDomain.p1
            p2 = self.outerDomain.p2
            domain =  geom.add_box(p1,np.array(p2)-np.array(p1))
        elif isinstance(self.outerDomain,Entity.DomainSphere):
            c = self.outerDomain.center
            r = self.outerDomain.radius
            domain = geom.add_ball(c,r)
        else :
            raise DomainTypeError(type(self.outerDomain))
        entities = []
        path = kwargs["path"] if "path" in kwargs else ""
        for i in self.entityList:
            r = i["entity"].radius
            c = i["entity"].center
            ball = geom.add_ball(c,r)
            geom.add_physical(ball,label=i["entity"].name)
            entities.append(ball)
        if not kwargs["load"] or not os.path.isfile(path+".xdmf"):
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
            for i,o in enumerate(self.entityList):
                a = self.outerDomain.getSubDomains()[-1]["patch"]+1
                o["patch"] = i+a
                
        else:
            boundary_markers = fcs.MeshFunction("size_t",mesh, mesh.topology().dim() - 1)
            boundary_markers.set_all(0)
            print("boundaries marked")
            
            for i in self.outerDomain.getSubDomains():
                i["entity"].getSubDomain().mark(boundary_markers,i["patch"])
            print("outer domain set")
            for i,o in enumerate(self.entityList):
                a = self.outerDomain.getSubDomains()[-1]["patch"]+1
                o["entity"].getCompiledSubDomain().mark(boundary_markers,i+a)
                o["patch"] = i+a
            print("loop complete")
            if "load_subdomain" in kwargs:
                subPath = kwargs["load_subdomain"]
                with fcs.HDF5File(fcs.MPI.comm_world,subPath,"w") as f:
                    f.write(boundary_markers,"/boundaries")
        return mesh,boundary_markers
    def jsonDump(self,path):
        print("dumping json cell data")
        dumpList = []
        for i in self.entityList:
            dumpList.append({"patch":i["patch"],"center":i["entity"].center})
        dump = json.dumps(dumpList)
        with open(path+"cell_dump.json","w") as file:
            file.write(dump)