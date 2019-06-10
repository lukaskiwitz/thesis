#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 12:48:49 2019

@author: kiwitz
"""
import MySolver
import MeshGenerator as mshGen
import dolfin as dlf
import os 

class FieldProblem:
    fieldName = ""
    res = 32
    meshCached = ""
    registeredEntities = []
    outerDomain = {}
    solver = None
    fieldQuantity = ""
    p = {}
    
    def __init__(self):
        pass
    def addEntity(self,entity):
        self.registeredEntities.append({"entity":entity,"patch":0})
    def removeEntity(self,entity):
        pass
#        self.registeredEntities.remove(entity)
    def isRegistered(self,entity):
        pass
    def setSolver(self,solver):
        solver.fieldQuantity = self.fieldQuantity
        self.solver = solver
        
    def setOuterDomain(self,domain):
        self.outerDomain = {"entity":domain,"patch":0}
    def getBoundaryValues(self):
        for b in self.solver.boundary_markers:
            pass
    def generateMesh(self,**kwargs):
        
        
        meshGen = mshGen.MeshGenerator(outerDomain=self.outerDomain)
        meshGen.entityList = self.registeredEntities
        meshGen.dim = 3
        res = self.res
        if kwargs["cache"]:
            if self.meshCached == "":
                mesh, boundary_markers = meshGen.meshGen(res)
                self.meshCached = "./cache/meshCache_{res}_{field}.xml".format(res=res,field=self.fieldName)
                if not os.path.isdir("./cache"):    
                    os.mkdir("./cache")
                file=dlf.File(self.meshCached)
                file<< mesh
            else:
                mesh, boundary_markers = meshGen.meshGen(res,load=self.meshCached)
        else:
            mesh, boundary_markers = meshGen.meshGen(res)
        
        self.solver.mesh = mesh
        self.solver.boundary_markers = boundary_markers
        self.solver.p = self.p
        
    def updateBCs(self):
        
        self.solver.fieldQuantity = self.fieldQuantity
        self.outerDomain["entity"].updateBCs()
        
        for i in self.registeredEntities:
            i["entity"].updateBCs()
            
        subdomains = [self.outerDomain]
        for e in self.registeredEntities:
            subdomains.append(e)
        self.solver.subdomains = subdomains
        
    def updateSolver(self):
        """
        can be used to change Simulation objects if the mesh was not changed
        """
        self.updateBCs()
        
        
        self.solver.compileSolver()
        self.solver.solver.parameters["newton_solver"]["linear_solver"] = "gmres"
        self.solver.solver.parameters["newton_solver"]["preconditioner"] = "ilu"
    
    def step(self,dt):
        self.updateBCs()
        return self.solver.solve()
    def getFields(self):
        """TODO"""
        return self.solver.u
        
        