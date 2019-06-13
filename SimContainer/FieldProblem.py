#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 12:48:49 2019

@author: kiwitz
"""
import MySolver
import MeshGenerator as mshGen
import dolfin as dlf
import fenics as fcs
import os 

class FieldProblem:
    
    def __init__(self):
        self.fieldName = ""
        self.res = 32
        self.meshCached = ""
        self.registeredEntities = []
        self.outerDomain = {}
        self.solver = None
        self.fieldQuantity = ""
        self.p = {}
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
    def log(self):
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
        self.solver.solver.parameters["newton_solver"]["preconditioner"] = "amg"
    def computeBoundaryFlux(self):
        boundary_markers = self.solver.boundary_markers
        mesh = self.solver.mesh
        u = self.solver.u
        n = fcs.FacetNormal(mesh)
        
        ds = fcs.Measure("ds", domain=mesh, subdomain_data=boundary_markers)
        for i in self.registeredEntities:
            entity = i["entity"]
            patch = i["patch"]
#            print("patch: "+str(patch))
            flux_key = "flux_{f}".format(f=self.fieldName)
           
#            flux = (-fcs.dot(n,fcs.grad(u)))*ds(patch)+entity.p["q"]*ds(patch)
            flux = -fcs.dot(n,fcs.grad(u))*ds(patch)
            flux_n = fcs.assemble(flux)
            
#            print("flux: "+str(fcs.assemble(fcs.dot(n,fcs.grad(u))*ds(patch)))+"flux_n"+str(flux_n))
            
#            print("flux_n:"+str(flux_n))
            if flux_key not in entity.p:
                entity.p["flux_cytokine"] = 0
            entity.p[flux_key] = flux_n #if flux_n > 0 else 0
        
    def step(self,dt):
        return self.solver.solve()
    def getFields(self):
        """TODO"""
        return self.solver.u
        
        