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
from decimal import Decimal
import os 
import math

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
        di = {"type":str(type(self)),
              "fieldName":self.fieldName,
              "res":self.res,
              "meshCache":self.meshCached,
              "registeredEntities":[i["entity"].id for i in self.registeredEntities],
              "outerDomain":self.outerDomain["entity"].log(),
              "solver":self.solver.log(),
              "fieldQuantity":self.fieldQuantity,
              "p":self.p
                }
        return di
    def generateMesh(self,**kwargs):
        
        
            
        meshGen = mshGen.MeshGenerator(outerDomain=self.outerDomain)
        meshGen.entityList = self.registeredEntities
        meshGen.dim = 3
        res = self.res
        if not kwargs["cache"] or self.meshCached == "":
            if not os.path.isdir("./cache"):    
                os.mkdir("./cache")
            mesh, boundary_markers = meshGen.meshGen(res,load=False,path="./cache/meshCache_{field}".format(field=self.fieldName))
            self.meshCached = "./cache/meshCache_{field}".format(field=self.fieldName)
        else:
            if "load_subdomain" in kwargs:
                mesh, boundary_markers = meshGen.meshGen(res,path=self.meshCached,load=True,load_subdomain=kwargs["load_subdomain"])
            else:
                mesh, boundary_markers = meshGen.meshGen(res,path=self.meshCached,load=True)
        
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
        self.solver.solver.parameters["newton_solver"]["absolute_tolerance"] = 1e-35
        self.solver.solver.parameters["newton_solver"]["relative_tolerance"] = 1e-14

        
#        self.solver.solver.parameters["newton_solver"]["krylov_solver"]["absolute_tolerance"] = 1e-20
#        self.solver.solver.parameters["newton_solver"]["krylov_solver"]["relative_tolerance"] = 1e-5
        
#        for k,v in self.solver.solver.parameters["newton_solver"]["krylov_solver"].items():
#            print(k+": "+"%.2E"%Decimal(v))
#        print("-----------------------------------------------")
#        self.solver.solver.parameters["newton_solver"]["relative_tolerance"] = 1e-20
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
            flux = fcs.dot(n,fcs.grad(u))*ds(patch)
            flux_n = -fcs.assemble(flux)
            
#            print("flux: "+str(fcs.assemble(fcs.dot(n,fcs.grad(u))*ds(patch)))+"flux_n"+str(flux_n))
            
#            print("flux_n:"+str(flux_n))
            if flux_key not in entity.p:
                entity.p["flux_cytokine"] = 0
            entity.p[flux_key] = flux_n #if flux_n > 0 else 0
        
    def step(self,dT):
        return self.solver.solve()
    def getField(self):
        """TODO"""
        return self.solver.u
    def getSubDomains(self):
        return self.solver.boundary_markers
    def getSubDomainsVis(self,key="q"):
        mesh = self.solver.mesh
#        mesh= fcs.BoundaryMesh(mesh,"exterior")
        boundary_markers = fcs.MeshFunction("double",mesh, mesh.topology().dim()-1)
        boundary_markers.set_all(0)
        
        for o in self.registeredEntities: 
            e = o["entity"]
            e.getCompiledSubDomain().mark(boundary_markers,e.getState(key="R_il2"))
        
        return boundary_markers
    def getOuterDomainVis(self,key="q"):
        mesh = self.solver.mesh
        boundary_markers = fcs.MeshFunction("double",mesh, mesh.topology().dim() - 1)
        boundary_markers.set_all(0)
        
        e = self.outerDomain["entity"]
        e.getSubDomain().mark(boundary_markers,e.getState(key="R_il2"))
        return self.solver.boundary_markers
        
        