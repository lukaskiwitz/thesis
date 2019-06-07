#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 12:48:49 2019

@author: kiwitz
"""
import MySolver
import MeshGenerator as mshGen

class FieldProblem:
    registeredEntities = []
    solver = None
    p = {}
    
    def __init__(self):
        pass
    def addEntity(self,entity):
        self.registeredEntities.append(entity)
    def removeEntity(self,entity):
        self.registeredEntities.remove(entity)
    def isRegistered(self,entity):
        pass
    def makeObjectList(self):
        objectList = []
        for e in self.registeredEntities:
            objectList.append(e.getSubDomain())
        return objectList
    def setSolver(self,solver):
        self.solver = solver
    def generateMesh(self,outerBoundary,res):
        
        meshGen = mshGen.MeshGenerator(outerBoundary=outerBoundary.getSubDomain())
        cellList=self.makeObjectList(self)
        meshGen.cellList= cellList
        meshGen.dim = 3

        mesh, boundary_markers = meshGen.meshGen(res)
        
        self.solver.mesh = mesh
        self.solver.boundary_markers = boundary_markers
        self.solver.p = self.p
        
        self.solver.solver.parameters["newton_solver"]["linear_solver"] = "gmres"
        self.solver.solver.parameters["newton_solver"]["preconditioner"] = "ilu"
        
        
    def generateSubDomains(self,outerBoundary):
        subDomainList = self.registeredEntities
        subDomainList.prepend(outerBoundary)
        self.solver.subdomains = subDomainList
    def updateSolver(self):
        """
        can be used to change Simulation objects if the mesh was not changed
        """
        self.solver.compileSolver()
    def step(self,dt):
        self.solver.solve()
    def getFields(self):
        """TODO"""
        return self.solver.u
        
        