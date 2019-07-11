#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 12:22:13 2019

@author: kiwitz
"""

import MySubDomain as SD
import pandas as pd
import fenics as fcs

class Entity:
    def __init__(self):
        self.name = "default"
        self.fieldQuantities = []
        self.internalSolvers = []
        self.p = {}
        self.id = 0
    def addSolver(self,solver):
        self.internalSolvers.append(solver)
    def getBC(self,fieldQuantity):
        for i in self.bcList:
            if i.fieldQuantity == fieldQuantity:
                return i
    def updateBCs(self):
        self.fieldQuantities = []
        for i in self.bcList:
            fq = i.fieldQuantity
            self.fieldQuantities.append(fq)
            i.p = self.p
    def step(self,dT):
        for i in self.internalSolvers:
            i.step(dT,self.p)
            return i.timeToStateChange(self.p)
#            print(self.p["R"])
    def log(self):
        return {}
    
    def getState(self,key="q"):
        return self.p[key]
        
class Cell(Entity):
    def __init__(self,center,radius,bcList):
        self.center = center
        self.radius = radius
        self.bcList = bcList
        super().__init__()
    def getSubDomain(self):
        return SD.CellSubDomain(self.center,self.radius)
    def getCompiledSubDomain(self):
        return fcs.CompiledSubDomain("on_boundary && abs((sqrt(pow(x[0]-c0,2)+pow(x[1]-c1,2)+pow(x[2]-c2,2))-r) <= 10e-2)",
                              c0=self.center[0],c1=self.center[1],c2=self.center[2],r=self.radius)
    def log(self):        
        di = {"type":str(type(self)),
              "id":self.id,
              "name":self.name,
              "center":self.center,
              "radius":self.radius,
              "p":self.p
              }
        return di
    def getState(self, key="q"):
        return self.p[key]
    
class DomainEntity(Entity):
    def __init__(self,**kwargs):
        super().__init__()
    def log(self):
        return {"type":str(type(self))}
class DomainSphere(DomainEntity):
    def __init__(self,center,radius,bcList):
        self.center = center
        self.radius = radius
        self.bcList = bcList
        super().__init__()
    def getSubDomain(self):
        return SD.OuterSphere(self.center,self.radius)
class DomainCube(DomainEntity):
    def __init__(self,p1,p2,bcList):
        self.p1 = p1
        self.p2 = p2
        self.bcList = bcList
        super().__init__()
    def getSubDomain(self):
        return SD.OuterCube(self.p1,self.p2)