#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 12:22:13 2019

@author: kiwitz
"""

import MySubDomain as SD
import pandas as pd
import fenics as fcs
import BC as bc

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
        if key in self.p:
            return self.p[key]
        else:
            return 0
        
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
    def __init__(self,p1,p2,p,bcList):
        self.p1 = p1
        self.p2 = p2
        self.p = p
#        for i,bc in enumerate(bcList):
#            bcList[i].p = p
##            print(bcList[i].p)
        self.bcList = bcList
        self.subdomainDict = self.__compileSubdomains()
        super().__init__()
    def __compileSubdomains(self):
        subdomainDict = {}
        for i,o in enumerate(self.bcList):
            if isinstance(o,bc.outerBC):
                e = compiledEntity(o.expr,o,self.p)
                e.p = self.p
                e.fieldQuantity = o.fieldQuantity
                if o.expr not in subdomainDict.keys():
                    subdomainDict[o.expr] = [e]
                else:
                    subdomainDict[o.expr].append(e)
        return subdomainDict
        
    def getSubDomains(self,**kwargs):
        subdomains = []
        for i,o in enumerate(self.subdomainDict.values()):
            if "fieldQuantity" in kwargs:
                for e in o:
                    if e.fieldQuantity == kwargs["fieldQuantity"]:
                        subdomains.append({"entity":e,"patch":i+1})
            else:
                subdomains.append({"entity":o[0],"patch":i+1})
        return subdomains
    def getSubDomainGeometry(self):
        return SD.OuterCube(self.p1,self.p2)
class compiledEntity(Entity):
    def __init__(self,expr,bc,p):
        self.bc = bc
        self.p = p
        self.expr = expr
    def getSubDomain(self):
        box = "near(x[0],{dd}) || near(x[0],-{dd}) || near(x[1],{dd}) || near(x[1],-{dd}) || near(x[2],{dd}) || near(x[2],-{dd})".format(dd=self.p["dd"])
        return fcs.CompiledSubDomain(self.expr+"&&("+box+")")
    def getBC(self,fieldQuantity):
        self.bc.p = self.p
        return self.bc