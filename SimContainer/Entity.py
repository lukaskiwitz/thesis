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
from copy import copy,deepcopy


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
        p_out = self.p
        for i in p_out.keys():
            if "numpy.float64" == str(type(p_out[i])):
                p_out[i] = float(p_out[i])
            if "ndarray" in str(type(p_out[i])):
                p_out[i] = list(p_out[i])
        di = {"type":str(type(self)),
              "id":self.id,
              "name":self.name,
              "center":self.center,
              "radius":self.radius,
              "p":self.p
              }
        return di
    
class DomainEntity(Entity):
    def __init__(self,**kwargs):
        super().__init__()
    def log(self):
        return {"type":str(type(self))}
class DomainSphere(DomainEntity):
    def __init__(self,center,radius,p,bcList):
        self.center = center
        self.radius = radius
        self.p = p
        self.bcList = bcList
        self.subdomainDict = self.__compileSubdomains()
        super().__init__()
    def __compileSubdomains(self):
        subdomainDict = {}
        for i,o in enumerate(self.bcList):
            if isinstance(o,bc.outerBC):
                e = compiledSphere(o.expr,o,self.p)
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
    def getSubDomain(self):
        return SD.OuterSphere(self.center,self.radius)
class DomainCube(DomainEntity):
    def __init__(self,p1,p2,p,bcList):
        self.p1 = p1
        self.p2 = p2
        self.p = p

        self.bcList = bcList
        self.subdomainDict = self.__compileSubdomains()
        super().__init__()
    def __compileSubdomains(self):
        subdomainDict = {}
        for i,o in enumerate(self.bcList):
            if isinstance(o,bc.outerBC):
                e = compiledCube(o.expr,o,self.p)
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
    
class compiledCube(Entity):
    def __init__(self,expr,bc,p):
        self.bc = bc
        self.p = p
        self.expr = expr
    def getSubDomain(self):
        box = "near(x[0],{d_x}) || near(x[0],-{d_x}) || near(x[1],{d_y}) || near(x[1],-{d_y}) || near(x[2],{d_z}) || near(x[2],-{d_z})".format(d_x=self.p["d_x"],d_y=self.p["d_y"],d_z=self.p["d_z"])
        print("expr "+self.expr)
        return fcs.CompiledSubDomain(self.expr+"&&("+box+") && on_boundary")
    def getBC(self,fieldQuantity):
        self.bc.p = self.p
        return self.bc
    
class compiledSphere(Entity):
    def __init__(self,expr,bc,p):
        self.bc = bc
        self.p = p
        self.expr = expr
    def getSubDomain(self):
        box = "abs((sqrt(pow(x[0]-{c0},2)+pow(x[1]-{c1},2)+pow(x[2]-{c2},2))-{r}))<= 10e-2".format(c0=0,c1=0,c2=0,r=self.p["radius"])
        return fcs.CompiledSubDomain(self.expr+"&&("+box+") && on_boundary")
    def getBC(self,fieldQuantity):
        self.bc.p = self.p
        return self.bc