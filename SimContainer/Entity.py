#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 12:22:13 2019

@author: kiwitz
"""

import MySubDomain as SD

class Entity:
    def __init__(self):
        self.fieldQuantities = []
        self.internalSolvers = []
        self.p = {}
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
        
        
class Cell(Entity):
    def __init__(self,center,radius,bcList):
        self.center = center
        self.radius = radius
        self.bcList = bcList
        super().__init__()
    def getSubDomain(self):
        return SD.CellSubDomain(self.center,self.radius)
    
class DomainEntity(Entity):
    def __init__(self,**kwargs):
        super().__init__()
        
class DomainSphere(DomainEntity):
    def __init__(self,center,radius,bcList):
        self.center = center
        self.radius = radius
        self.bcList = bcList
        super().__init__()
    def getSubDomain(self):
        return SD.OuterSphere(self.center,self.radius)