#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 12:21:51 2019

@author: kiwitz
"""
#import FieldProblem as fp
#import Entity
import numpy as np
import pandas as pd
import fenics as fcs

class SimContainer:
    def __init__(self):    
        self.entityList = []
        self.fields = []
    
    def initialize(self):
        self.subdomainFiles = []
        self.domainFiles = []
        self.fieldFiles = []
        
        for i in self.fields:
            self.subdomainFiles.append(fcs.File("./sol/subdomain_"+i.fieldName+".pvd"))
            self.domainFiles.append(fcs.File("./sol/domain_"+i.fieldName+".pvd"))
            self.fieldFiles.append(fcs.File("./sol/field_"+i.fieldName+".pvd"))
            
        for field in self.fields:
            fq = field.fieldQuantity
            for entity in self.entityList:
                entity.updateBCs()
                if fq in entity.fieldQuantities:
                    field.addEntity(entity)

            
            field.generateMesh(cache=True)
            field.updateSolver()
    def log(self):
        lst = []
        for f in self.fields:
            lst.append(f.log())
        for e in self.entityList:
            lst.append(e.log())
        return lst
    def getEntityByName(self,name):
        for i in self.entityList:
            if i.name == name:
                return i
    def step(self,dT_min):
        
#        print("next time step"+str(self.nextTimeStep))
        
        for i,entity in enumerate(self.entityList):
            entity.step(dT_min)
#        print("steps: "+str(nextTimeSteps))
        for field in self.fields:
            field.updateSolver()
            field.step(dT_min)
            field.computeBoundaryFlux()
    def addEntity(self,entity):
        if len(self.entityList) > 0:
            entity.id = self.entityList[-1].id+1
        else:
            entity.id = 0
        self.entityList.append(entity)
    def addField(self,field):
        self.fields.append(field)
    def saveSubdomains(self):
        for o,i in enumerate(self.fields):
            self.subdomainFiles[o] << i.getSubDomainsVis("R")
    def saveDomain(self):
        for o,i in enumerate(self.fields):
            self.domainFiles[o] << self.fields[0].getOuterDomainVis("R")
    def saveFields(self):
        for o,i in enumerate(self.fields):
            self.fieldFiles[o] << i.getField()
    