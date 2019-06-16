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

class SimContainer:
    def __init__(self):    
        self.entityList = []
        self.fields = []
    
    def initialize(self):
        for field in self.fields:

            fq = field.fieldQuantity
            for entity in self.entityList:
                entity.updateBCs()
                if fq in entity.fieldQuantities:
                    field.addEntity(entity)

            
            field.generateMesh(cache=True)
            field.updateSolver()
    def log(self):
        dataDict = {}
        for i,entity in enumerate(self.entityList):
            dataDict["entity_"+str(i)] = entity.log()
        return dataDict
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
#            field.log()

    def addEntity(self,entity):
        self.entityList.append(entity)
    def addField(self,field):
        self.fields.append(field)
    