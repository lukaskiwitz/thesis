#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 12:21:51 2019

@author: kiwitz
"""
#import FieldProblem as fp
#import Entity

class SimContainer:
    entityList = []
    fields = []
    
    def initialize(self):
        for field in self.fields:

            fq = field.fieldQuantity
            for entity in self.entityList:
                entity.updateBCs()
                if fq in entity.fieldQuantities:
                    field.addEntity(entity)

            
            field.generateMesh(cache=True)
            field.updateSolver()
                        
    def step(self,dT):
        for entity in self.entityList:
            entity.step(dT)
        for field in self.fields:
            field.step(dT)
    def addEntity(self,entity):
        self.entityList.append(entity)
    def addField(self,field):
        self.fields.append(field)
    