#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 12:21:51 2019

@author: kiwitz
"""
#import FieldProblem as fp
#import Entity
import numpy as np
class SimContainer:
    def __init__(self):    
        self.entityList = []
        self.fields = []
        self.nextTimeStep = 1
    
    def initialize(self):
        for field in self.fields:

            fq = field.fieldQuantity
            for entity in self.entityList:
                entity.updateBCs()
                if fq in entity.fieldQuantities:
                    field.addEntity(entity)

            
            field.generateMesh(cache=True)
            field.updateSolver()
                        
    def step(self,dT_min):
        
        nextTimeSteps = []
        print("next time step"+str(self.nextTimeStep))
        for entity in self.entityList:
            nextStep = entity.step(self.nextTimeStep)
            if nextStep:
                nextTimeSteps.append(nextStep)
        print("steps: "+str(nextTimeSteps))
        if len(nextTimeSteps) > 0:
            self.nextTimeStep = np.min(nextTimeSteps)
        else:
            self.nextTimeStep = dT_min
        
        for field in self.fields:
            field.updateSolver()
            field.step(self.nextTimeStep)
            field.computeBoundaryFlux()

    def addEntity(self,entity):
        self.entityList.append(entity)
    def addField(self,field):
        self.fields.append(field)
    