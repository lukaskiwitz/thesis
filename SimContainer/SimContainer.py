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
import dolfin as dlf
import json
import os

class SimContainer:
    def __init__(self):    
        self.entityList = []
        self.fields = []
        self.path = "./"
    
    def initialize(self,**kwargs):
        self.subdomainFiles = []
        self.domainFiles = []
        self.fieldFiles = []
        self.boundaryMarkers = []
        
        
        os.makedirs(self.path,exist_ok=True)
        
        for i in self.fields:
            self.subdomainFiles.append(fcs.XDMFFile(fcs.MPI.comm_world,self.path+"cache/subdomain_"+i.fieldName+".xdmf"))
            self.domainFiles.append(fcs.XDMFFile(fcs.MPI.comm_world,self.path+"cache/domain_"+i.fieldName+".xdmf"))
            self.fieldFiles.append("field_"+i.fieldName)
            
        for field in self.fields:
            fq = field.fieldQuantity
            for entity in self.entityList:
                entity.updateBCs()
                if fq in entity.fieldQuantities:
                    field.addEntity(entity)
            if "load_subdomain" in kwargs:
                field.generateMesh(cache=True,path_prefix=self.path,**kwargs)
            else:
                field.generateMesh(cache=True,path_prefix=self.path,**kwargs)
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
#            field.computeBoundaryFlux()
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
            self.subdomainFiles[o].write(i.getSubDomainsVis(key="q_il2"))
    def saveDomain(self):
        for o,i in enumerate(self.fields):
            self.domainFiles[o].write(self.fields[0].getOuterDomainVis("R"))
    def saveFields(self,t):
        if not os.path.isdir(self.path+"sol/distplot"):
            try:
                os.mkdir(self.path+"sol/distplot")
            except:
                pass
            
        for o,i in enumerate(self.fields):
#            self.fieldFiles[o].write_checkpoint(i.getField(),"il2",t,dlf.cpp.io.XDMFFile.Encoding.HDF5)
            with fcs.HDF5File(fcs.MPI.comm_world,self.path+"sol/distplot/"+self.fieldFiles[o]+"_"+str(t)+"_distPlot.h5","w") as f:
                f.write(i.getField(),i.fieldName)
            with fcs.XDMFFile(fcs.MPI.comm_world,self.path+"sol/"+self.fieldFiles[o]+"_"+str(t)+".xdmf") as f:
                f.write(i.getField(),t)
    