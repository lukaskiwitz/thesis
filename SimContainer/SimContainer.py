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
        self.T = 0
    def initLogs(self):
        if os.path.exists(self.path):
            for i in os.listdir(self.path):
                if i.endswith(".log"):
                    print("removing old logs {log}".format(log=i))
                    os.remove(self.path+i)
    def initXdmfFiles(self):
        self.subdomainFiles = []
        self.domainFiles = []
        self.fieldFiles = []
        
        os.makedirs(self.path,exist_ok=True)
        for i in self.fields:
            self.subdomainFiles.append(fcs.XDMFFile(fcs.MPI.comm_world,self.path+"cache/subdomain_"+i.fieldName+".xdmf"))
            self.domainFiles.append(fcs.XDMFFile(fcs.MPI.comm_world,self.path+"cache/domain_"+i.fieldName+".xdmf"))
            self.fieldFiles.append("field_"+i.fieldName)
    def initialize(self,**kwargs):
        
        self.boundaryMarkers = []
        self.initXdmfFiles()
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
    def timeStepLog(self,dT):
        logFile = ""
        newLog = {"timestep":dT,"time":self.T}
#        newLog = [1]
        os.makedirs(self.path,exist_ok=True)
        for i in os.listdir(self.path):
            if i.endswith(".log"):
                logFile = self.path+i
                break
        if logFile == "":
            with open(self.path+"timeStep.log","w") as f:
                f.write(json.dumps(newLog)+",")
        else:
            with open(logFile,"a") as f:
                f.write(json.dumps(newLog)+",")
            
    def getEntityByName(self,name):
        for i in self.entityList:
            if i.name == name:
                return i
    def step(self,dT):
        
#        print("next time step"+str(self.nextTimeStep))
        
        self.T = self.T + dT
        
        
        for i,entity in enumerate(self.entityList):
            entity.step(dT)
            
#        print("steps: "+str(nextTimeSteps))
        for field in self.fields:
            field.updateSolver()
            field.step(dT)
#            field.computeBoundaryFlux()
        self.timeStepLog(dT)
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
            print("writing to {f}".format(f="{path}cache".format(path=self.path)))
            self.subdomainFiles[o].write(i.getSubDomainsVis(key="type"))
    def saveDomain(self):
        for o,i in enumerate(self.fields):
            self.domainFiles[o].write(self.fields[0].getOuterDomainVis("R"))
    def saveFields(self,t):
        os.makedirs(self.path+"sol/distplot",exist_ok=True)
        result = {}
        for o,i in enumerate(self.fields):
            distplot = self.path+"sol/distplot/"+self.fieldFiles[o]+"_"+str(t)+"_distPlot.h5"
            sol = self.path+"sol/"+self.fieldFiles[o]+"_"+str(t)+".xdmf"
            with fcs.HDF5File(fcs.MPI.comm_world,distplot,"w") as f:
                f.write(i.getField(),i.fieldName)
            with fcs.XDMFFile(fcs.MPI.comm_world,sol) as f:
                f.write(i.getField(),t)
            result[i.fieldName] = (distplot,sol)
        return result
        
    