#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  8 14:17:07 2019

@author: lukas
"""
import fenics as fcs
class BC:
    
    def __init__(self,**kwargs):
        self.fieldQuantity = ""
        if "fieldQuantity" in kwargs:
            self.fieldQuantity = kwargs["fieldQuantity"]
class Integral(BC):
    def __init__(self,q,**kwargs):
        self.value = None
        self.p = {}
        self.q = q
        super().__init__(**kwargs)
    def getBC(self,u):
      return self.q(u,self.p)
        
class DirichletBC(BC):
    
    def __init__(self,value,**kwargs):
        self.degree = 1
        self.value = value
        super().__init__(**kwargs)
        
    def getBC(self,V,boundary_markers,patch):
        value = fcs.Expression(str(self.value),degree=self.degree)
        bc = fcs.DirichletBC(V, value ,boundary_markers,patch)
        return bc
    