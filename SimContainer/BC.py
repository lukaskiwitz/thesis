#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  8 14:17:07 2019

@author: lukas
"""
import fenics as fcs
class BC:
    fieldQuantity = ""
    pass
    def __init__(self,**kwargs):
        if "fieldQuantity" in kwargs:
            self.fieldQuantity = kwargs["fieldQuantity"]
class Integral(BC):
    value = None
    p = {}
    q = lambda u,p : 0
    def __init__(self,q,**kwargs):
        self.q = q
        super().__init__(**kwargs)
    def getBC(self,u):
      return self.q(u,self.p)
        
class DirichletBC(BC):
    degree = 1
    value = None
    def __init__(self,**kwargs):
        super.__init__(**kwargs)
        
    def getBC(self,V,boundary_markers,patch):
        value = fcs.Expression(str(self.value),degree=self.degree)
        bc = fcs.DirichletBC(V, value ,self.boundary_markers,patch)
        return bc