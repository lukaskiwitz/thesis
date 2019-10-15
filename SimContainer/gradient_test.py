#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 15:30:40 2019

@author: kiwitz
"""

import fenics as fcs
import dolfin as dlf
import numpy as np
import matplotlib.pyplot as plt
r = 50
mesh = fcs.UnitSquareMesh(r,r)
y = []
x = []
V = fcs.FunctionSpace(mesh,"P",1)
expr = fcs.Expression("5*exp(-10*x[0]) + exp(-S*(sqrt(pow(x[0]-0.1,2)+pow(x[1]-0.5,2)))) + exp(-S*(sqrt(pow(x[0]-0.75,2)+pow(x[1]-0.5,2))))",degree=2,S=1)
#for s in np.arange(1,100,1):
#    
#    expr.S = s
#    u = fcs.project(expr,V)
#    #fcs.plot(u)
#    V_vec = fcs.VectorFunctionSpace(mesh,"P",1)
#    grad = fcs.project(-fcs.grad(u),V_vec,solver_type="gmres")
#    norm = fcs.project(fcs.sqrt(fcs.dot(grad,grad)),V,solver_type="gmres")
#    y.append(fcs.assemble(fcs.sqrt(fcs.dot(grad,grad))*fcs.dX))
#    x.append(s)
##    plt.figure()
##    fcs.plot(norm)
#plt.figure()
#plt.plot(x,y)

expr.S = 10

u = fcs.project(expr,V)

fcs.plot(u)
V_vec = fcs.VectorFunctionSpace(mesh,"P",1)
grad = fcs.project(-fcs.grad(u),V_vec,solver_type="gmres")
div = fcs.project(fcs.div(grad),V,solver_type="gmres")
norm = fcs.project(fcs.sqrt(fcs.dot(grad,grad)),V,solver_type="gmres")
y.append(fcs.assemble(fcs.sqrt(fcs.dot(grad,grad))*fcs.dX))
#x.append(1)
plt.figure()
fcs.plot(grad)
plt.figure()
fcs.plot(div)
plt.figure()
fcs.plot(norm)