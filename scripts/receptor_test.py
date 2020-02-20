#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 14:15:09 2019

@author: kiwitz
"""

import fenics as fcs
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from my_debug import message

class Left_boundary(fcs.SubDomain):
    def inside(self, x, on_boundary):
        return fcs.near(x[0], 0) and on_boundary


class Right_boundary(fcs.SubDomain):
    def inside(self, x, on_boundary):
        return fcs.near(x[0], 1) and on_boundary


class topBottom(fcs.SubDomain):
    def inside(self, x, on_boundary):
        return (fcs.near(x[1], 1) or fcs.near(x[1], 0)) and on_boundary


mesh = fcs.UnitSquareMesh(20, 20)
# fcs.plot(mesh)

P1 = fcs.FiniteElement("P", mesh.ufl_cell(), 1)
M = fcs.FunctionSpace(mesh, P1)

v = fcs.TestFunction(M)
u = fcs.TrialFunction(M)

f = fcs.Constant(0)

boundary_markers = fcs.MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
boundary_markers.set_all(0)

left = Left_boundary()
right = Right_boundary()
tb = topBottom()

left.mark(boundary_markers, 1)
right.mark(boundary_markers, 2)
tb.mark(boundary_markers, 3)

ds = fcs.Measure("ds", domain=mesh, subdomain_data=boundary_markers)

u_0 = fcs.Expression('0', degree=1)
u_n = fcs.Function(M)
fcs.assign(u_n, fcs.interpolate(u_0, M))

q1 = fcs.Expression("q",degree=1,q=1)
r1 = fcs.Expression("r",degree=1,r=1)

q2 = fcs.Expression("q",degree=1,q=1)
r2 = fcs.Expression("r",degree=1,r=1)

D = fcs.Constant(1)
kd = fcs.Constant(0)

k = fcs.Constant(2)  # binding rate
dt = fcs.Constant(1)

F = -D * (u * v * fcs.dx +
          dt * fcs.dot(fcs.grad(u), fcs.grad(v)) * fcs.dx) -\
        dt*kd*u*v*fcs.dx +\
    u_n * v * fcs.dx +\
    D*dt*(
        q1*v*ds(1) - (k*u*r1)*v*ds(1) +
        q2*v*ds(2) - (k*u*r2)*v*ds(2)
    )
# u = fcs.Function(M)
a = fcs.lhs(F)
L = fcs.rhs(F)

bc = []

t = 0

T = []

def ode(y,t,i):

    a = 0.2
    b = 0.2

    return [
        -b*i*y[0] + a,
        0#-a*y[1] + -b*i*y[0]
    ]

receptors = []
I =  []
y0 = [[0,1],[1,0]]
n = fcs.FacetNormal(mesh)
u = fcs.Function(M)
vtkfile = fcs.File('./receptor_test.pvd')

for n in range(50):

    i = [
        fcs.assemble(u*ds(1)),
        fcs.assemble(u*ds(2))
    ]
    I.append(i)
    y1 = odeint(ode,y0[0],[0,dt.values()[0]],args=(i[0],))
    y2 = odeint(ode, y0[1], [0, dt.values()[0]], args=(i[1],))
    y0 = [y1[1],y2[1]]
    receptors.append(y0)

    r1.r = y0[0][0]
    r2.r = y0[1][0]
    # if n > 0 and n < 100:
    #     r1.r = 0.1
    # else:
    #     r1.r = 10

    problem = fcs.LinearVariationalProblem(a, L, u)
    solver = fcs.LinearVariationalSolver(problem)
    # solver.parameters["linear_solver"] = "gmres"
    # solver.parameters["preconditioner"] = "amg"
    # solver.parameters["krylov_solver"]["relative_tolerance"] = 1e-10
    # solver.parameters["newton_solver"]["absolute_tolerance"] = 1
    # solver.parameters["newton_solver"]["maximum_iterations"] = 100
    t += dt

    solver.solve()
    fcs.assign(u_n, u)
    u.rename("u","u")
    vtkfile << u
    # plt.figure()
    # plt.colorbar(fcs.plot(u))
    # plt.show()

plt.figure()
plt.colorbar(fcs.plot(u))
plt.show()

receptors = np.array(receptors)

plt.figure()
plt.subplot(2,2,1)
plt.title("cell 1")
plt.plot(receptors[:,0])

# tot = receptors[:,0,0]+receptors[:,0,1]
# plt.plot(tot)

plt.subplot(2,2,2)
plt.title("cell 2")
plt.plot(receptors[:,1])

plt.subplot(2,2,3)
plt.title("I")
plt.plot(I)
plt.show()