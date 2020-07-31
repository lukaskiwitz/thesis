#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 14:15:09 2019

@author: kiwitz
"""

import fenics as fcs
import matplotlib.pyplot as plt

q1 = fcs.Constant(1)
# r1 = fcs.Constant(1)

q2 = fcs.Constant(1)
# r2 = fcs.Constant(1)

D = fcs.Constant(1)
kd = fcs.Constant(1)
class Left_boundary(fcs.SubDomain):
    def inside(self,x,on_boundary):
        return fcs.near(x[0],0) and on_boundary
class Right_boundary(fcs.SubDomain):
    def inside(self,x,on_boundary):
        return fcs.near(x[0],1) and on_boundary
class topBottom(fcs.SubDomain):
    def inside(self,x,on_boundary):
        return (fcs.near(x[1],1) or fcs.near(x[1],0)) and on_boundary
    

mesh = fcs.UnitSquareMesh(10,10)
#fcs.plot(mesh)

P1 = fcs.FiniteElement("P",mesh.ufl_cell(),1)
R_1 = fcs.FiniteElement("R",mesh.ufl_cell(),0)
R_2 = fcs.FiniteElement("R",mesh.ufl_cell(),0)

M = fcs.FunctionSpace(mesh,(P1*(R_1*R_2)))

(v,v_r1,v_r2) = fcs.TrialFunction(M)
(u,u_r1,u_r2) = fcs.TestFunction(M)

#b = fcs.Expression("sin(10*(x[0]+x[1]))",degree=1)
f = fcs.Constant(0)

boundary_markers = fcs.MeshFunction("size_t",mesh,mesh.topology().dim() - 1)
boundary_markers.set_all(0)

left = Left_boundary()
right = Right_boundary()
tb = topBottom()

left.mark(boundary_markers,1)
right.mark(boundary_markers,2)
tb.mark(boundary_markers,3)

ds = fcs.Measure("ds",domain=mesh,subdomain_data=boundary_markers)


sol = fcs.Function(M)
u,r= fcs.split(sol)
u_r1,u_r2 = fcs.split(r)


sol0 = fcs.Expression(('0','1','0'),degree=1)
sol_n = fcs.Function(M)
fcs.assign(sol_n,fcs.interpolate(sol0,M))
u_n,r_n = fcs.split(sol_n)
r1_n, r2_n = fcs.split(r_n)


c = fcs.Constant(1)#binding rate
d = fcs.Constant(0.001)# receptor production
dt = fcs.Constant(0.1)

F = -D*(u*v*fcs.dx + dt*fcs.dot(fcs.grad(u), fcs.grad(v))*fcs.dx)  +  u_n*v*fcs.dx + D*dt*(v*ds(1) - v*ds(2))#D*((q1-c*u*u_r1)*v*ds(1) + (q2-c*u*u_r2)*v*ds(2))

F_r1 = (-c*u_r1*dt)*r1_n*v_r1*fcs.dx #- u_r1*v_r1*fcs.dx - u*v_r1*ds(1)

F_r2 = 0#-dt*c*u_r2*v_r2*fcs.dx + dt*d*v_r2*fcs.dx - u_r2*v_r2*fcs.dx - u*v_r2*ds(2)


bc = []

F = F + F_r1 + F_r2
t = 0

r1_l = []
r2_l = []
T = []
assigner = fcs.FunctionAssigner(M,M)

for n in range(4):
    problem = fcs.NonlinearVariationalProblem(F,sol,J=fcs.derivative(F,sol))
    solver= fcs.NonlinearVariationalSolver(problem)

    solver.parameters["newton_solver"]["linear_solver"] = "gmres"
    solver.parameters["newton_solver"]["preconditioner"] = "amg"
    solver.parameters["newton_solver"]["relative_tolerance"] = 1e-5
    #solver.parameters["newton_solver"]["absolute_tolerance"] = 1
    #solver.parameters["newton_solver"]["maximum_iterations"] = 100
    t += dt

    solver.solve()


    # plt.figure()
    # plt.colorbar(fcs.plot(u))
    # plt.show()

    fcs.assign(sol_n,sol)

    u_r1_v = u_r1.ufl_operands[0].vector()
    r1_l.append(u_r1_v[121])
    r2_l.append(u_r1_v[122])
    #print(r1_l[-1])

plt.figure()
plt.colorbar(fcs.plot(u))
plt.show()

plt.plot(r1_l)
plt.plot(r2_l)
plt.show()
# (u,u_r1,u_r2) = sol.split()
#






