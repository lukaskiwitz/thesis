import os

import fenics as fcs
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import N_A, pi
from scipy.integrate import odeint

# os.environ['LD_LIBRARY_PATH'] = '/home/kiwitz/anaconda3/envs/fenics/lib:$LD_LIBRARY_PATH'

# q1 = fcs.Constant(1)
# r1 = fcs.Constant(1)

# q2 = fcs.Constant(1)
# r2 = fcs.Constant(1)



distance = 10

class Left_boundary(fcs.SubDomain):
    def inside(self, x, on_boundary):
        return fcs.near(x[0], 0) and on_boundary


class Right_boundary(fcs.SubDomain):
    def inside(self, x, on_boundary):
        return fcs.near(x[0], distance) and on_boundary


class topBottom(fcs.SubDomain):
    def inside(self, x, on_boundary):
        return (fcs.near(x[1], 1) or fcs.near(x[1], 0)) and on_boundary


mesh = fcs.RectangleMesh(fcs.Point(0,0),fcs.Point(distance,10),10,10)
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

u0 = fcs.Expression("0",degree=1)
#
u_n = fcs.Function(M)
fcs.assign(u_n,fcs.interpolate(u0,M))




# d = fcs.Constant(0.001)# receptor production






t = 0

x_list = []
u1_l = []
u2_l = []
A = (5*4*pi**2)
l = np.sqrt(10)
T = []
# assigner = fcs.FunctionAssigner(M,M)

k_on  = (111.6)/(60**2) * 1e9 * 10e15 * 1/N_A
k_iR = 0.5
k_off = 1
k_rec = 0.1
k_iC = 10
k_deg = 0
v_trans = 15000/(60**2)/A*l

D = fcs.Constant(10)#mu² per s
q = 10/A*l

r1 = 1000/A*l
r2 = 0
c1 = 0
c2 = 0

e1 = 0
e2 = 1000/A*l

def solve_odes(u1, u2, r1, c1, e1, r2, c2, e2, t, dt):
    def ode(x,t,u1,u2,T):
        r1, c1, e1, r2, c2, e2 = x
        os = 1#np.sin(0.1*T)*1
        x_n = [
            v_trans -(os*k_on * u1  + k_iR)* r1 + k_off*c1 + k_rec*e1,
            os*k_on*u1*r1 - (k_off + k_iC)*c1,
            k_iC*c1 - (k_rec + k_deg)*e1,
            v_trans -(k_on * u2 + k_iR) * r2 + k_off * c2 + k_rec * e2,
            k_on * u2 * r2 - (k_off + k_iC) * c2,
            k_iC * c2 - (k_rec + k_deg) * e2
        ]
        return x_n


    return odeint(ode, [r1, c1, e1, r2, c2, e2], [0,dt],args=(u1,u2,t))[1]

os.makedirs("/home/kiwitz/solution",exist_ok=True)
vtk = fcs.File("/home/kiwitz/solution/sol.pvd")

DT = 1
dt = fcs.Constant(DT)
timespan = 50
u = fcs.Function(M)

for n in range(int(timespan/DT)):


    u1 = fcs.assemble(u * ds(1))
    u2 = fcs.assemble(u * ds(2))
    u1_l.append(u1*(1e9*10e15/N_A))
    u2_l.append(u2*(1e9*10e15/N_A))
    ode_result = solve_odes(
        u1,
        u2,
        r1,
        c1,
        e1,
        r2,
        c2,
        e2,
        t, DT)

    r1, c1, e1, r2, c2, e2 = ode_result
    x_list.append(ode_result)


    R1 = fcs.Constant(r1)
    R2 = fcs.Constant(r2)
    C1 = fcs.Constant(c1)
    C2 = fcs.Constant(c2)
    K_on = fcs.Constant(k_on)
    Q = fcs.Constant(q)

    u = fcs.TrialFunction(M)
    F = -D * (u * v * fcs.dx + dt * fcs.dot(fcs.grad(u), fcs.grad(v)) * fcs.dx) + u_n * v * fcs.dx + dt * D * (
            (Q - K_on*R1 * u + C1*k_off) * v * ds(1)
            + (Q- K_on*R2 * u + C2*k_off) * v * ds(2)
    )

    u = fcs.Function(M)

    problem = fcs.LinearVariationalProblem(fcs.lhs(F), fcs.rhs(F), u, [])
    solver = fcs.LinearVariationalSolver(problem)

    solver.parameters["krylov_solver"]["relative_tolerance"] = 1e-10
    solver.parameters["linear_solver"] = "gmres"
    solver.solve()


    t += DT
    u.rename("u","u")
    vtk << u
    fcs.assign(u_n, u)
    # plt.figure()
    # plt.colorbar(fcs.plot(u_n))

plt.subplot(2,2,1)
plt.colorbar(fcs.plot(u_n))

plt.subplot(2,2,2)
for x in np.transpose(x_list)[0:3]:
    plt.plot(x)
plt.legend(["R1","C1","E1"],loc=1)
plt.subplot(2,2,4)
for x in np.transpose(x_list)[3:6]:
    plt.plot(x)
plt.legend(["R2","C2","E2"],loc=1)
plt.subplot(2,2,3)
plt.plot(u1_l)
plt.plot(u2_l)
plt.legend(["u1","u2"])
plt.show()




# y_list  = []
# for i in range(10000):
#     ode_result = solve_odes(
#             1,
#             1,
#             r1, c1, e1, r2, c2, e2,
#             t, DT)
#
#     r1, c1, e1, r2, c2, e2 = ode_result
#     y_list.append(ode_result)
#
# for x in np.transpose(y_list)[0:3]:
#     plt.plot(x)
# plt.show()
