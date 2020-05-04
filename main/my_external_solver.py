import os
import pickle as pl
import sys
from sys import argv

import dolfin as dlf
import fenics as fcs
from mpi4py import MPI

from ParameterSet import ParameterSet

comm = MPI.COMM_WORLD
rank = comm.rank

# import multiprocessing as mp
pickle_loc = argv[1]
import lxml.etree as ET

with open(pickle_loc + "dirichlet", "rb") as f:
    dirichlet = pl.load(f)
with open(pickle_loc + "integral", "rb") as f:
    integral = pl.load(f)
with open(pickle_loc + "mesh_file", "rb") as f:
    mesh_file = pl.load(f)
with open(pickle_loc + "markers_path", "rb") as f:
    markers_path = pl.load(f)
with open(pickle_loc + "field_quantity", "rb") as f:
    field_quantity = pl.load(f)
with open(pickle_loc + "result_path", "rb") as f:
    result_path = pl.load(f)
with open(pickle_loc + "patch_list", "rb") as f:
    patch_list = pl.load(f)

# import lxml.etree as ET
# p_xml = ET.XML(p_xml)
p = ParameterSet("dummy", [])
p_xml = ET.parse(pickle_loc + "p_xml").getroot()

p.deserialize_from_xml(p_xml)

mesh = dlf.Mesh()
with dlf.XDMFFile(dlf.MPI.comm_world, mesh_file) as f:
    f.read(mesh)

boundary_markers = fcs.MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
with fcs.HDF5File(fcs.MPI.comm_world, markers_path, "r") as f:
    f.read(boundary_markers, "/boundaries")

V = fcs.FunctionSpace(mesh, "P", 1)
u = fcs.TrialFunction(V)
v = fcs.TestFunction(V)

ds = fcs.Measure("ds", domain=mesh, subdomain_data=boundary_markers)

dirichlet_bc = []
integral_bc = []

for bc in dirichlet:
    dirichlet_bc.append(fcs.DirichletBC(V, bc[0], boundary_markers, bc[1]))
for bc in integral:
    linear = fcs.Constant(bc[0])
    billinear = fcs.Constant(bc[1])

    integral_bc.append((linear + billinear * u) * v * ds(bc[2]))

D = fcs.Constant(p.get_physical_parameter_by_field_quantity("D", field_quantity).get_in_sim_unit())
kd = fcs.Constant(p.get_physical_parameter_by_field_quantity("kd", field_quantity).get_in_sim_unit())

F = -D * (fcs.dot(fcs.grad(u), fcs.grad(v)) * fcs.dx) - u * kd * v * fcs.dx + D * (sum(integral_bc))

a = fcs.lhs(F)
L = fcs.rhs(F)
u = fcs.Function(V)
problem = fcs.LinearVariationalProblem(a, L, u, dirichlet_bc)

# instantiates fenics solver
solver = fcs.LinearVariationalSolver(problem)

solver.parameters["linear_solver"] = p.get_misc_parameter(
    "linear_solver", "numeric").get_in_sim_unit().replace("\"", "")

solver.parameters["preconditioner"] = p.get_misc_parameter(
    "preconditioner", "numeric").get_in_sim_unit().replace("\"", "")

solver.parameters["krylov_solver"]["absolute_tolerance"] = float(p.get_misc_parameter(
    "krylov_atol", "numeric").get_in_sim_unit())
solver.parameters["krylov_solver"]["relative_tolerance"] = float(p.get_misc_parameter(
    "krylov_rtol", "numeric").get_in_sim_unit())

if rank == 0:
    signal = sys.stdin.readline()

else:
    signal = None

signal = comm.bcast(signal, root=0)

if signal == "START":

    solver.solve()

    os.makedirs(result_path, exist_ok=True)

    file = result_path + "tmp"

    with fcs.HDF5File(comm, file + ".h5", "w") as f:
        f.write(u, "field")
    # comm.Barrier()

    if not rank == 0:
        exit(0)
    else:
        sys.stdout.write(file)
        exit(0)
