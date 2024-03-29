import os
import pickle as pl
import sys
from sys import argv

import dill
import dolfin as dlf
import fenics as fcs
from mpi4py import MPI

from thesis.main.ParameterSet import ParameterSet

comm = MPI.COMM_WORLD
rank = comm.rank
LINEAR_FENICS_SOLVERS = ["default", "mumps", "petsc", "umfpack"]

def prepare_solver():
    global fq

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
        fq = pl.load(f)
    with open(pickle_loc + "result_path", "rb") as f:
        result_path = pl.load(f)
    with open(pickle_loc + "patch_list", "rb") as f:
        patch_list = pl.load(f)

    p_xml = ET.parse(pickle_loc + "p_xml").getroot()
    p = ParameterSet.deserialize_from_xml(p_xml)

    mesh = dlf.Mesh()
    with dlf.XDMFFile(dlf.MPI.comm_world, mesh_file) as f:
        f.read(mesh)

    boundary_markers = fcs.MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    with fcs.HDF5File(fcs.MPI.comm_world, markers_path, "r") as f:
        f.read(boundary_markers, "/boundaries")

    V = fcs.FunctionSpace(mesh, "CG", 1)
    u = fcs.TrialFunction(V)
    v = fcs.TestFunction(V)

    ds = fcs.Measure("ds", domain=mesh, subdomain_data=boundary_markers)


    if p.get_misc_parameter("linear", "numeric").get_in_sim_unit(type=bool):
        solver, u = linear_solver(u, v, V, p, ds, integral, dirichlet, boundary_markers, fq=fq)
    else:
        try:
            const_prod = p.get_misc_parameter("constant_production", "numeric").get_in_sim_unit(type=bool)
        except AttributeError:
            const_prod = False
        if const_prod:
            solver, u = production_non_linear_solver(u, v, V, p, ds, integral, dirichlet, boundary_markers, fq=fq)
        else:
            solver, u = non_linear_solver(u, v, V, p, ds, integral, dirichlet, boundary_markers, fq=fq)

    return solver, result_path, u


def linear_solver(u, v, V, p, ds, integral, dirichlet, boundary_markers, fq=None):
    dirichlet_bc = []
    integral_bc = []

    for bc in dirichlet:
        dirichlet_bc.append(fcs.DirichletBC(V, bc[0], boundary_markers, bc[1]))

    for bc in integral:

        p_bc = bc["p"]
        f = dill.loads(bc["f"])
        area = bc["area"]
        patch = bc["patch"]

        # field_quantity = bc["field_quantity"]

        for k, value in p_bc.items():
            try:
                p_bc[k] = fcs.Constant(value)
            except:
                pass

        r = f(u, p_bc, fq, area=area)

        integral_bc.append(r * v * ds(patch))

    D = fcs.Constant(p.get_physical_parameter_by_field_quantity("D", fq).get_in_sim_unit())
    kd = fcs.Constant(p.get_physical_parameter_by_field_quantity("kd", fq).get_in_sim_unit())

    F = -D * (fcs.dot(fcs.grad(u), fcs.grad(v)) * fcs.dx) - u * kd * v * fcs.dx + D * (sum(integral_bc))

    a = fcs.lhs(F)
    L = fcs.rhs(F)
    u = fcs.Function(V)
    problem = fcs.LinearVariationalProblem(a, L, u, dirichlet_bc)

    # instantiates fenics solver
    solver = fcs.LinearVariationalSolver(problem)

    solver.parameters["linear_solver"] = p.get_misc_parameter(
        "linear_solver", "numeric").get_in_sim_unit()
    if (solver.parameters["linear_solver"] in LINEAR_FENICS_SOLVERS):
        solver.parameters["preconditioner"] = p.get_misc_parameter(
            "preconditioner", "numeric").get_in_sim_unit()

    solver.parameters["krylov_solver"]["absolute_tolerance"] = p.get_misc_parameter(
        "krylov_atol", "numeric").get_in_sim_unit(type=float)
    solver.parameters["krylov_solver"]["relative_tolerance"] = p.get_misc_parameter(
        "krylov_rtol", "numeric").get_in_sim_unit(type=float)

    return solver, u


def production_non_linear_solver(u, v, V, p, ds, integral, dirichlet, boundary_markers, fq=None):
    dirichlet_bc = []
    integral_bc = []

    v = fcs.TestFunction(V)
    u = fcs.Function(V)
    du = fcs.TrialFunction(V)

    for bc in dirichlet:
        dirichlet_bc.append(fcs.DirichletBC(V, bc[0], boundary_markers, bc[1]))

    for bc in integral:

        p_bc = bc["p"]
        f = dill.loads(bc["f"])
        area = bc["area"]
        patch = bc["patch"]

        # field_quantity = bc["field_quantity"]

        for k, value in p_bc.items():
            try:
                p_bc[k] = fcs.Constant(value)
            except:
                pass

        r = f(u, p_bc, fq, area=area)
        integral_bc.append(r * v * ds(patch))

    D = fcs.Constant(p.get_physical_parameter_by_field_quantity("D", fq).get_in_sim_unit())
    kd = fcs.Constant(p.get_physical_parameter_by_field_quantity("kd", fq).get_in_sim_unit())
    alpha = fcs.Constant(p.get_physical_parameter_by_field_quantity("alpha", fq).get_in_sim_unit())

    F = -D * (fcs.dot(fcs.grad(u), fcs.grad(v)) * fcs.dx) - u * kd * v * fcs.dx + D * (sum(integral_bc)) + v * alpha * fcs.dx

    problem = fcs.NonlinearVariationalProblem(F, u, dirichlet_bc, J=fcs.derivative(F, u, du))

    # instantiates fenics solver
    solver = fcs.NonlinearVariationalSolver(problem)

    # solver.parameters['newton_solver']['relaxation_parameter'] = 1.5
    solver.parameters["newton_solver"]["linear_solver"] = p.get_misc_parameter(
        "linear_solver", "numeric").get_in_sim_unit()
    if not (solver.parameters["newton_solver"]["linear_solver"] in LINEAR_FENICS_SOLVERS):
        solver.parameters["newton_solver"]["preconditioner"] = p.get_misc_parameter(
            "preconditioner", "numeric").get_in_sim_unit()

    solver.parameters["newton_solver"]["absolute_tolerance"] = p.get_misc_parameter(
        "newton_atol", "numeric").get_in_sim_unit(type=float)

    solver.parameters["newton_solver"]["relative_tolerance"] = p.get_misc_parameter(
        "newton_rtol", "numeric").get_in_sim_unit(type=float)

    solver.parameters["newton_solver"]["krylov_solver"]["absolute_tolerance"] = p.get_misc_parameter(
        "krylov_atol", "numeric").get_in_sim_unit(type=float)
    solver.parameters["newton_solver"]["krylov_solver"]["relative_tolerance"] = p.get_misc_parameter(
        "krylov_rtol", "numeric").get_in_sim_unit(type=float)

    return solver, u

def non_linear_solver(u, v, V, p, ds, integral, dirichlet, boundary_markers, fq=None):
    dirichlet_bc = []
    integral_bc = []

    v = fcs.TestFunction(V)
    u = fcs.Function(V)
    du = fcs.TrialFunction(V)

    for bc in dirichlet:
        dirichlet_bc.append(fcs.DirichletBC(V, bc[0], boundary_markers, bc[1]))

    for bc in integral:

        p_bc = bc["p"]
        f = dill.loads(bc["f"])
        area = bc["area"]
        patch = bc["patch"]

        # field_quantity = bc["field_quantity"]

        for k, value in p_bc.items():
            try:
                p_bc[k] = fcs.Constant(value)
            except:
                pass

        r = f(u, p_bc, fq, area=area)
        integral_bc.append(r * v * ds(patch))

    D = fcs.Constant(p.get_physical_parameter_by_field_quantity("D", fq).get_in_sim_unit())
    kd = fcs.Constant(p.get_physical_parameter_by_field_quantity("kd", fq).get_in_sim_unit())

    F = -D * (fcs.dot(fcs.grad(u), fcs.grad(v)) * fcs.dx) - u * kd * v * fcs.dx + D * (sum(integral_bc))

    problem = fcs.NonlinearVariationalProblem(F, u, dirichlet_bc, J=fcs.derivative(F, u, du))

    # instantiates fenics solver
    solver = fcs.NonlinearVariationalSolver(problem)

    # solver.parameters['newton_solver']['relaxation_parameter'] = 1.5
    solver.parameters["newton_solver"]["linear_solver"] = p.get_misc_parameter(
        "linear_solver", "numeric").get_in_sim_unit()
    if not (solver.parameters["newton_solver"]["linear_solver"] in LINEAR_FENICS_SOLVERS):
        solver.parameters["newton_solver"]["preconditioner"] = p.get_misc_parameter(
            "preconditioner", "numeric").get_in_sim_unit()

    solver.parameters["newton_solver"]["absolute_tolerance"] = p.get_misc_parameter(
        "newton_atol", "numeric").get_in_sim_unit(type=float)

    solver.parameters["newton_solver"]["relative_tolerance"] = p.get_misc_parameter(
        "newton_rtol", "numeric").get_in_sim_unit(type=float)

    solver.parameters["newton_solver"]["krylov_solver"]["absolute_tolerance"] = p.get_misc_parameter(
        "krylov_atol", "numeric").get_in_sim_unit(type=float)
    solver.parameters["newton_solver"]["krylov_solver"]["relative_tolerance"] = p.get_misc_parameter(
        "krylov_rtol", "numeric").get_in_sim_unit(type=float)

    return solver, u


import signal as sig


def alarm_handler(signum, frame):
    print("Solver timed out")
    raise TimeoutError


timeout = 60 ** 2

if __name__ == "__main__":

    solver, result_path, u = prepare_solver()

    if rank == 0:
        sig.signal(sig.SIGALRM, alarm_handler)
        sig.alarm(timeout)
        try:
            signal = sys.stdin.readline()
        except TimeoutError as e:
            signal = "TERMINATE"
        sig.alarm(0)

    else:
        signal = None

    signal = comm.bcast(signal, root=0)

    if signal == "START":

        file = result_path + "tmp_{fq}".format(fq=fq)

        try:
            solver.solve()
        except Exception as e:
            sys.stdout.write(str(e))
            sys.stdout.write("solution_failed")
            exit(0)

        if rank == 0:
            os.makedirs(result_path, exist_ok=True)

        with fcs.HDF5File(comm, file + ".h5", "w") as f:
            f.write(u, "field")

        if rank == 0:
            sys.stdout.write(file)
            exit(0)
        else:
            pass
            exit(0)

    if signal == "TERMINATE":
        if rank == 0:
            sys.stdout.write("TIMEOUT: Terminating solver thread")
            exit(0)
        else:
            pass
            exit(0)
