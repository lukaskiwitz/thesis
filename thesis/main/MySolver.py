#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 13:56:50 2019

@author: Lukas Kiwitz
"""
from __future__ import print_function

import getpass
import os
import pickle as pl
import signal
import subprocess as sp
from abc import *
from typing import List, Any

import dill
import fenics as fcs
import lxml.etree as ET
import mpi4py.MPI as MPI

import thesis.main.BC as BC
from thesis.main.MySubDomain import MySubDomain
from thesis.main.ParameterSet import ParameterSet
from thesis.main.my_debug import message, warning

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


class MySolver(ABC):

    def __init__(self):
        self.p: ParameterSet = ParameterSet("solver_dummy", [])
        self.field_quantity: str = ""

    @abstractmethod
    def solve(self):pass

    @abstractmethod
    def compile(self, tmp_path: str):pass

    def kill(self):pass

    def get_solution(self) -> Any: pass


class MyDiffusionSolver(MySolver):
    """
    class to solve the stationary diffusion equation in 3D with nonlinear boundary condition

    :var solver: instance of fenics linear solver
    :vartype solver: fcs.LinearVariationalSolver

    :var p: stores modell parameters: R,q,k_on,D,rho
    :vartype p: Dict

    :var u: stores solution
    :vartype u: fcs.Function

    :var dirichlet: stores dirichlet BCs; set by self.compile()
    :vartype dirichlet: List[BC.DirichletBC]

    :var integralBC: stores nonlinear BCs; set by self.compile()
    :vartype integralBC: List[BC.Integral]

    :var dim: mesh dimension
    :vartype dim: int

    :var mesh: domain mesh; intended to be genrated by myMeshing.MeshGenerator
    :vartype mesh: fenics.Mesh

    :var V: Function Space
    :vartype V: fcs.FunctionSpace

    """

    def __init__(self):

        self.dirichlet: List[BC.DirichletBC] = []
        self.integralBC: List[BC.Integral] = []
        self.subdomains: List[MySubDomain] = []
        self.p: ParameterSet = ParameterSet("solver_dummy", [])
        self.field_quantity: str = ""
        self.mesh = None
        self.boundary_markers = None
        self.timeout = None

        super().__init__()

    def compile(self, tmp_path: str):

        def sig_handler(signum, frame):
            if self.process is not None:
                self.process.kill()
                exit(0)

        signal.signal(signal.SIGINT, sig_handler)
        signal.signal(signal.SIGTERM, sig_handler)

        self.V = fcs.FunctionSpace(self.mesh, "P", 1)
        self.u = fcs.Function(self.V)
        self.tmp_path = tmp_path + "{fq}/".format(fq=self.field_quantity)

        self.dirichlet = []
        self.integralBC = []

        dirichlet_bc_test = []
        integral_bc_test = []
        patch_list = []

        for i in self.subdomains:
            e = i["entity"]
            patch = i["patch"]
            bc = e.get_interaction(self.field_quantity)

            patch_list.append([i["patch"], i["entity"].get_surface_area()])

            if isinstance(bc, BC.DirichletBC) or isinstance(bc, BC.OuterDirichletBC):
                dirichlet_bc_test.append([bc.value, patch])

            if isinstance(bc, BC.Integral) or isinstance(bc, BC.OuterIntegral):
                # linear, billinear = bc.get_BC(1, area=e.get_surface_area())
                # integral_bc_test.append([linear, billinear, patch])

                d = {
                    "f": dill.dumps(bc.q),
                    "area": e.get_surface_area(),
                    "patch": patch,
                    "field_quantity": self.field_quantity,
                    "p": bc.p.get_as_dictionary(in_sim=True, with_collection_name=False,
                                                field_quantity=self.field_quantity)
                }
                integral_bc_test.append(d)

        markers_path = self.tmp_path + "boundary_markers.h5"
        mesh_file = self.tmp_path + "mesh.xdmf"

        with fcs.HDF5File(fcs.MPI.comm_world, markers_path, "w") as f:
            f.write(self.boundary_markers, "/boundaries")

        with fcs.XDMFFile(mesh_file) as f:
            f.write(self.mesh)

        p_xml = ET.ElementTree(self.p.serialize_to_xml())

        pickle_loc = self.tmp_path + "pickle/"
        os.makedirs(pickle_loc, exist_ok=True)

        with open(pickle_loc + "patch_list", "wb") as f:
            pl.dump(patch_list, f)
        with open(pickle_loc + "dirichlet", "wb") as f:
            pl.dump(dirichlet_bc_test, f)
        with open(pickle_loc + "integral", "wb") as f:
            pl.dump(integral_bc_test, f)
        with open(pickle_loc + "mesh_file", "wb") as f:
            pl.dump(mesh_file, f)
        with open(pickle_loc + "markers_path", "wb") as f:
            pl.dump(markers_path, f)
        p_xml.write(pickle_loc + "p_xml")
        with open(pickle_loc + "field_quantity", "wb") as f:
            pl.dump(self.field_quantity, f)
        with open(pickle_loc + "result_path", "wb") as f:
            pl.dump(self.tmp_path, f)

        if hasattr(self, "process"):
            if not self.process == None:
                self.process.kill()

        max_nodes = self.p.get_misc_parameter("max_mpi_nodes", "numeric").get_in_sim_unit(type=int)
        dofs_per_node = self.p.get_misc_parameter("dofs_per_node", "numeric").get_in_sim_unit(type=int)

        user = getpass.getuser()
        dof_n = self.V.dim()

        mpi_nodes = int(dof_n / dofs_per_node)

        if mpi_nodes == 0:
            mpi_nodes = 1

        mpi_nodes = mpi_nodes if mpi_nodes <= max_nodes else max_nodes

        dofs_per_node = int(dof_n / mpi_nodes)

        message("Launching {n} mpi process(es) to solve system of {dn} dofs with {dofs_p_n} per node".format(
            n=mpi_nodes, dn=dof_n, dofs_p_n=dofs_per_node))

        from thesis.main import __path__ as main_path

        ext_solver_path = main_path._path[0] + "/my_external_solver.py"

        if hasattr(self, "process"):
            self.process.kill()

        cpu_list = "0-0{e}".format(e=mpi_nodes - 1) if mpi_nodes < 11 else "0:{e}".format(e=mpi_nodes - 1)

        p = sp.Popen(
            ["nice", "-n", "19", "mpiexec", "-n", str(mpi_nodes), "python", ext_solver_path,
             pickle_loc], stdin=sp.PIPE, stdout=sp.PIPE)

        self.process = p

    def kill(self):
        if hasattr(self, "process"):
            self.process.kill()

    def solve(self) -> fcs.Function:

        from sys import exit
        def sig_handler(signum, frame):
            if self.process is not None:
                self.process.kill()
                exit(0)

        signal.signal(signal.SIGINT, sig_handler)
        signal.signal(signal.SIGTERM, sig_handler)

        class SolutionFailedError(Exception):
            pass

        limit = 1

        for o in range(limit + 1):

            if o == limit:
                raise SolutionFailedError

            try:
                signal_out, signal_err = self.process.communicate(b"START", timeout=self.timeout)
                for i in signal_out.decode("utf-8").split("\n"):
                    message(i)

                file = str(signal_out).split("\\n")[-1].replace("'", "")
                if file == "solution_failed":
                    raise SolutionFailedError
                    return -1
                if os.path.exists(file + ".h5"):
                    message("loading solution from {f}".format(f=file))
                    with fcs.HDF5File(comm, file + ".h5", "r") as f:
                        f.read(self.u, "field")
                    self.u.rename(self.field_quantity, self.field_quantity)
                    return self.u
                else:
                    message("Something went wrong. Solution file {f}.h5 doesn't exist".format(f=file))
                    message(o)
                    continue


            except (FileNotFoundError, sp.TimeoutExpired, ValueError) as e:
                self.process.kill()
                warning("External solver timed out or crashed, restarting worker threads")
                self.compile(self.tmp_path)

    def get_solution(self) -> Any:

        return self.u

class MyMeanFieldSolver(MySolver):

    def solve(self):
        pass

    def compile(self, tmp_path: str):
        pass

    def get_solution(self) -> Any:

        return 1e-17
