#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 13:56:50 2019

@author: Lukas Kiwitz
"""
from __future__ import print_function

from typing import List

import fenics as fcs

import BC
from MySubDomain import MySubDomain
from ParameterSet import ParameterSet
from my_debug import message


class MySolver:
    pass


class MyLinearSoler(MySolver):
    """
    class to solve the stationary diffusion equation in 2D/3D with nonlinear boundary condition

    :var solver: instance of fenics linear solver
    :vartype solver: fcs.LinearVariationalSolver

    :var p: stores modell parameters: R,q,k_on,D,rho
    :vartype p: Dict

    :var u: stores solution
    :vartype u: fcs.Function

    :var dirichlet: stores dirichlet BCs; set by self.compileSolver()
    :vartype dirichlet: List[BC.DirichletBC]

    :var integralBC: stores nonlinear BCs; set by self.compileSolver()
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
        self.p: ParameterSet = ParameterSet("solver_dummy",[])
        self.mesh = None
        self.boundary_markers = None
        self.field_quantity: str = ""

        super().__init__()

    def compileSolver(self):

        self.V = fcs.FunctionSpace(self.mesh, "P", 1)

        # define trialfunction u, testfunction v from functionspace V

        u = fcs.TrialFunction(self.V)
        v = fcs.TestFunction(self.V)

        # checks for parameters
        D = fcs.Constant(self.p.get_physical_parameter_by_field_quantity("D",self.field_quantity).get_in_sim_unit())
        kd = fcs.Constant(self.p.get_physical_parameter_by_field_quantity("kd", self.field_quantity).get_in_sim_unit())

        # defines measure on boundary
        ds = fcs.Measure("ds", domain=self.mesh, subdomain_data=self.boundary_markers)

        # iterates over subdomain list to set boundary conditions

        self.dirichlet = []
        self.integralBC = []
        for i in self.subdomains:
            e = i["entity"]
            patch = i["patch"]
            bc = e.get_BC(self.field_quantity)
            if isinstance(bc, BC.OuterBC):
                pass
            if type(bc) == BC.DirichletBC or type(bc) == BC.OuterDirichletBC:
                # Dirichlet BC
                message("patch" + str(patch))
                self.dirichlet.append(bc.get_BC(self.V, self.boundary_markers, patch))
            if type(bc) == BC.Integral or type(bc) == BC.OuterIntegral:
                # p_update = {"kd":self.p["kd"],"D":self.p["D"]}
                self.integralBC.append(bc.get_BC(u,area = e.get_surface_area()) * v * ds(patch))

        F = -D * (fcs.dot(fcs.grad(u), fcs.grad(v)) * fcs.dx) - u * kd * v * fcs.dx + D * (sum(self.integralBC))

        # Defines variational problem as F == 0 with
        # with respect to u
        a = fcs.lhs(F)
        L = fcs.rhs(F)
        u = fcs.Function(self.V)
        problem = fcs.LinearVariationalProblem(a, L, u, self.dirichlet)

        # instantiates fenics solver
        self.solver = fcs.LinearVariationalSolver(problem)

        self.solver.parameters["linear_solver"] = self.p.get_misc_parameter(
            "linear_solver", "numeric").get_in_sim_unit()
        self.solver.parameters["preconditioner"] = self.p.get_misc_parameter(
            "preconditioner", "numeric").get_in_sim_unit()
        self.solver.parameters["krylov_solver"]["absolute_tolerance"] = self.p.get_misc_parameter(
            "krylov_atol","numeric").get_in_sim_unit()
        self.solver.parameters["krylov_solver"]["relative_tolerance"] = self.p.get_misc_parameter(
            "krylov_rtol","numeric").get_in_sim_unit()

        # sets field u to be trialfunction
        self.u = u

    def solve(self) -> fcs.Function:
        """
        Wrapper around fenics solve()
        """
        # calls fenics solver; renames u for proper vtk output and returns solution u
        self.solver.solve()
        self.u.rename(self.field_quantity, self.field_quantity)

        return self.u
