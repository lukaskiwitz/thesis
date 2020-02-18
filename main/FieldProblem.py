#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 12:48:49 2019

@author: Lukas Kiwitz
"""

import os

import fenics as fcs

import MeshGenerator as mshGen
from Entity import Entity, DomainEntity
from MySolver import MySolver
import numpy as np


class FieldProblem:
    """

    :var field_name: string to identify this field
    :vartype field_name: str

    :var res: desired mesh resolution (currently not used)
    :vartype res: int

    :var mesh_cache: path for mesh storage (optional)
    :vartype mesh_cache: str

    :var registered_entities: entities that interact with this field
    :vartype registered_entities: List[Entity]

    :var outer_domain: outside domain
    :vartype outer_domain: DomainEntity

    :var solver: solver to be used for this field
    :vartype solver: MySolver

    :var field_quantity: physical properties associated with this field
    :vartype field_quantity: str

    :var ext_cache: path to folder containing pre generated mesh
    :vartype ext_cache: str

    :var p: parameter dict to be handed down
    :vartype p: Dict


    """
    def __init__(self) -> None:

        self.field_name = ""
        self.res = 32
        self.mesh_cached = ""
        self.registered_entities = []
        self.outer_domain = None
        self.solver = None
        self.field_quantity = ""
        self.ext_cache = ""
        self.p = {}

    def add_entity(self, entity: Entity) -> None:
        """
        :param entity:  entity to be added

        """
        self.registered_entities.append({"entity": entity, "patch": 0})

    def remove_entity(self, entity) -> None:
        """
        TODO
        """
        pass

    #        self.registered_entities.remove(entity)
    def is_registered(self, entity) -> None:
        """TODO"""
        pass

    def set_solver(self, solver: MySolver) -> None:
        """

        sets solver
        :param solver:

        """

        solver.field_quantity = self.field_quantity
        self.solver = solver

    def set_outer_domain(self, domain: DomainEntity) -> None:
        """sets outer domain
        :param domain: domain object
        """
        self.outer_domain = domain

    # def log(self):
    #     di = {"type": str(type(self)),
    #           "field_name": self.field_name,
    #           "res": self.res,
    #           "meshCache": self.mesh_cached,
    #           "registered_entities": [i["entity"].id for i in self.registered_entities],
    #           "outer_domain": self.outer_domain.log(),
    #           "solver": self.solver.log(),
    #           "field_quantity": self.field_quantity,
    #           "p": self.p
    #           }
    #     return di

    def generate_mesh(self, path_prefix: str, **kwargs) -> None:
        """
        generates mesh

        :param path_prefix: path to store mesh data

        """
        mesh_gen = mshGen.MeshGenerator(outer_domain=self.outer_domain,**kwargs)
        mesh_gen.entityList = self.registered_entities
        mesh_gen.dim = 3
        res = self.res
        mesh_path = self.ext_cache if not (self.ext_cache == "") else path_prefix
        if not (self.ext_cache == ""):
            os.makedirs(self.ext_cache, exist_ok=True)
        if not kwargs["cache"] or (self.mesh_cached == "" and self.ext_cache == ""):

            os.makedirs(path_prefix, exist_ok=True)
            mesh, boundary_markers = mesh_gen.meshGen(res, load=False,
                                                      path=mesh_path + "meshCache_{field}".format(field=self.field_name),
                                                      **kwargs)
            self.mesh_cached = "meshCache_{field}".format(field=self.field_name)
        else:
            mesh, boundary_markers = mesh_gen.meshGen(res, path=mesh_path, load=True, **kwargs)
        self.solver.mesh = mesh
        self.solver.boundary_markers = boundary_markers
        self.solver.p = self.p

    def update_bcs(self) -> None:
        """
        updates boundary conditions for all child objects

        """
        self.solver.field_quantity = self.field_quantity
        self.outer_domain.update_bcs()
        for i in self.registered_entities:
            i["entity"].update_bcs()
        subdomains = self.outer_domain.getSubDomains(field_quantity=self.field_quantity)
        for e in self.registered_entities:
            subdomains.append(e)
        self.solver.subdomains = subdomains

    def update_solver(self) -> None:
        """
        updates solver
        """
        self.update_bcs()
        self.solver.compileSolver()

        self.solver.solver.parameters["linear_solver"] = "gmres"
        self.solver.solver.parameters["preconditioner"] = "hypre_amg"
        self.solver.solver.parameters["krylov_solver"]["absolute_tolerance"] = 1e-35
        self.solver.solver.parameters["krylov_solver"]["relative_tolerance"] = 1e-5


    def get_boundary_concentrations(self) -> None:

        """

        computes average concentration over each entity and stores results in p["surf_c_{field}"]


        """

        ds = fcs.Measure("ds", domain=self.solver.mesh, subdomain_data=self.solver.boundary_markers)
        u = self.get_field()

        for i in self.registered_entities:
            i["entity"].p["surf_c_{f}".format(f=self.field_name)] = fcs.assemble(
                u * ds(i["patch"])
            )/(4 * np.pi * i["entity"].p["rho"] ** 2)

    def get_boundary_gradients(self) -> None:

        """

        computes average gradient over each entity and stores results in p["surf_g_{field}"]


        """

        ds = fcs.Measure("ds", domain=self.solver.mesh, subdomain_data=self.solver.boundary_markers)
        u = self.get_field()
        n = fcs.FacetNormal(self.solver.mesh)
        for i in self.registered_entities:
            i["entity"].p["surf_g_{f}".format(f=self.field_name)] = fcs.assemble(
                fcs.dot(fcs.grad(u), n) * ds(i["patch"])
            )/(4 * np.pi * i["entity"].p["rho"] ** 2)

    # noinspection PyPep8Naming
    def step(self, dT: float) -> None:
        """
        runs timestep

        :param dT: delta t for this timestep
        """
        return self.solver.solve()

    def get_field(self) -> fcs.Function:
        return self.solver.u

    def get_sub_domains(self) -> fcs.MeshFunction:
        return self.solver.boundary_markers

    def get_sub_domains_vis(self, key="q") -> fcs.MeshFunction:
        """

        gets MeshFunction that has subdomains surfaces labeled by entity.getState(key)
        (gets very slow)
        :param key: key to pass to getState()
        """

        mesh = self.solver.mesh
        boundary_markers = fcs.MeshFunction("double", mesh, mesh.topology().dim() - 1)
        boundary_markers.set_all(0)

        for o in self.registered_entities:
            e = o["entity"]
            #            if e.getState(key=key) > 0:
            #            print(e.getState(key=key))
            st = e.getState(key=key)
            if st == 0:
                e.getCompiledSubDomain().mark(boundary_markers, -1)
            else:
                e.getCompiledSubDomain().mark(boundary_markers, st)

        return boundary_markers

    def get_outer_domain_vis(self, key="q") -> fcs.MeshFunction:
        """

        gets MeshFunction that has outer domain surfaces labeled by entity.getState(key)

        :param key: key to pass to getState()
        """
        mesh = self.solver.mesh
        boundary_markers = fcs.MeshFunction("double", mesh, mesh.topology().dim() - 1)
        boundary_markers.set_all(0)

        e = self.outer_domain
        for o in e.getSubDomains():
            o["entity"].getSubDomain().mark(boundary_markers, e.getState(key=key))
        return self.solver.boundary_markers
