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
from ParameterSet import ParameterSet, ParameterCollection, PhysicalParameter
from my_debug import message


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

        self.mesh_path = ""
        self.path_prefix = ""
        self.field_name = ""
        self.mesh_cached = ""
        self.registered_entities = []
        self.outer_domain = None
        self.solver = None
        self.field_quantity = ""
        self.ext_cache = ""
        self.p: ParameterSet = ParameterSet("field_dummy", [])
        # self.unit_length_exponent: int = 1

    def add_entity(self, entity: Entity) -> None:
        """
        :param entity:  entity to be added

        """
        self.registered_entities.append({"entity": entity, "patch": 0})

    def get_mesh_path(self):

        result = ""
        file = self.file_name.format(field = self.field_name)+".xdmf"
        if not self.ext_cache == "":
            result = self.ext_cache
        else:
            result = self.path_prefix

        # message(self.path_prefix)
        # message(self.mesh_cached)



        return result +file

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


    def generate_mesh(self, path = "", path_prefix = "",file_name = "mesh_{field}", **kwargs) -> None:
        """
        generates mesh

        :param path_prefix: path to store mesh data

        """
        mesh_gen = mshGen.MeshGenerator(outer_domain=self.outer_domain, **kwargs)
        mesh_gen.entityList = self.registered_entities
        mesh_gen.dim = 3
        mesh_path = path + self.ext_cache if not (self.ext_cache == "") else path + path_prefix
        mesh_path += "/"
        self.mesh_path = mesh_path
        self.path_prefix = path_prefix
        self.file_name = file_name

        os.makedirs(mesh_path, exist_ok=True)

        if not kwargs["cache"] or (self.mesh_cached == "" and self.ext_cache == ""):


            mesh, boundary_markers = mesh_gen.meshGen(self.p,
                                                      load=False,
                                                      path=mesh_path,
                                                      file_name = file_name.format(field=self.field_name))
            self.mesh_cached = file_name.format(field=self.field_name)
        else:
            mesh, boundary_markers = mesh_gen.meshGen(self.p,
                                                      load=True,
                                                      path=mesh_path,
                                                      file_name=file_name.format(field=self.field_name))
        self.solver.mesh = mesh
        self.solver.boundary_markers = boundary_markers
        self.solver.p = self.p

    def apply_sample(self, sample):

        self.update_bcs(sample.p)

        self.outer_domain.apply_sample(sample.outer_domain_parameter_dict)

    def update_parameter_set(self, p):
        self.p.update(p)
        self.solver.p.update(p)
        self.outer_domain.update_bcs(p)


    def update_bcs(self, p = None) -> None:
        """
        updates boundary conditions for all child objects

        """
        self.solver.field_quantity = self.field_quantity
        self.outer_domain.update_bcs(p = p)
        for i in self.registered_entities:
            i["entity"].update_bcs()
        subdomains = self.outer_domain.getSubDomains(field_quantity=self.field_quantity)
        for e in self.registered_entities:
            subdomains.append(e)
        self.solver.subdomains = subdomains

    def update_solver(self, p = None) -> None:
        """
        updates solver
        """
        self.update_bcs(p = p)
        self.solver.compileSolver()


    def get_boundary_concentrations(self) -> None:

        """

        computes average concentration over each entity and stores results in p["surf_c_{field}"]


        """

        ds = fcs.Measure("ds", domain=self.solver.mesh, subdomain_data=self.solver.boundary_markers)
        u = self.get_field()

        from PostProcess import get_concentration_conversion
        f = get_concentration_conversion(self.p.get_misc_parameter("unit_length_exponent","numeric").get_in_sim_unit())

        for i in self.registered_entities:
            value = fcs.assemble(
                u * ds(i["patch"])
            ) / (i["entity"].get_surface_area())



            physical = PhysicalParameter("surf_c", 0, to_sim=1/f)
            physical.set_in_sim_unit(value)



            i["entity"].p.update(
                ParameterSet("update_dummy", [ParameterCollection("{f}".format(f=self.field_name), [
                    physical
                ], field_quantity=self.field_quantity)])
            ,override = True)

    def get_boundary_gradients(self) -> None:

        """

        computes average gradient over each entity and stores results in p["surf_g_{field}"]


        """

        ds = fcs.Measure("ds", domain=self.solver.mesh, subdomain_data=self.solver.boundary_markers)
        u = self.get_field()
        n = fcs.FacetNormal(self.solver.mesh)
        from PostProcess import get_gradient_conversion
        f = get_gradient_conversion(self.p.get_misc_parameter("unit_length_exponent","numeric").get_in_sim_unit())

        for i in self.registered_entities:
            value = fcs.assemble(
                fcs.dot(fcs.grad(u), n) * ds(i["patch"])
            ) / (i["entity"].get_surface_area())

            physical = PhysicalParameter("surf_g", 0, to_sim=1/f)
            physical.set_in_sim_unit(value)

            i["entity"].p.update(
                ParameterSet("update_dummy", [ParameterCollection("{f}".format(f=self.field_name), [
                    physical
                ], field_quantity=self.field_quantity)]),
            override = True)

    # noinspection PyPep8Naming
    def step(self, dt: float) -> None:
        """
        runs timestep

        :param dT: delta t for this timestep
        """

        self.solver.solve()
        self.get_boundary_concentrations()
        self.get_boundary_gradients()


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
