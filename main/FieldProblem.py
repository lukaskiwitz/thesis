#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 12:48:49 2019

@author: Lukas Kiwitz
"""

import multiprocessing as mp
import os
import time

import fenics as fcs
import numpy as np

import MeshGenerator as mshGen
from Entity import Entity, DomainEntity
from MySolver import MySolver
from ParameterSet import ParameterSet, ParameterCollection, PhysicalParameter
from PostProcess import get_concentration_conversion
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
                                                      load_subdomain=False,
                                                      path=mesh_path,
                                                      file_name=file_name.format(field=self.field_name))
            self.mesh_cached = file_name.format(field=self.field_name)
        else:
            mesh, boundary_markers = mesh_gen.meshGen(self.p,
                                                      load=True,
                                                      load_subdomain=True,
                                                      path=mesh_path,
                                                      file_name=file_name.format(field=self.field_name))
        self.solver.mesh = mesh
        self.solver.boundary_markers = boundary_markers
        self.solver.p = self.p

    def apply_sample(self, sample):

        self.update_bcs(sample.p)
        self.p.update(sample.p)
        self.solver.p.update(sample.p)
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

    def update_solver(self, tmp_path, p=None) -> None:
        """
        updates solver
        """
        self.update_bcs(p=p)
        self.solver.compileSolver(tmp_path)

    def get_boundary_concentrations(self, tmp_path: str) -> None:

        def init(_mesh, _markers, _solution):

            global mesh
            mesh = _mesh
            global boundary_markers
            boundary_markers = _markers
            global u
            u = _solution
            global ds
            ds = fcs.Measure("ds", domain=mesh, subdomain_data=boundary_markers)

            global vertex_values
            vertex_values = u.compute_vertex_values(_mesh)

        ds = fcs.Measure("ds", domain=self.solver.mesh, subdomain_data=self.solver.boundary_markers)
        entity_list = [[e["patch"], self.get_entity_surface_area(e)] for e in self.registered_entities]
        pn = os.cpu_count()

        cs = int(len(entity_list) / pn)

        cells_per_worker = self.p.get_misc_parameter("cells_per_worker", "numeric").get_in_sim_unit(type=int)
        max_pool_size = self.p.get_misc_parameter("max_pool_size", "numeric").get_in_sim_unit(type=int)

        pool_size = int(len(entity_list) / cells_per_worker)
        pool_size = pool_size if pool_size >= 1 else 1
        pool_size = pool_size if pool_size < max_pool_size else max_pool_size

        chunksize = int(len(entity_list) / pool_size)

        message("Extracting boundary concentrations. Poolsize: {pool_size}, cell_per_worker: {cells_per_worker}".format(
            pool_size=pool_size,
            cells_per_worker=chunksize
        ))
        start = time.time()
        with mp.Pool(processes=pn, initializer=init,
                     initargs=(self.solver.mesh, self.solver.boundary_markers, self.solver.u)) as pool:
            result = pool.map(target, entity_list, chunksize=chunksize)
        f = get_concentration_conversion(
            self.p.get_misc_parameter("unit_length_exponent", "numeric").get_in_sim_unit(type=int))
        end = time.time()
        from my_debug import total_time
        total_time(end - start, pre="Cell concentration extracted in ")

        for i, (patch, sa, value) in enumerate(result):

            re = self.registered_entities[i]
            assert re["patch"] == patch
            if self.get_entity_surface_area(re) == None:
                re["surface_area"] = sa

            physical = PhysicalParameter("surf_c", 0, to_sim=1 / f)
            physical.set_in_sim_unit(value)

            self.registered_entities[i]["entity"].p.update(
                ParameterSet("update_dummy", [ParameterCollection("{f}".format(f=self.field_name), [
                    physical
                ], field_quantity=self.field_quantity)])
                , override=True)

        message("done pool map with chunksize {}".format(cs))

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
    def step(self, dt: float, tmp_path: str) -> None:
        """
        runs timestep

        :param dT: delta t for this timestep
        """
        from my_debug import message
        message("Solving Field Problem")
        self.solver.solve()
        message("Computing Boundary Concentrations")
        self.get_boundary_concentrations(tmp_path)
        # message("Computing Boundary Gradients")
        # self.get_boundary_gradients()

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

    def get_entity_surface_area(self, e):

        if not "surface_area" in e.keys():
            return None
        return e["surface_area"]

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


def target(entity):
    patch_index = entity[0]
    sa = entity[1]
    if sa == None:
        sa = fcs.assemble(1 * ds(patch_index))

    log = fcs.get_log_level()
    fcs.set_log_level(50)

    facets = np.array([i for i in fcs.SubsetIterator(boundary_markers, patch_index)])

    fcs.set_log_level(log)

    verts = np.unique(np.array(list(map(lambda x: x.entities(0), facets))))

    vertex_sum = vertex_values.take(verts).sum() / len(verts)

    result = [patch_index, sa, vertex_sum]

    return result
