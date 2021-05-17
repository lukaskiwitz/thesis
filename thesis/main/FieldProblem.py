#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 12:48:49 2019

@author: Lukas Kiwitz
"""

import multiprocessing as mp
import os
import time

import dolfin as dlf
import fenics as fcs
import numpy as np

import thesis.main.MeshGenerator as mshGen
from thesis.main.Entity import Entity, DomainEntity
from thesis.main.MySolver import MySolver
from thesis.main.ParameterSet import ParameterSet, ParameterCollection, PhysicalParameter
from thesis.main.PostProcess import get_concentration_conversion
from thesis.main.my_debug import message, total_time, warning


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
        self.moving_mesh = False
        self.ale = False
        self.remesh = False

    def add_entity(self, entity: Entity) -> None:
        """
        :param entity:  entity to be added

        """
        self.registered_entities.append({"entity": entity, "patch": 0})

    def get_mesh_path(self, time_index=None, local=False, full = True):

        result = ""

        dynamic_mesh = self.moving_mesh or self.remesh

        if (not time_index is None) and (dynamic_mesh):
            file = self.file_name + "_{t}" + ".xdmf"
        else:
            file = self.file_name + ".xdmf"

        if (not self.ext_cache == "") and ((local == False) or (not dynamic_mesh)):
            result = self.ext_cache
            os.makedirs(os.path.join(self.path,result),exist_ok=True)
        else:
            if local:
                result = self.path_prefix
            else:
                result = os.path.join(self.path,self.path_prefix)

        os.makedirs(os.path.join(self.path,self.path_prefix), exist_ok=True)

        if full:
            return (os.path.join(result,file)).format(field=self.field_name, t=time_index)
        else:
            return file.format(field=self.field_name, t=time_index)

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

    def generate_mesh(self, path="", path_prefix="", file_name="mesh_{field}", time_index = None, **kwargs) -> None:
        """
        generates mesh

        :param path_prefix: path to store mesh data

        """
        self.file_name = file_name
        self.path_prefix = path_prefix
        self.path = path


        mesh_gen = mshGen.MeshGenerator(outer_domain=self.outer_domain, **kwargs)
        mesh_gen.entityList = self.registered_entities
        mesh_gen.dim = 3

        if self.moving_mesh or self.ale or self.remesh:
            self.ext_cache = ""
            mesh_path = os.path.join(path, self.get_mesh_path(time_index=time_index,local=True))
        else:
            mesh_path = os.path.join(path, self.get_mesh_path(time_index=time_index))

        # self.mesh_path = mesh_path
        if not kwargs["cache"] or (self.mesh_cached == "" and self.ext_cache == ""):

            mesh, boundary_markers = mesh_gen.meshGen(self.p,
                                                      load=False,
                                                      load_subdomain=False,
                                                      path=mesh_path,
                                                     )
            self.mesh_cached = file_name.format(field=self.field_name)
        else:
            mesh, boundary_markers = mesh_gen.meshGen(self.p,
                                                      load=True,
                                                      load_subdomain=True,
                                                      path=mesh_path,
                                                    )
        self.solver.mesh = mesh
        self.solver.boundary_markers = boundary_markers
        self.solver.p = self.p

    def apply_sample(self, sample):

        self.update_bcs(sample.p)
        self.p.update(sample.p, override=True)
        self.solver.p.update(sample.p, override = True)
        self.outer_domain.apply_sample(sample.outer_domain_parameter_dict)

    def update_parameter_set(self, p):
        self.p.update(p)
        if self.solver is not None:
            self.solver.p.update(p)
        self.outer_domain.update_bcs(p = p)

    def update_bcs(self, p=None) -> None:
        """
        updates boundary conditions for all child objects

        """

        self.solver.field_quantity = self.field_quantity
        self.outer_domain.update_bcs(p=p)
        for i in self.registered_entities:
            i["entity"].update_bcs()
        subdomains = self.outer_domain.get_subdomains(field_quantity=self.field_quantity)
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

    # noinspection PyPep8Naming
    def step(self, dt: float, time_index: int, tmp_path: str) -> None:
        """
        runs timestep

        :param dT: delta t for this timestep
        """

        message("Solving Field Problem")
        try:
            self.solver.solve()
        except Exception as e:
            self.u = None
            self.solver.u = None
            return None

        message("Computing Boundary Concentrations")
        self.get_boundary_concentrations(tmp_path)
        self.compute_boundary_term()
        self.save_mesh(time_index,self.path)

        # message("Computing Boundary Gradients")
        # self.get_boundary_gradients()

    def load_mesh(self, path):

        with dlf.XDMFFile(path) as f:
            f.read(self.solver.mesh)

    def save_mesh(self, time_index: int, path: str):

        path = path + self.get_mesh_path(time_index, local=True)
        with dlf.XDMFFile(path) as f:
            f.write(self.solver.mesh)

    def ale(self, dt, tmp_path):

        mesh = self.solver.mesh

        b_mesh = fcs.BoundaryMesh(mesh, "exterior")

        b_map = b_mesh.entity_map(0)
        b_coords = b_mesh.coordinates()

        for i, entity in enumerate(self.registered_entities):

            patch = entity["patch"]
            cell = entity["entity"]
            mc = cell.p.get_physical_parameter("mc", "motility")
            message("mc: {mc}".format(mc=mc))
            cell.velocity = np.random.normal(0, np.sqrt(2 * 1 * 60 * dt), 3)
            cell.move(dt)

            # message("patch {p}".format(p = patch))
            verts = np.array([], dtype=int)
            for facet in fcs.SubsetIterator(self.solver.boundary_markers, patch):
                for vertex_index in facet.entities(0):
                    verts = np.insert(verts, 0, b_map.where_equal(vertex_index)[0])

            for i in np.unique(verts):
                x = b_coords[i]

                x[0] += cell.offset[0]
                x[1] += cell.offset[1]
                x[2] += cell.offset[2]

            # mesh.snap_boundary(cell.get_compiled_subdomain())

        dlf.ALE.move(mesh, b_mesh)
        # mesh.smooth(10)

        # with fcs.XDMFFile(tmp_path + "/b_mesh.xdmf") as f:
        #     f.write(b_mesh)

    def get_field(self) -> fcs.Function:
        return self.solver.u

    def get_sub_domains(self) -> fcs.MeshFunction:
        return self.solver.boundary_markers

    def get_sub_domains_vis(self, marker_key, lookup=None) -> fcs.MeshFunction:
        """

        save meshfunction that labels cells by type_name

        :param {type_name:label}

        """
        lookup = {} if lookup is None else lookup

        mesh = self.solver.mesh
        boundary_markers = self.solver.boundary_markers
        marker = fcs.MeshFunction("double", mesh, mesh.topology().dim() - 1)
        marker.set_all(0)

        for o in self.registered_entities:
            e = o["entity"]
            p = o["patch"]

            if hasattr(e,marker_key):
                value = getattr(e,marker_key)
            elif marker_key in e.p.get_as_dictionary():
                value = e.p.get_as_dictionary()[marker_key]
            else:
                continue

            from numbers import Number
            if isinstance(value,Number):
                value = float(value)
            elif value in lookup:
                value = lookup[value]
            else:
                lookup[value] = max(lookup.values()) + 1 if len(lookup.values()) > 0 else 1
                value = lookup[e.type_name]

            for v in boundary_markers.where_equal(p):
                marker.set_value(v, value)

        return marker, lookup

    def get_entity_surface_area(self, e):

        if not "surface_area" in e.keys():
            return None
        return e["surface_area"]

    def get_outer_domain_vis(self, parameter_name="q") -> fcs.MeshFunction:
        """

        gets MeshFunction that has outer domain surfaces labeled by entity.getState(key)

        :param key: key to pass to getState()
        """
        mesh = self.solver.mesh
        boundary_markers = fcs.MeshFunction("double", mesh, mesh.topology().dim() - 1)
        boundary_markers.set_all(0)

        e = self.outer_domain
        for o in e.get_subdomains():
            o["entity"].get_subdomain().mark(boundary_markers, e.getState(parameter_name=parameter_name,
                                                                          field_quantity=self.field_quantity))
        return self.solver.boundary_markers

    def compute_boundary_term(self):
        from thesis.main.BC import BC

        for e in self.registered_entities:
            entity = e["entity"]
            for bc in entity.interactions:

                if not isinstance(bc,BC):
                    continue

                fq = bc.field_quantity
                try:
                    u = entity.p.get_physical_parameter_by_field_quantity("surf_c", fq).get_in_sim_unit()
                except:
                    pass

                g = bc.q(u, entity.p.get_as_dictionary(in_sim=True, with_collection_name=False), fq,
                         1) * entity.p.get_physical_parameter_by_field_quantity("D", fq).get_in_sim_unit()
                ule = self.p.get_misc_parameter("unit_length_exponent", "numeric").get_in_sim_unit(type=int)

                rate_conversion = entity.p.get_physical_parameter_by_field_quantity("q", fq).to_sim

                boundary = PhysicalParameter("boundary", 0, to_sim=rate_conversion)
                boundary.set_in_sim_unit(g)

                entity.p.get_collections_by_field_quantity(fq)[0].set_parameter(boundary, override=True)


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

    try:
        values = vertex_values.take(verts)
    except TypeError as e:
        warning("Failed to extract vertex values for patch index {pi}. Check your mesh settings".format(pi = patch_index))
        raise e

    vertex_sum = values.sum() / len(verts)
    result = [patch_index, sa, vertex_sum]

    return result
