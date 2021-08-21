#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 12:21:51 2019

@author: kiwitz
"""

import os
import time
from typing import *

#
import fenics as fcs
import lxml.etree as et
import pandas as pd
import numpy as np

from thesis.main.Entity import Entity
from thesis.main.EntityType import CellType, EntityType
from thesis.main.FieldProblem import FieldProblem
from thesis.main.InternalSolver import InternalSolver
from thesis.main.ParameterSet import ParameterSet, MiscParameter
from thesis.main.SimComponent import SimComponent
from thesis.main.my_debug import message, total_time, debug
from thesis.main.TaskRecord import TaskRecord, ClassRecord


class SimContainer(SimComponent):
    """

    Sim container

    :ivar p:
    :ivar entity_list:
    :ivar fields:
    :ivar path:
    :ivar scan_path:
    :ivar worker_sub_path:
    :ivar subdomain_files:
    :ivar domain_files:
    :ivar field_files:
    :ivar boundary_marker:
    :ivar entity_templates:
    :ivar marker_lookup:
    :ivar internal_solvers:
    :ivar t:
    :ivar _time_log_df:
    :ivar record:
    :ivar orig_stdout:
    :ivar markers: properties to mark

    :vartype p: Parameterset
    :vartype entity_list: List[Entity]
    :vartype fields: List[FieldProblem]
    :vartype path: str
    :vartype scan_path: str
    :vartype worker_sub_path: str
    :vartype subdomain_files: List
    :vartype domain_files: List
    :vartype field_files: List
    :vartype boundary_marker: List[fcs.MeshFunction]
    :vartype entity_templates: List[EntityType]
    :vartype marker_lookup: Mapping[str,int]
    :vartype internal_solvers: List[InternalSolver]
    :vartype t: int
    :vartype _time_log_df: pd.DataFrame
    :vartype record: ClassRecord
    :vartype orig_stdout:
    :vartype markers: List[str]
    """

    def __init__(self, parameter_set: ParameterSet) -> None:
        self.p: ParameterSet = parameter_set
        self.entity_list: List[Entity] = []
        self.fields: List[FieldProblem] = []
        self.path: str = "./"
        self.scan_path: str = self.path

        self.worker_sub_path: str = ""
        self.subdomain_files: List = []
        self.domain_files: List = []
        self.field_files: List = []
        self.boundary_markers: List = []
        self.entity_templates: List[EntityType] = []
        self.marker_lookup: Mapping[str, int] = {}
        self.internal_solvers: List[InternalSolver] = []
        self.t: float = 0
        self._time_log_df: pd.DataFrame = pd.DataFrame()
        self.record: ClassRecord = ClassRecord("SimContainer")
        self.orig_stdout = None
        self.markers: List[str] = ["type_name", "IL-2_surf_c"]

        # self.unit_length_exponent: int = 1

    def init_logs(self) -> None:

        """
        Cleans old logs and writes directory structure for log files.

        """

        if os.path.exists(self.get_current_path()):
            for i in os.listdir(self.get_current_path()):
                if i.endswith(".log"):
                    message("removing old logs {log}".format(log=i))
                    os.remove(self.get_current_path() + i)

    def init_xdmf_files(self) -> None:

        """

        Writes directory structure for XDMF and sets up file lists (subdomain_files, domain_files, field_files)
        """

        os.makedirs(self.path, exist_ok=True)
        for i in self.fields:
            self.subdomain_files.append(
                fcs.XDMFFile(fcs.MPI.comm_world, self.get_current_path() + "cache/subdomain_" + i.field_name + ".xdmf"))
            self.domain_files.append(
                fcs.XDMFFile(fcs.MPI.comm_world, self.get_current_path() + "cache/domain_" + i.field_name + ".xdmf"))
            self.field_files.append("field_" + i.field_name)

    def initialize(self, **kwargs) -> None:

        """

        loads subdomain; generates mesh; adds entities to Fields
        """


        self.init_xdmf_files()

        self.ext_boundary_markers = kwargs["load_subdomain"] if "load_subdomain" in kwargs else self.path + "cache/"

        self.register_entites()

        for field in self.fields:

            field.update_parameter_set(self.p)
            field.generate_mesh(cache=True, path=self.path, path_prefix="cache", **kwargs)
            field.update_parameter_set(self.p)
            field.update_solver(self.get_tmp_path())

    def set_ext_cache(self, ext_cache: str):
        for field in self.fields:
            field.ext_cache = ext_cache

    def register_entites(self):

        for field in self.fields:
            fq = field.field_quantity
            for entity in self.entity_list:
                entity.update_bcs()
                if fq in [i.field_quantity for i in entity.interactions]:
                    field.add_entity(entity)

            field.set_outer_domain(field.domain_template.get_domain(self.p, field.registered_entities))

    def get_entity_by_name(self, name: str) -> Entity:

        """

        returns first entity with name
        :param name:

        """

        for i in self.entity_list:
            if i.name == name:
                return i

    def get_entity_by_type(self, type: CellType) -> Entity:

        """
        TODO
        returns all entities of type
        maybe add EntityType <- CellType, ... classes

        :param type:
        """

        raise NotImplementedError()

    def run(self, T: List[float]) -> None:

        run_task = self.record.start_child("run")


        self.t = T[0]
        for time_index, t in enumerate(T[1:]):


            run_task.info.update({"time_index": time_index})
            run_task.update_child_info()

            self._pre_step(self, time_index + 1, T[time_index + 1], T)  # internal use

            self.pre_step(self, time_index + 1, T[time_index + 1], T)  # user defined

            dt = t - T[time_index]

            self.move_cells(time_index, dt)

            run_task.start_child("step")
            self.step(dt, time_index)
            run_task.stop_child("step")

            time_index += 1
            assert t == self.t

            self._post_step(self, time_index, t, T)  # internal use

            self.post_step(self, time_index, t, T)  # user defined

        run_task.stop()

    def move_cells(self, time_index: int, dt: float) -> None:

        from time import sleep
        for field in self.fields:
            tmp_path = self.get_tmp_path()

            if field.moving_mesh == True and field.remesh == True:
                for i, entity in enumerate(field.registered_entities):
                    patch = entity["patch"]
                    cell = entity["entity"]

                    mc = cell.p.get_physical_parameter("mc", "motility")
                    if mc is None:
                        mc = 0
                    else:
                        mc = mc.get_in_sim_unit()

                    cell.velocity = np.random.normal(0, np.sqrt(2 * mc * dt), 3)

                    cell.move_real(dt, field.outer_domain)

                field.generate_mesh(cache=False, path=self.get_current_path(), path_prefix="cache", file_name=
                "mesh_{field}", time_index=time_index)

            elif field.remesh:
                field.generate_mesh(cache=False, path=self.get_current_path(), path_prefix="cache",
                                    file_name="mesh_{field}", time_index=time_index)

            elif field.moving_mesh == True:

                message("Moving Cells")
                if field.ale:
                    field.load_mesh(field.mesh_path + field.file_name + ".xdmf")
                    field.ale(dt, tmp_path)
                    sleep(2)
                    field.solver.compileSolver(tmp_path)
                # field.save_mesh(time_index, self.path)

    def get_tmp_path(self) -> str:

        tmp_path = os.path.join(
            os.path.join(self.path, self.worker_sub_path),
            "solver_tmp")

        return tmp_path

    def get_current_path(self) -> str:

        if self.scan_path == "./":
            return self.path
        else:
            return self.scan_path

    def _pre_step(self, sc: 'SimContainer', time_index: int, t: float, T: List[float]) -> None:
        return None

    def _post_step(self, sc: 'SimContainer', time_index: int, t: float, T: List[float]) -> None:
        return None

    def pre_step(self, sc: 'SimContainer', time_index: int, t: float, T: List[float]) -> None:
        return None

    def post_step(self, sc: 'SimContainer', time_index: int, t: float, T: List[float]) -> None:
        return None

    def apply_type_changes(self) -> None:

        for i, entity in enumerate(self.entity_list):
            entity.get_surface_area()
            if not entity.change_type == "":
                debug("changing type for entity {id}".format(id=entity.id))
                entity_type = self.get_entity_type_by_name(entity.change_type)
                assert entity_type is not None

                internal_solver = self.get_internal_solver_by_name(entity_type.internal_solver)
                entity.set_cell_type(entity_type, internal_solver)
                entity.change_type = ""
            elif not entity.type_name == "":  # todo very bad solution
                pass
                # entity_type = self.get_entity_type_by_name(entity.type_name)
                # assert entity_type is not None
                #
                # internal_solver = self.get_internal_solver_by_name(entity_type.internal_solver)
                # entity.set_cell_type(entity_type, internal_solver)
                # entity.change_type = ""

    def step(self, dt: float, time_index: int) -> None:

        """
        advanches simulation by dt

        :param dt:

        """


        self.apply_type_changes()

        for field in self.fields:
            tmp_path = self.get_tmp_path()
            field.path = self.path
            # field.unit_length_exponent = self.unit_length_exponent
            field.update_solver(tmp_path, p=self.p)
            field.step(dt, time_index, tmp_path)

        for i, entity in enumerate(self.entity_list):
            entity.step(self.t, dt)

        self.t = self.t + dt

    def get_number_of_entites(self) -> Dict[int,str]:

        result = {}
        for e in self.entity_list:

            if hasattr(e, "type_name"):
                name = e.change_type if e.type_name == "" else e.type_name
            elif hasattr(e, "change_type"):
                name = e.change_type
            else:
                name = "none"

            if name in result.keys():
                result[name] += 1
            else:
                result[name] = 1
        return result

    def add_entity(self, entity: Entity) -> None:

        """
        adds entity to simulation

        :param entity:
        """

        entity.p.update(self.p)

        if len(self.entity_list) > 0:
            entity.id = self.entity_list[-1].id + 1
        else:
            entity.id = 0

        entity.p.update(MiscParameter("id", entity.id), override = True)
        self.entity_list.append(entity)

    def add_entity_type(self, template: EntityType) -> None:
        """
        adds entity template to simulation

        :param template: the entity type register with simcontainer
        """
        # if template not in self.entity_templates:
        #     self.entity_templates.append(template)
        # else:
        #     i = self.entity_templates.index(template)
        #     self.entity_templates[i] = template

        if not self.get_entity_type_by_name(template.name):
            self.entity_templates.append(template)
        else:
            i = self.entity_templates.index(self.get_entity_type_by_name(template.name))
            #todo quick fix!
            template.interactions = self.entity_templates[i].interactions

            self.entity_templates[i] = template

    def add_internal_solver(self, internal_solver: InternalSolver) -> None:

        if not self.get_internal_solver_by_name(internal_solver.name):
            self.internal_solvers.append(internal_solver)
        else:
            i = self.internal_solvers.index(self.get_internal_solver_by_name(internal_solver.name))
            self.internal_solvers[i] = internal_solver

    def get_internal_solver_by_name(self, name: str) -> None:

        for s in self.internal_solvers:
            if s.name == name:
                return s

        debug("could not find internal solver {n}".format(n=name))
        return None

    def get_entity_type_by_name(self, name: str) -> EntityType:
        for e in self.entity_templates:
            if e.name == name:
                return e
        debug("could not find entity type {n}".format(n = name))
        return None

    def add_problem(self, field: FieldProblem) -> None:

        """

        adds field to simulation

        :param field:

        """

        # field.update_parameter_set(self.p)
        # field.unit_length_exponent = self.unit_length_exponent
        self.fields.append(field)

    def save_subdomains(self) -> None:

        """

        saves Mesh Function with subdomains domains colored according to field.get_sub_domain_vis

        :return:

        """

        for o, i in enumerate(self.fields):
            message("writing to {f}".format(f="{path}cache".format(path=self.get_current_path())))
            self.subdomain_files[o].write(i.get_sub_domains_vis())

    def save_domain(self) -> None:

        """
        saves Mesh Function with outer domains colored according to field.get_outer_domain_vis

        :return:
        """
        for o, i in enumerate(self.fields):
            self.domain_files[o].write(self.fields[0].get_outer_domain_vis("type_name"))

    def save_markers(self, time_index: int) -> Mapping[str,str]:

        path_dict = {}
        top_dir = os.path.join(self.get_current_path(), "entity_markers")
        for marker_key in self.markers:
            marker_dir = os.path.join(top_dir, marker_key)
            os.makedirs(marker_dir, exist_ok=True)

            marker, new_lookup = self.fields[0].get_sub_domains_vis(marker_key, lookup=self.marker_lookup)
            self.marker_lookup.update(new_lookup)
            marker_path = os.path.join(marker_dir, "marker_{n}.xdmf".format(n=time_index))
            with fcs.XDMFFile(fcs.MPI.comm_world, marker_path) as f:
                f.write(marker)
            path_dict[marker_key] = marker_path

        return path_dict

    def save_fields(self, n: int) -> Mapping[str,Tuple[str,str,int]]:

        """

        saves results for each element in self.fields as xdmf file with index n

        :param n:

        """
        os.makedirs(self.get_current_path() + "sol/distplot", exist_ok=True)
        result: Dict = {}

        for o, i in enumerate(self.fields):

            # markers = i.get_sub_domains_vis(marker_lookup = self.marker_lookup)
            # with fcs.XDMFFile(fcs.MPI.comm_world, self.get_current_path() + "markers_{n}.xdmf".format(n=n)) as f:
            #     f.write(markers[0])
            #
            # with fcs.XDMFFile(fcs.MPI.comm_world, self.get_current_path() + "il2_{n}.xdmf".format(n=n)) as f:
            #     f.write(markers[1])

            distplot = "sol/distplot/" + self.field_files[o] + "_" + str(n) + "_distPlot.h5"

            sol = "sol/" + self.field_files[o] + "_" + str(n) + ".xdmf"
            u = i.get_field()
            if u is not None:
                u.rename(i.field_quantity, i.field_quantity)

                with fcs.HDF5File(fcs.MPI.comm_world, self.get_current_path() + distplot, "w") as f:
                    f.write(i.get_field(), i.field_name)
                with fcs.XDMFFile(fcs.MPI.comm_world, self.get_current_path() + sol) as f:
                    f.write(i.get_field(), n)
                result[i.field_name] = (distplot, sol, o)
            else:
                result[i.field_name] = (None, None, o)
        return result

    def to_xml(self) -> et.Element:
        raise NotImplementedError
        pass

    def from_xml(self, e: et.Element):
        raise NotImplementedError
        pass

    def set_parameters(self, set: ParameterSet):

        self.p = set
