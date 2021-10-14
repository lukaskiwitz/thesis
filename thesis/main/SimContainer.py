#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 12:21:51 2019

@author: kiwitz
"""

import os
from typing import *

#
import lxml.etree as et
import numpy as np
import pandas as pd

from thesis.main.Entity import Entity
from thesis.main.EntityType import EntityType
from thesis.main.FieldProblem import GlobalProblem
from thesis.main.InternalSolver import InternalSolver
from thesis.main.ParameterSet import ParameterSet, MiscParameter
from thesis.main.ScanContainer import ScanSample
from thesis.main.TaskRecord import ClassRecord
from thesis.main.my_debug import message, debug

Element = et.Element


class InternalSolverNotRegisteredError(Exception): pass


class EntityTypeNotRegisteredError(Exception): pass


class SimContainer:
    """

    Container to hold all simulation object. Can be used as standalone or aggregate with
    StateManager.

    :ivar p: parmeters set, that get handed down to global problems.
    :ivar entity_list: List of all entities in the simulation
    :ivar global_problems: list of global problems in this simulation

    :ivar top_path: path to the top level directory to serve as reference for ext_cache
    :ivar path: (sub) path to simulation results
    :ivar relative_path:
    :ivar worker_sub_path:

    :ivar entity_templates: list of available entity templates
    :ivar internal_solvers: list of available internal solvers
    :ivar marker_lookup: lookup dict to manually set the mesh function value for exported properties
    :ivar markers: list of entity properties, which should be exported as mesh functions
    :ivar t: current time
    :ivar record: Class record object to store execution timing data

    :ivar default_sample: scan sample to reset sim object before each parameter scan

    :vartype p: Parameterset
    :vartype entity_list: List[Entity]
    :vartype global_problems: List[GlobalProblem]
    :vartype relative_path: str
    :vartype scan_path: str
    :vartype worker_sub_path: str
    :vartype entity_templates: List[EntityType]
    :vartype internal_solvers: List[InternalSolver]
    :vartype marker_lookup: Mapping[str,int]
    :vartype markers: List[str]
    :vartype t: float
    :vartype record: ClassRecord
    """

    def __init__(self, parameter_set: ParameterSet) -> None:

        self.p: ParameterSet = parameter_set
        self.entity_list: List[Entity] = []
        self.global_problems: List[GlobalProblem] = []

        self.top_path: str = "./"
        self.path: str = "./"
        self.relative_path: str = "./"
        self.worker_sub_path: str = ""

        self.entity_templates: List[EntityType] = []
        self.internal_solvers: List[InternalSolver] = []
        self.marker_lookup: Dict[Union[str, int, float], int] = {}
        self.t: float = 0
        self._time_log_df: pd.DataFrame = pd.DataFrame()
        self.record: ClassRecord = ClassRecord("SimContainer")

        self.markers: List[str] = ["type_name", "IL-2_surf_c"]
        self.default_sample: ScanSample = None

    def init_logs(self) -> None:

        """
        Cleans old logs and writes directory structure for log files.

        """

        if os.path.exists(self.get_current_path()):
            for i in os.listdir(self.get_current_path()):
                if i.endswith(".log"):
                    message("removing old logs {log}".format(log=i))
                    os.remove(self.get_current_path() + i)

    def initialize(self) -> None:

        """

        Setup task before simulation.
        Registers entities with global problems.
        """

        self.register_entites()

    def set_ext_cache(self, ext_cache: str) -> None:
        """
        Sets the external cache directory for all global_problems.
        TODO some abstraction to differentiate between spatial and mean field classes

        :param ext_cache: external cache directory path
        """

        for field in self.global_problems:
            field.ext_cache = ext_cache

    def register_entites(self) -> None:

        """

        Registers entities with global problems based on matching field_quantity values
        between GlobalProblem and Entity interactions.

        """
        for field in self.global_problems:
            fq = field.field_quantity
            for entity in self.entity_list:
                entity.update_bcs()
                if fq in [i.field_quantity for i in entity.interactions]:
                    field._add_entity(entity)

            if hasattr(field, "domain_template"):
                field.set_outer_domain(field.domain_template.get_domain(self.p, field.registered_entities))

    def get_entity_by_name(self, name: str) -> Entity:

        """

        Returns first entity with that matches given name


        """

        for i in self.entity_list:
            if i.name == name:
                return i

    def get_entities_by_type(self, entity_type: EntityType) -> Entity:

        """
        returns all entities of type

        :param entity_type:
        """

        raise NotImplementedError

    def run(self, T: List[float], number_of_replicats: int = 1) -> None:
        """
        Runs the simulation of given time range.
        :param T: time range; dt is calculated from the differences of successive elements.
        """

        run_task = self.record.start_child("run")

        for replicat_index in range(number_of_replicats):

            # self.path = os.path.join(original_path, "replicat_{}".format(replicat_index))
            for field in self.global_problems:
                field.initialize_run(self.p, self.top_path, self.path, self.get_tmp_path())

            self._pre_replicat(self, 0 + 1, replicat_index, T[0 + 1], T)  # internal use
            self.pre_replicat(self, 0 + 1, replicat_index, T[0 + 1], T)  # user defined

            self.t = T[0]
            for time_index, t in enumerate(T[1:]):
                run_task.info.update({"time_index": time_index})
                run_task.update_child_info()

                self._pre_step(self, time_index + 1, replicat_index, T[time_index + 1], T)  # internal use

                self.pre_step(self, time_index + 1, replicat_index, T[time_index + 1], T)  # user defined

                dt = t - T[time_index]

                # self.move_cells(time_index, dt)

                run_task.start_child("step")
                self.step(dt, time_index, replicat_index)
                run_task.stop_child("step")

                # time_index += 1
                assert t == self.t

                self._post_step(self, time_index, replicat_index, t, T)  # internal use

                self.post_step(self, time_index, replicat_index, t, T)  # user defined

            self._post_replicat(self, T[-1], replicat_index, t, T)  # internal use
            self.post_replicat(self, T[-1], replicat_index, t, T)  # user defined

        run_task.stop()

    def move_cells(self, time_index: int, dt: float) -> None:

        return None
        """
        TODO: the code below needs to be reworked, because some of the field have been removed

        """
        from time import sleep
        for field in self.global_problems:
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
                    field.solver.compile(tmp_path)
                # field.save_mesh(time_index, self.path)

    def get_tmp_path(self) -> str:
        """
        gets current temporary directory
        :return:tmp dir path
        """
        tmp_path = os.path.join(
            os.path.join(self.path, self.worker_sub_path),
            "solver_tmp")

        return tmp_path

    def get_current_path(self) -> str:

        """Returns current path"""

        return self.path

    def _pre_replicat(self, sc: 'SimContainer', time_index: int, replicat_index: int, t: float, T: List[float]) -> None:
        pass

    def _post_replicat(self, sc: 'SimContainer', time_index: int, replicat_index: int, t: float,
                       T: List[float]) -> None:
        pass

    def pre_replicat(self, sc: 'SimContainer', time_index: int, replicat_index: int, t: float, T: List[float]) -> None:
        pass

    def post_replicat(self, sc: 'SimContainer', time_index: int, replicat_index: int, t: float, T: List[float]) -> None:
        pass

    def _pre_step(self, sc: 'SimContainer', time_index: int, replicat_index: int, t: float, T: List[float]) -> None:
        pass

    def _post_step(self, sc: 'SimContainer', time_index: int, replicat_index: int, t: float, T: List[float]) -> None:
        pass

    def pre_step(self, sc: 'SimContainer', time_index: int, replicat_index: int, t: float, T: List[float]) -> None:

        """
        This method can be overwritten to call custom code before each time step. Must have this signature.

        :param sc: this SimContainer
        :param time_index: current time index
        :param t: current time
        :param T: complete time range
        """
        return None

    def post_step(self, sc: 'SimContainer', time_index: int, replicat_index: int, t: float, T: List[float]) -> None:
        """
        This method can be overwritten to call custom code after each time step. Must have this signature.

        :param sc: this SimContainer
        :param time_index: current time index
        :param t: current time
        :param T: complete time range
        """
        return None

    def apply_type_changes(self, replicat_index) -> None:

        """
        Applies pending entity type changes. Should be executed between time steps.

        """
        for i, entity in enumerate(self.entity_list):
            entity.get_surface_area()
            if not entity.change_type == "":
                debug("changing type for entity {id}".format(id=entity.id))
                entity_type = self.get_entity_type_by_name(entity.change_type)
                assert entity_type is not None

                internal_solver = self.get_internal_solver_by_name(entity_type.internal_solver)
                entity.set_cell_type(entity_type, internal_solver, replicat_index)
                entity.change_type = ""

    def reapply_entity_types(self, replicat_index):

        for i, entity in enumerate(self.entity_list):
            entity.get_surface_area()

            debug("changing type for entity {id}".format(id=entity.id))
            entity_type = self.get_entity_type_by_name(entity.type_name)
            assert entity_type is not None

            internal_solver = self.get_internal_solver_by_name(entity_type.internal_solver)
            entity.set_cell_type(entity_type, internal_solver, replicat_index)

    def step(self, dt: float, time_index: int, replicat_index: int) -> None:

        """
        advanches simulation by dt

        :param dt: delta t for this timestep
        :param time_index: time index for this timestep
        :param replicat_index: index for this replicat

        """

        self.apply_type_changes(replicat_index)

        for field in self.global_problems:
            tmp_path = self.get_tmp_path()
            # field.path = self.relative_path
            field.update_step(self.p, self.path, time_index, tmp_path)
            field.step(self.t, dt, time_index, tmp_path)

        for i, entity in enumerate(self.entity_list):
            entity.step(self.t, dt)

        self.t = self.t + dt

    def get_number_of_entities(self) -> Dict[str, int]:

        """
        Return the number of entities per entity type

        :return: dictionary with entity type names as keys
        """
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

    def add_entity(self, entity: Entity) -> Entity:

        """
        Adds entity to simulation

        """

        entity.p.update(self.p)

        if len(self.entity_list) > 0:
            entity.id = self.entity_list[-1].id + 1
        else:
            entity.id = 0

        entity.p.update(MiscParameter("id", entity.id), overwrite=True)
        self.entity_list.append(entity)
        return entity

    def add_entity_type(self, template: EntityType) -> EntityType:
        """
        Adds entity template to simulation

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
            # todo quick fix!
            template.interactions = self.entity_templates[i].interactions

            self.entity_templates[i] = template

        return template

    def add_internal_solver(self, internal_solver: InternalSolver) -> InternalSolver:

        """
        Registers and internal solver object

        :return:
        """
        if not self.get_internal_solver_by_name(internal_solver.name):
            self.internal_solvers.append(internal_solver)
        else:
            i = self.internal_solvers.index(self.get_internal_solver_by_name(internal_solver.name))
            self.internal_solvers[i] = internal_solver

        return internal_solver

    def get_internal_solver_by_name(self, name: str) -> Union[InternalSolver, None]:

        """Returns first internal solver that matches name"""

        for s in self.internal_solvers:
            if s.name == name:
                return s
        return None
        # raise InternalSolverNotRegisteredError("could not find internal solver {n}".format(n=name))

    def get_entity_type_by_name(self, name: str) -> Union[EntityType, None]:
        for e in self.entity_templates:
            if e.name == name:
                return e
        return None
        # raise EntityTypeNotRegisteredError("could not find entity type {n}".format(n=name))

    def add_problem(self, field: GlobalProblem) -> GlobalProblem:

        """

        Adds global problem to simulation

        """

        self.global_problems.append(field)
        return field

    def save_subdomains(self) -> None:

        """

        saves Mesh Function with subdomains domains colored according to field.get_sub_domain_vis

        :return:

        """

        for o, i in enumerate(self.global_problems):
            message("writing to {f}".format(f="{path}cache".format(path=self.get_current_path())))
            self.subdomain_files[o].write(i.get_sub_domains_vis())

    def save_domain(self) -> None:

        """
        Saves Mesh Function with outer domains colored according to field.get_outer_domain_vis

        """

        for o, i in enumerate(self.global_problems):
            self.domain_files[o].write(self.global_problems[0].get_outer_domain_vis("type_name"))

    # def save_markers(self, time_index: int) -> Dict[str, str]:
    #
    #     raise NotImplementedError
    #
    #     """TODO dirty solution. this should be on the FieldProblem"""
    #
    #     path_dict = {}
    #     top_dir = os.path.join(self.get_current_path(), "entity_markers")
    #     for marker_key in self.markers:
    #         marker_dir = os.path.join(top_dir, marker_key)
    #         os.makedirs(marker_dir, exist_ok=True)
    #
    #         marker, new_lookup = self.global_problems[0].get_sub_domains_vis(marker_key, lookup=self.marker_lookup)
    #         self.marker_lookup.update(new_lookup)
    #         marker_path = os.path.join(marker_dir, "marker_{n}.xdmf".format(n=time_index))
    #         with fcs.XDMFFile(fcs.MPI.comm_world, marker_path) as f:
    #             f.write(marker)
    #         path_dict[marker_key] = marker_path
    #
    #     return path_dict

    def to_xml(self) -> Element:
        raise NotImplementedError
        pass

    def from_xml(self, e: Element):
        raise NotImplementedError
        pass

    def set_parameters(self, p: ParameterSet):

        """
        Sets parameter set for sim container
        :return:
        """
        self.p = p
