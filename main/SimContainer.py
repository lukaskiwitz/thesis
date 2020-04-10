#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 12:21:51 2019

@author: kiwitz
"""

import os
import time
from typing import Dict

#
import fenics as fcs
import lxml.etree as et

from Entity import Entity
from EntityType import CellType, EntityType
from FieldProblem import FieldProblem
from InternalSolver import InternalSolver
from ParameterSet import ParameterSet
from SimComponent import SimComponent
from my_debug import message, total_time, debug


class SimContainer(SimComponent):
    """

    Sim container


    :var entity_list: stores entities in simulation
    :vartype entity_list: List[Entity]

    :var fields: field in simulation
    :vartype fields: List[FieldProblem]

    :var path: path to store simulation results
    :vartype path: str

    :var subdomain_files:
    :vartype subdomain_files: List[fcs.XDMFFile]

    :var domain_files:
    :vartype domain_files: List[fcs.XDMFFile]

    :var field_files:
    :vartype field_files: List[str]

    :var boundary_markers:
    :vartype boundary_markers: List[fcs.MeshFunction]

    :var ext_boundary_markers:
    :vartype ext_boundary_markers: fcs.MeshFunction

    :var T: current time
    :vartype T: float

    :var entity_templates: list of entity templates
    :vartyp entity_templates: List[EntityTemplate]

    """

    def __init__(self, parameter_set: ParameterSet) -> None:
        self.p = parameter_set
        self.entity_list = []
        self.fields = []
        self.path = "./"
        self.T = 0
        self.subdomain_files = []
        self.domain_files = []
        self.field_files = []
        self.boundary_markers = []
        self.entity_templates = []
        self.internal_solvers = []
        #self.unit_length_exponent: int = 1

    def init_logs(self) -> None:

        """
        Cleans old logs and writes directory structure for log files.

        """

        if os.path.exists(self.path):
            for i in os.listdir(self.path):
                if i.endswith(".log"):
                    message("removing old logs {log}".format(log=i))
                    os.remove(self.path + i)

    def init_xdmf_files(self) -> None:

        """

        Writes directory structure for XDMF and sets up file lists (subdomain_files, domain_files, field_files)
        """

        os.makedirs(self.path, exist_ok=True)
        for i in self.fields:
            self.subdomain_files.append(
                fcs.XDMFFile(fcs.MPI.comm_world, self.path + "cache/subdomain_" + i.field_name + ".xdmf"))
            self.domain_files.append(
                fcs.XDMFFile(fcs.MPI.comm_world, self.path + "cache/domain_" + i.field_name + ".xdmf"))
            self.field_files.append("field_" + i.field_name)

    def initialize(self, **kwargs) -> None:

        """

        loads subdomain; generates mesh; adds entities to Fields
        """

        self.init_xdmf_files()

        self.ext_boundary_markers = kwargs["load_subdomain"] if "load_subdomain" in kwargs else self.path + "cache/"

        for field in self.fields:
            field.update_parameter_set(self.p)

            fq = field.field_quantity
            for entity in self.entity_list:
                entity.update_bcs()
                if fq in entity.fieldQuantities:
                    field.add_entity(entity)
            field.generate_mesh(cache=True, path= self.path, path_prefix="cache/", **kwargs)
            field.update_solver()

    def get_entity_by_name(self, name: str) -> Entity:

        """

        returns first entity with name
        :param name:

        """

        for i in self.entity_list:
            if i.name == name:
                return i

    def get_entity_by_type(self,type: CellType) -> Entity:

        """
        TODO
        returns all entities of type
        maybe add EntityType <- CellType, ... classes

        :param type:
        """

        raise NotImplementedError()

    def run(self, T, dt):

        for time_index, t in enumerate(T):
            self._pre_step(self, time_index, t, T)  # internal use
            self.pre_step(self, time_index, t, T)  # user defined

            start = time.process_time()
            self.step(dt)
            end = time.process_time()
            total_time(end - start, pre = "Total time ", post=" for step number {n}".format(n=time_index))

            self._post_step(self,time_index,t,T)# internal use
            self.post_step(self, time_index, t, T)# user defined

    def _pre_step(self,sc,time_index,t,T):
        return None

    def _post_step(self, sc, time_index,t,T):
        return None

    def pre_step(self,sc,time_index,t,T):
        return None

    def post_step(self, sc, time_index,t,T):
        return None

    def step(self, dt: float) -> None:

        """
        advanches simulation by dt

        :param dt:

        """

        for i, entity in enumerate(self.entity_list):
            entity.get_surface_area()
            if not entity.change_type == "":
                debug("changing type for entity {id}".format(id=entity.id))
                entity_type = self.get_entity_type_by_name(entity.change_type)
                internal_solver = self.get_internal_solver_by_name(entity_type.internal_solver)
                entity.set_cell_type(entity_type, internal_solver)
                entity.change_type = ""

        for field in self.fields:
            # field.unit_length_exponent = self.unit_length_exponent
            field.update_solver(p = self.p)
            field.step(dt)

        for i, entity in enumerate(self.entity_list):
            entity.step(self.T, dt)



        self.T = self.T + dt

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
        self.entity_list.append(entity)

    def add_entity_type(self,template: EntityType):
        """
        adds entity template to simulation

        :param template
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
            self.entity_templates[i] = template

    def add_internal_solver(self, internal_solver: InternalSolver):

        if not self.get_internal_solver_by_name(internal_solver.name):
            self.internal_solvers.append(internal_solver)
        else:
            i = self.internal_solvers.index(self.get_internal_solver_by_name(internal_solver.name))
            self.internal_solvers[i] = internal_solver

    def get_internal_solver_by_name(self,name: str):


        for s in self.internal_solvers:
            if s.name == name:
                return  s

        message("could not find internal solver {n}".format(n = name))
        return None

    def get_entity_type_by_name(self,name: str)->EntityType:
        for e in self.entity_templates:
            if e.name == name:
                return e
        message("could not find entity type {n}".format(n = name))
        return None

    def add_field(self, field: FieldProblem) -> None:

        """

        adds field to simulation

        :param field:

        """

        field.update_parameter_set(self.p)
        # field.unit_length_exponent = self.unit_length_exponent
        self.fields.append(field)

    def save_subdomains(self) -> None:

        """

        saves Mesh Function with subdomains domains colored according to field.get_sub_domain_vis

        :return:

        """

        for o, i in enumerate(self.fields):
            message("writing to {f}".format(f="{path}cache".format(path=self.path)))
            self.subdomain_files[o].write(i.get_sub_domains_vis(key="type_int"))

    def save_domain(self) -> None:

        """
        saves Mesh Function with outer domains colored according to field.get_outer_domain_vis

        :return:
        """
        for o, i in enumerate(self.fields):
            self.domain_files[o].write(self.fields[0].get_outer_domain_vis("R"))

    def save_fields(self, n: int) -> Dict:

        """

        saves results for each element in self.fields as xdmf file with index n

        :param n:

        """
        os.makedirs(self.path + "sol/distplot", exist_ok=True)
        result: Dict = {}
        for o, i in enumerate(self.fields):
            distplot = "sol/distplot/" + self.field_files[o] + "_" + str(n) + "_distPlot.h5"

            sol = "sol/" + self.field_files[o] + "_" + str(n) + ".xdmf"

            with fcs.HDF5File(fcs.MPI.comm_world, self.path+distplot, "w") as f:
                f.write(i.get_field(), i.field_name)
            with fcs.XDMFFile(fcs.MPI.comm_world, self.path+sol) as f:
                f.write(i.get_field(), n)
            result[i.field_name] = (distplot, sol, o)
        return result

    def to_xml(self) -> et.Element:

        pass

    def from_xml(self, e: et.Element):

        pass

    def set_parameters(self, set: ParameterSet):

        self.p = set
