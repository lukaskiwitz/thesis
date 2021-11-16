#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 12:22:13 2019

@author: Lukas Kiwitz
"""

from abc import *
from typing import List, Dict

import fenics as fcs
import numpy as np
import thesis.main.BC as bc
import thesis.main.MySubDomain as SD
from thesis.main.BC import BC, OuterBC
from thesis.main.EntityType import CellType
from thesis.main.InternalSolver import InternalSolver
from thesis.main.MySubDomain import CellSubDomain
from thesis.main.ParameterSet import ParameterSet
from thesis.main.TaskRecord import TaskRecord


class Entity:
    """

    :var: name: entity name
    :vartype name: str

    :var: fieldQuantities: list of field.field_quantity this entity interacts with
    :vartype fieldQuantities: List[str]

    :var: internal_solver: internal_solver for this entity
    :vartype internal_solver: InternalSolver

    :var: p: parameter Dict for this entity
    :vartype p: Dict

    :var: id: numeric identifier (uniqueness not strictly enforced)
    :vartype id: int

    :var bc_list: list of boundary contitions for this entity
    :vartype bc_list: List[BC]

    """

    def __init__(self):
        self.name: str = "default"
        self.fieldQuantities: List[str] = []
        self.internal_solver: InternalSolver = None
        # self.bc_list = []
        self.id: int = 0
        self.change_type = ""
        self.type_name = None
        self.task_record = TaskRecord("Entity")
        self.interactions = []
        self.p = ParameterSet("Entity_dummy", [])

    def move(self, dt):
        pass

    def set_internal_solver(self, solver: InternalSolver) -> None:

        """

        set internal sovler for this entity

        :param solver:
        :return:

        """

        self.internal_solver = solver

    def get_interaction(self, field_quantity: str) -> BC:

        """
        returns boundary condition for field_quantity

        :param field_quantity:
        :return:

        """

        for i in self.interactions:
            if i.field_quantity == field_quantity:
                return i

    def update_bcs(self, p=None) -> None:

        """

        updates elements of self.bc_list
        :return:

        """
        if p or hasattr(self, "p"):

            for i in self.interactions:
                fq = i.field_quantity

                if p:
                    bc.p = p
                else:
                    if hasattr(i, "p"):
                        i.p.update(self.p, overwrite=True)
                    else:
                        i.p = self.p

    def step(self, t: float, dt: float) -> None:

        """
        advances internal solver by dt
        :param t: absolute time
        :param dt: time step length
        :return:
        """

        if not self.internal_solver == None:
            self.internal_solver.step(t, t + dt, dt, self.p, entity=self)

    def getState(self, parameter_name="q", field_quantity="il2", in_post=True) -> float:

        return 0
        """

        returns values of self.p[key]
        TODO throw exception
        :param key:
        :return:

        """
        parameter = self.p.get_physical_parameter_by_field_quantity(parameter_name, field_quantity)
        if parameter is None:
            return float(0)
        else:
            if in_post:
                return parameter.get_in_post_unit()
            else:
                return parameter.get_in_sim_unit()

        # if key in self.p and (type(self.p[key]) == float):
        #     return self.p[key]
        # else:
        #     return float(0)


class Cell(Entity):
    """

    :var: center: cell center
    :vartype center: List[float]

    :var radius: cell radius
    :vartype radius: float

    :var cell_type: reference to instance of CellType
    :vartype CellType


    """

    def __init__(self, center: List[float], radius: float, bc_list: List[BC]) -> None:
        """

        :param center:
        :param radius:
        :param bc_list:
        """
        super().__init__()
        self.p = ParameterSet("Cell_dummy", [])
        self.interactions = []
        self.center: List[float] = center
        self.offset = [0, 0, 0]

        self.velocity = np.zeros(np.shape(center))

        self.radius = radius
        # self.bc_list: List[BC] = bc_list

    def move(self, dt):

        self.offset[0] += self.velocity[0] * dt
        self.offset[1] += self.velocity[1] * dt
        self.offset[2] += self.velocity[1] * dt

    def move_real(self, dt, bouding_box):

        from thesis.main.my_debug import message

        def inside(p1, p2, i, x, m):

            r = [p1[i], p2[i]]
            x = x[i]

            if (r[0] + m < x) and (r[1] - m > x):
                return True
            else:
                message("cell outside of bounding box")
                return False

        from copy import deepcopy

        old = deepcopy(self.center)

        p1 = bouding_box.p1
        p2 = bouding_box.p2
        m = self.radius * 2

        self.center[0] += self.velocity[0] * dt
        if not inside(p1, p2, 0, self.center, m):
            message("setting x to " + str(old[0]))
            self.center[0] = old[0]

        self.center[1] += self.velocity[1] * dt
        if not inside(p1, p2, 1, self.center, m):
            message("setting y to " + str(old[1]))
            self.center[1] = old[1]

        self.center[2] += self.velocity[2] * dt
        if not inside(p1, p2, 2, self.center, m):
            message("setting z to " + str(old[2]))
            self.center[2] = old[2]

    def get_surface_area(self):

        rho = self.p.get_physical_parameter("rho", "rho").get_in_sim_unit()
        return (4 * np.pi * rho ** 2)

    def get_subdomain(self) -> CellSubDomain:
        """
        returns cell subdomain
        :return:
        """
        return SD.CellSubDomain(self.center, self.radius)

    def get_compiled_subdomain(self) -> fcs.CompiledSubDomain:
        """
        gets compiled subdomain

        :return:
        """
        return fcs.CompiledSubDomain(
            "on_boundary && abs((sqrt(pow(x[0]-c0,2)+pow(x[1]-c1,2)+pow(x[2]-c2,2))-r) <= 10e-2)",
            c0=self.center[0], c1=self.center[1], c2=self.center[2], r=self.radius)

    def set_cell_type(self, cell_type: CellType, internal_solver: InternalSolver, replicat_index: int) -> None:
        """

        sets this cells type from template object
        :param cell_type: template object
        :param kwargs:
        :return:

        """
        self.type_name = cell_type.name
        self.interactions = []
        if internal_solver is None:
            self.set_internal_solver(None)
        else:
            if issubclass(internal_solver, InternalSolver):
                self.set_internal_solver(internal_solver())
            else:
                self.set_internal_solver(None)

            self.internal_solver.on_type_change(self.p, replicat_index, entity=self)
        # miscs = ParameterCollection("misc", [
        #     MiscParameter("type_name", self.type_name),
        #     # MiscParameter("center", self.center)
        # ])
        self.p.update(cell_type.p, overwrite=True)
        for interaction_template in cell_type.interactions:
            interaction = interaction_template.get_interaction()
            self.interactions.append(interaction)

        # self.p.update(ParameterSet("dummy", [miscs]),override=True)

    def change_entity_type(self, type_name: str):
        """
        schedules type change for next timestep
        """
        self.change_type = type_name


class DomainEntity(Entity):
    """

    base class

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def getState(self, parameter_name="q", field_quantity="il2", in_post=True) -> float:
        return 0

    @abstractmethod
    def apply_sample(self, outer_domain_dict) -> None:
        raise NotImplementedError

    @abstractmethod
    def get_subdomains(self, **kwargs) -> List:
        raise NotImplementedError

    def update_bcs(self, p=None) -> None:

        if p or hasattr(self, "p"):
            for k, v in self.subdomain_dict.items():
                for i in v:
                    if p:
                        for ii in i.interactions:
                            if hasattr(ii, "p"):
                                ii.p.update(p, overwrite=False)
                            else:
                                ii.p = p

    def apply_sample(self, outer_domain_dict) -> None:

        for k, v in self.subdomain_dict.items():
            for i in v:
                for ii in i.interactions:
                    name = ii.name
                    if name in outer_domain_dict.keys():
                        p = outer_domain_dict[name]
                        if hasattr(ii, "p"):
                            ii.p.update(p, overwrite=True)
                        else:
                            ii.p = p


class DomainSphere(DomainEntity):
    """

    define spherical outside domain

    :var center:
    :vartype center: List[float]

    :var radius:
    :vartype radius: float

    :var p:
    :vartype p: Dict

    :var bc_list:
    :vartype bc_list: List[OuterBC]

    :var subdomain_dict: dict constructed from bc_list
    :vartype subdomain_dict: Dict
    """

    def __init__(self, center: List[float], radius: float, bc_list: List[OuterBC]):
        self.center = center
        self.radius = radius

        self.bc_list = bc_list
        self.subdomain_dict = self.__compile_subdomains()
        super().__init__()

    def __compile_subdomains(self) -> Dict:

        subdomain_dict = {}
        for i, o in enumerate(self.bc_list):
            if isinstance(o, bc.OuterBC):
                e = CompiledSphere(o.expr, o, self)
                e.field_quantity = o.field_quantity
                if o.expr not in subdomain_dict.keys():
                    subdomain_dict[o.expr] = [e]
                else:
                    subdomain_dict[o.expr].append(e)
        return subdomain_dict

    def get_subdomains(self, **kwargs):
        subdomains = []
        for i, o in enumerate(self.subdomain_dict.values()):
            if "field_quantity" in kwargs:
                for e in o:
                    if e.field_quantity == kwargs["field_quantity"]:
                        subdomains.append({"entity": e, "patch": i + 1})
            else:
                subdomains.append({"entity": o[0], "patch": i + 1})
        return subdomains




class CompiledSphere(DomainSphere, Entity):

    def __init__(self, expr, bc, parent):
        self.interactions = [bc]
        self.parent = parent
        self.expr = expr

    def get_subdomain(self):
        c = self.parent.center

        box = "abs((sqrt(pow(x[0]-{c0},2)+pow(x[1]-{c1},2)+pow(x[2]-{c2},2))-{r}))<= 10e-2".format(
            c0=c[0], c1=c[1], c2=c[2],
            r=self.parent.radius)
        return fcs.CompiledSubDomain(self.expr + "&&(" + box + ") && on_boundary")

    def get_BC(self, field_quantity):

        for i in self.interactions:
            if i.field_quantity == field_quantity:
                return i

    def get_surface_area(self):
        p = self.interactions[0].p.get_physical_parameter("norm_area", "geometry")
        if p:
            return p.get_in_sim_unit()
        else:
            return None



class DomainCube(DomainEntity):

    def __init__(self, p1, p2, interactions, **kwargs):
        self.p1 = p1
        self.p2 = p2

        self.interactions = interactions

        self.subdomain_dict = self.__compile_subdomains()
        super().__init__(**kwargs)

    def __compile_subdomains(self):
        subdomain_dict = {}

        for i, o in enumerate(self.interactions):
            if isinstance(o, bc.OuterBC):

                e = CompiledCube(o.expr, o, self)

                e.field_quantity = o.field_quantity
                if o.expr not in subdomain_dict.keys():
                    subdomain_dict[o.expr] = [e]
                else:
                    subdomain_dict[o.expr].append(e)
        return subdomain_dict

    def get_subdomains(self, **kwargs):
        subdomains = []
        for i, o in enumerate(self.subdomain_dict.values()):
            if "field_quantity" in kwargs:
                for e in o:
                    if e.field_quantity == kwargs["field_quantity"]:
                        subdomains.append({"entity": e, "patch": i + 1})
            else:
                subdomains.append({"entity": o[0], "patch": i + 1})
        return subdomains

    def get_subdomain_geometry(self):
        return SD.OuterCube(self.p1, self.p2)


class CompiledCube(DomainCube, Entity):

    def __init__(self, expr, bc, parent):
        self.interactions = [bc]
        self.parent = parent
        self.expr = expr

    def get_subdomain(self):

        box = "near(x[0],{p1x}) || near(x[0],{p0x}) || near(x[1],{p1y}) || near(x[1],{p0y}) || near(x[2],{p1z}) || near(x[2],{p0z})".format(
            p0x=self.parent.p1[0],
            p1x=self.parent.p2[0],
            p0y=self.parent.p1[1],
            p1y=self.parent.p2[1],
            p0z=self.parent.p1[2],
            p1z=self.parent.p2[2],
        )
        expr = self.expr.format(
            p0x=self.parent.p1[0],
            p1x=self.parent.p2[0],
            p0y=self.parent.p1[1],
            p1y=self.parent.p2[1],
            p0z=self.parent.p1[2],
            p1z=self.parent.p2[2],
        )
        return fcs.CompiledSubDomain(expr + "&&(" + box + ") && on_boundary")

    def get_BC(self, field_quantity):

        for i in self.interactions:
            if i.field_quantity == field_quantity:
                return i

    def get_surface_area(self):
        p = self.interactions[0].p.get_physical_parameter("norm_area", "geometry")
        if p:
            return p.get_in_sim_unit()
        else:
            return None
