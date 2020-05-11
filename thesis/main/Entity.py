#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 12:22:13 2019

@author: Lukas Kiwitz
"""

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
        self.bc_list = []
        self.id: int = 0
        self.change_type = ""

    def set_internal_solver(self, solver: InternalSolver) -> None:

        """

        set internal sovler for this entity

        :param solver:
        :return:

        """

        self.internal_solver = solver

    def get_BC(self, field_quantity: str) -> BC:

        """
        returns boundary condition for field_quantity

        :param field_quantity:
        :return:

        """

        for i in self.bc_list:
            if i.field_quantity == field_quantity:
                return i

    def update_bcs(self, p = None) -> None:

        """

        updates elements of self.bc_list
        :return:

        """
        if p or hasattr(self,"p"):
            self.fieldQuantities = []
            for i in self.bc_list:
                fq = i.field_quantity
                self.fieldQuantities.append(fq)

                if p:
                    bc.p = p
                else:
                    if hasattr(i, "p"):
                        i.p.update(self.p, override = True)
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
            self.internal_solver.step(t, dt, self.p, entity=self)

    def getState(self, key="q") -> float:

        """

        returns values of self.p[key]
        TODO throw exception
        :param key:
        :return:

        """

        if key in self.p and (type(self.p[key]) == float):
            return self.p[key]
        else:
            return float(0)


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
        self.center: List[float] = center
        self.radius = radius
        self.bc_list: List[BC] = bc_list

    def get_surface_area(self):

        rho = self.p.get_physical_parameter("rho", "rho").get_in_sim_unit()
        return (4 * np.pi * rho ** 2)

    def getSubDomain(self) -> CellSubDomain:
        """
        returns cell subdomain
        :return:
        """
        return SD.CellSubDomain(self.center, self.radius)

    def getCompiledSubDomain(self) -> fcs.CompiledSubDomain:
        """
        gets compiled subdomain

        :return:
        """
        return fcs.CompiledSubDomain(
            "on_boundary && abs((sqrt(pow(x[0]-c0,2)+pow(x[1]-c1,2)+pow(x[2]-c2,2))-r) <= 10e-2)",
            c0=self.center[0], c1=self.center[1], c2=self.center[2], r=self.radius)

    def set_cell_type(self, cell_type: CellType, internal_solver: InternalSolver) -> None:
        """

        sets this cells type from template object
        :param cell_type: template object
        :param kwargs:
        :return:

        """
        self.type_name = cell_type.name
        if cell_type.internal_solver:
            self.set_internal_solver(internal_solver())
        else:
            self.set_internal_solver(None)

        # miscs = ParameterCollection("misc", [
        #     MiscParameter("type_name", self.type_name),
        #     # MiscParameter("center", self.center)
        # ])
        self.p.update(cell_type.p, override=True)
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
        pass

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

    def __init__(self, center: List[float], radius: float, p: Dict, bc_list: List[OuterBC]):
        self.center = center
        self.radius = radius
        self.p = p
        self.bc_list = bc_list
        self.subdomain_dict = self.__compileSubdomains()
        super().__init__()

    def __compileSubdomains(self) -> Dict:
        subdomain_dict = {}
        for i, o in enumerate(self.bc_list):
            if isinstance(o, bc.OuterBC):
                e = CompiledSphere(o.expr, o, self.p)
                e.p = self.p
                e.field_quantity = o.field_quantity
                if o.expr not in subdomain_dict.keys():
                    subdomain_dict[o.expr] = [e]
                else:
                    subdomain_dict[o.expr].append(e)
        return subdomain_dict

    def getSubDomains(self, **kwargs):
        subdomains = []
        for i, o in enumerate(self.subdomain_dict.values()):
            if "field_quantity" in kwargs:
                for e in o:
                    if e.field_quantity == kwargs["field_quantity"]:
                        subdomains.append({"entity": e, "patch": i + 1})
            else:
                subdomains.append({"entity": o[0], "patch": i + 1})
        return subdomains

    def getSubDomain(self):
        return SD.OuterSphere(self.center, self.radius)


class DomainCube(DomainEntity):

    def __init__(self, p1, p2, bc_list, **kwargs):
        self.p1 = p1
        self.p2 = p2
        # self.p = p

        self.bc_list = bc_list
        self.subdomainDict = self.__compileSubdomains()
        super().__init__(**kwargs)

    def __compileSubdomains(self):
        subdomainDict = {}
        for i, o in enumerate(self.bc_list):
            if isinstance(o, bc.OuterBC):
                e = CompiledCube(o.expr, o, self)
                e.field_quantity = o.field_quantity
                if o.expr not in subdomainDict.keys():
                    subdomainDict[o.expr] = [e]
                else:
                    subdomainDict[o.expr].append(e)
        return subdomainDict

    def getSubDomains(self, **kwargs):
        subdomains = []
        for i, o in enumerate(self.subdomainDict.values()):
            if "field_quantity" in kwargs:
                for e in o:
                    if e.field_quantity == kwargs["field_quantity"]:
                        subdomains.append({"entity": e, "patch": i + 1})
            else:
                subdomains.append({"entity": o[0], "patch": i + 1})
        return subdomains

    def getSubDomainGeometry(self):
        return SD.OuterCube(self.p1, self.p2)

    def update_bcs(self, p = None) -> None:

        if p or hasattr(self, "p"):
            for k, v in self.subdomainDict.items():
                for i in v:
                    if p:
                        if hasattr(i.bc, "p"):
                            i.bc.p.update(p, override = False)
                        else:
                            i.bc.p = p

    def apply_sample(self, outer_domain_dict) -> None:

        for k, v in self.subdomainDict.items():
            for i in v:
                name = i.bc.name
                if name in outer_domain_dict.keys():

                    p = outer_domain_dict[name]

                    if hasattr(i.bc, "p"):
                        i.bc.p.update(p, override = True)
                    else:
                        i.bc.p = p


class CompiledCube(DomainCube, Entity):

    def __init__(self, expr, bc, parent):
        self.bc = bc
        self.parent = parent
        # self.p = self.bc.p
        self.expr = expr

    def getSubDomain(self):
        box = "near(x[0],{p1x}) || near(x[0],{p0x}) || near(x[1],{p1y}) || near(x[1],{p0y}) || near(x[2],{p1z}) || near(x[2],{p0z})".format(
            p1x=self.parent.p1[0],
            p0x=self.parent.p2[0],
            p1y=self.parent.p1[1],
            p0y=self.parent.p2[1],
            p1z=self.parent.p1[2],
            p0z=self.parent.p2[2],
        )
        return fcs.CompiledSubDomain(self.expr + "&&(" + box + ") && on_boundary")

    def get_BC(self, field_quantity):
        # self.bc.p = self.p
        return self.bc

    def get_surface_area(self):
        p = self.bc.p.get_physical_parameter("norm_area", "geometry")
        if p:
            return p.get_in_sim_unit()
        else:
            return None


class CompiledSphere(DomainSphere, Entity):

    def __init__(self, expr, bc, p):
        self.bc = bc
        self.p = p
        self.expr = expr

    def getSubDomain(self):
        box = "abs((sqrt(pow(x[0]-{c0},2)+pow(x[1]-{c1},2)+pow(x[2]-{c2},2))-{r}))<= 10e-2".format(c0=0, c1=0, c2=0,
                                                                                                   r=self.p["radius"])
        return fcs.CompiledSubDomain(self.expr + "&&(" + box + ") && on_boundary")

    def get_BC(self, field_quantity):
        self.bc.p = self.p
        return self.bc
