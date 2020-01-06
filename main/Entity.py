#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 12:22:13 2019

@author: Lukas Kiwitz
"""

from copy import deepcopy
from typing import List, Dict

import fenics as fcs

import BC as bc
import MySubDomain as SD
from BC import BC, OuterBC
from CellType import CellType
from InternalSolver import InternalSolver
from MySubDomain import CellSubDomain


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
        self.name: st = "default"
        self.fieldQuantities: List[str] = []
        self.internal_solver: InternalSolver = None
        self.bc_list = []
        self.p: Dict = {}
        self.id: int = 0

    def set_internal_solver(self, solver: InternalSolver) -> None:

        """

        set internal sovler for this entity

        :param solver:
        :return:

        """

        self.internalSolvers = solver

    def get_BC(self, fieldQuantity: str) -> BC:

        """
        returns boundary condition for fieldQuantity

        :param fieldQuantity:
        :return:

        """

        for i in self.bcList:
            if i.fieldQuantity == fieldQuantity:
                return i

    def update_bcs(self) -> None:

        """

        updates elements of self.bc_list
        :return:

        """

        self.fieldQuantities = []
        for i in self.bc_list:
            fq = i.fieldQuantity
            self.fieldQuantities.append(fq)
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

    # def log(self):
    #     return {}

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

    def __init__(self, center: List[float], radius: float, bcList: List[BC]) -> None:
        """

        :param center:
        :param radius:
        :param bcList:
        """

        self.center: List[float] = center
        self.radius = radius
        self.bc_list: List[BC] = bcList
        super().__init__()

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

    def set_cell_type(self, cell_type: CellType, p=None) -> None:
        """

        sets this cells type from template object
        :param cell_type: template object
        :param kwargs:
        :return:

        """
        self.cell_type = cell_type
        self.set_internal_solver(cell_type.internal_solver())
        self.p = deepcopy(p) if (not p == None) else deepcopy(cell_type.p)
        self.p["center"] = self.center

    # def log(self):
    #     p_out = self.p
    #     for i in p_out.keys():
    #         if "numpy.float64" == str(type(p_out[i])):
    #             p_out[i] = float(p_out[i])
    #         if "ndarray" in str(type(p_out[i])):
    #             p_out[i] = list(p_out[i])
    #     di = {"type": str(type(self)),
    #           "id": self.id,
    #           "name": self.name,
    #           "center": self.center,
    #           "radius": self.radius,
    #           "p": self.p
    #           }
    #     return di


class DomainEntity(Entity):
    """

    base class

    """

    def __init__(self, **kwargs):
        super().__init__()

    # def log(self):
    #     return {"type": str(type(self))}


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
                e.fieldQuantity = o.fieldQuantity
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
                    if e.fieldQuantity == kwargs["field_quantity"]:
                        subdomains.append({"entity": e, "patch": i + 1})
            else:
                subdomains.append({"entity": o[0], "patch": i + 1})
        return subdomains

    def getSubDomain(self):
        return SD.OuterSphere(self.center, self.radius)


class DomainCube(DomainEntity):

    def __init__(self, p1, p2, p, bcList):
        self.p1 = p1
        self.p2 = p2
        self.p = p

        self.bcList = bcList
        self.subdomainDict = self.__compileSubdomains()
        super().__init__()

    def __compileSubdomains(self):
        subdomainDict = {}
        for i, o in enumerate(self.bcList):
            if isinstance(o, bc.OuterBC):
                e = CompiledCube(o.expr, o, self)
                e.p = self.p
                e.fieldQuantity = o.fieldQuantity
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
                    if e.fieldQuantity == kwargs["field_quantity"]:
                        subdomains.append({"entity": e, "patch": i + 1})
            else:
                subdomains.append({"entity": o[0], "patch": i + 1})
        return subdomains

    def getSubDomainGeometry(self):
        return SD.OuterCube(self.p1, self.p2)


class CompiledCube(Entity):

    def __init__(self, expr, bc, parent):
        self.bc = bc
        self.parent = parent
        self.p = self.parent.p
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
        print("expr " + self.expr)
        return fcs.CompiledSubDomain(self.expr + "&&(" + box + ") && on_boundary")

    def get_BC(self, fieldQuantity):
        self.bc.p = self.p
        return self.bc


class CompiledSphere(Entity):

    def __init__(self, expr, bc, p):
        self.bc = bc
        self.p = p
        self.expr = expr

    def getSubDomain(self):
        box = "abs((sqrt(pow(x[0]-{c0},2)+pow(x[1]-{c1},2)+pow(x[2]-{c2},2))-{r}))<= 10e-2".format(c0=0, c1=0, c2=0,
                                                                                                   r=self.p["radius"])
        return fcs.CompiledSubDomain(self.expr + "&&(" + box + ") && on_boundary")

    def get_BC(self, fieldQuantity):
        self.bc.p = self.p
        return self.bc
