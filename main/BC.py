#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  8 14:17:07 2019

@author: Lukas Kiwitz
"""
from typing import Dict, Callable

import fenics as fcs


class BC:

    def __init__(self, **kwargs: Dict) -> None:

        self.field_quantity = ""
        if "field_quantity" in kwargs:
            self.field_quantity = kwargs["field_quantity"]


class Integral(BC):
    def __init__(self, q: Callable, **kwargs: Dict) -> None:
        self.value = None
        self.p = {}
        self.q = q
        super().__init__(**kwargs)

    def get_BC(self, u: fcs.Function, p_update: Dict) -> object:
        p_temp = self.p
        p_temp.update(p_update)
        return self.q(u, p_temp)


class DirichletBC(BC):

    def __init__(self, value: object, **kwargs: Dict) -> None:
        self.degree = 1
        self.value = value
        super().__init__(**kwargs)

    def get_BC(self, V: fcs.FunctionSpace, boundary_markers: fcs.MeshFunction, patch: int) -> fcs.DirichletBC:
        value = fcs.Expression(str(self.value), degree=self.degree)
        bc = fcs.DirichletBC(V, value, boundary_markers, patch)
        return bc


class OuterBC(BC):
    def __init__(self, **kwargs: Dict) -> None:
        super().__init__(**kwargs)


class OuterDirichletBC(OuterBC,DirichletBC):

    def __init__(self, value: object, expr: str, **kwargs: Dict) -> None:
        self.degree = 1
        self.expr = expr
        self.value = value
        super().__init__(**kwargs)

    def get_BC(self, V: fcs.FunctionSpace, boundary_markers: fcs.MeshFunction, patch: int) -> fcs.DirichletBC:
        value = fcs.Expression(str(self.value), degree=self.degree)
        bc = fcs.DirichletBC(V, value, boundary_markers, patch)
        return bc


class OuterIntegral(Integral, OuterBC):

    def __init__(self, q: Callable, expr: str, **kwargs) -> None:
        self.value = None
        self.expr = expr
        self.q = q
        super().__init__(q,**kwargs)


