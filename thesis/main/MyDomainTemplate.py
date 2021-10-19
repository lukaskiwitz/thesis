from typing import List

import numpy as np

from thesis.main.Entity import DomainCube, DomainSphere, Entity
from thesis.main.ParameterSet import ParameterSet


class BoundingBoxError(Exception): pass


class MyDomainTemplate:

    def __init__(self):
        pass

    def get_domain(self, p: ParameterSet, entity_list: List[Entity]):
        pass


class MyBoxDomainTemplate(MyDomainTemplate):

    def __init__(self,p1,p2):
        self.p1 = p1
        self.p2 = p2

        self.bc_list = []

    def get_domain(self, p: ParameterSet, entity_list: List[Entity]):
        domain = DomainCube(self.p1, self.p2, self.bc_list)
        return domain


class MyBoundingBoxTemplate(MyDomainTemplate):

    def __init__(self):
        self.bc_list = []

    def get_domain(self, p: ParameterSet, entity_list: List[Entity]):
        cell_position = np.array([c["entity"].center for c in entity_list])
        margin = p.get_misc_parameter("margin", "geometry").get_in_sim_unit(type=float)

        if len(cell_position) == 0:
            raise BoundingBoxError

        p1 = np.min(cell_position, axis=0) - margin
        p2 = np.max(cell_position, axis=0) + margin

        domain = DomainCube(p1, p2, self.bc_list)

        return domain


class MySphereDomainTemplate(MyDomainTemplate):

    def __init__(self,c,r):

        assert len(c) ==3

        self.c = c
        self.r = r
        self.bc_list = []

    def get_domain(self, p: ParameterSet, entity_list: List[Entity]):
        domain = DomainSphere(self.c, self.r, self.bc_list)

        return domain
