from thesis.main.Entity import DomainCube
import numpy as np

class MyDomainTemplate:

    def __init__(self):
        pass

    def get_domain(self):
        pass


class MyBoxDomainTemplate:

    def __init__(self):

        self.bc_list = []

    def get_domain(self, p, p1,p2):

        domain = DomainCube(p1, p2, self.bc_list)
        return domain

class BoundingBoxError(Exception): pass

class MyBoundingBoxTemplate:

    def __init__(self):

        self.bc_list = []

    def get_domain(self, p, entity_list):


        cell_position = np.array([c["entity"].center for c in entity_list])
        margin = p.get_misc_parameter("margin","geometry").get_in_sim_unit(type = float)

        if len(cell_position) == 0:
            raise BoundingBoxError

        p1 = np.min(cell_position, axis=0) - margin
        p2 = np.max(cell_position, axis=0) + margin


        domain = DomainCube(p1, p2, self.bc_list)

        return domain

