from copy import deepcopy
from typing import Dict


class EntityType():

    def __init__(self, p: Dict, name: str):

        self.p = deepcopy(p)
        self.name = name


class CellType(EntityType):

    def __init__(self, p: Dict, name: str, solver_name: str):

        self.p = deepcopy(p)
        self.name = name
        self.internal_solver = solver_name


