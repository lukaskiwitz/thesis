from typing import Dict
from copy import deepcopy
import lxml.etree as ET

class EntityType():
    pass


class CellType(EntityType):

    def __init__(self, p: Dict, name: str, solver_name: type):

        self.p = deepcopy(p)
        self.name = name
        self.internal_solver = solver_name


