from typing import Dict


class CellType:

    def __init__(self, p: Dict, solver_class: type, type_name: str):
        if not type(solver_class) == type:
            raise Exception("solver_class must be a class name")
        self.internal_solver: type = solver_class
        self.p = p
        self.name = type_name

