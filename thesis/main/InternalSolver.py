"""
@author: Lukas Kiwitz
"""

from abc import ABC, abstractmethod

from thesis.main.ParameterSet import ParameterSet


class InternalSolver(ABC):
    name = "InternalSolver"

    @abstractmethod
    def on_type_change(self, p: ParameterSet, replicat_index, entity=None): pass

    @abstractmethod
    def step(self, t1: float, t2: float, dt: float, p: ParameterSet, entity=None, **kwargs) -> ParameterSet:
        return p
