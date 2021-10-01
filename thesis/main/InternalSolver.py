"""
@author: Lukas Kiwitz
"""

from abc import ABC, abstractmethod

from thesis.main.ParameterSet import ParameterSet


class InternalSolver(ABC):

    @abstractmethod
    def __init__(self):
        self.name = "InternalSolver"
        pass

    @abstractmethod
    def on_type_change(self, p: ParameterSet, replicat_index, entity=None): pass

    @abstractmethod
    def step(self, t1: float, t2: float, dt: float, p: ParameterSet, entity=None, **kwargs):
        return p
