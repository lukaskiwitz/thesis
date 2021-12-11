"""
@author: Lukas Kiwitz
"""

from abc import ABC, abstractmethod

from thesis.main.ParameterSet import ParameterSet
from thesis.main.SimComponent import SimComponent


# module_logger = logging.getLogger(__name__)


class InternalSolver(ABC, SimComponent):
    name = "InternalSolver"

    def __init__(self):
        super(InternalSolver, self).__init__()

    @abstractmethod
    def on_type_change(self, p: ParameterSet, replicat_index, entity=None): pass

    @abstractmethod
    def step(self, t1: float, t2: float, dt: float, p: ParameterSet, entity=None, **kwargs) -> ParameterSet:
        return p
