"""
@author: Lukas Kiwitz
"""


class InternalSolver:
    def __init__(self):
        self.name = "InternalSolver (dummy)"
        pass

    def step(self, t1, t2, dt, p, **kwargs):
        return p
