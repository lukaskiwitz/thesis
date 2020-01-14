"""
@author: Lukas Kiwitz
"""


class InternalSolver:
    def __init__(self):
        self.name = "InternalSolver (dummy)"
        pass
    def step(self,t,dt,p,**kwargs):
        return p


