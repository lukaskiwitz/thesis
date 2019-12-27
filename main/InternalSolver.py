"""
@author: Lukas Kiwitz
"""
from scipy.integrate import odeint
from scipy.constants import N_A
import cell_types
import numpy as np

class InternalSolver:
    def __init__(self):
        pass
    def step(self,t,dt,p,**kwargs):
        return p

class ODE_Solver(InternalSolver):
    def __init__(self):
        pass
    def f(self,x,t,p):

        C = p["k_on"]*x[0]*p["surf_c_il2"]
        R = x[0]

        K = 1000 * N_A**-1*10e9
        v1 = 3000 *N_A**-1*10e9
        v0 = 150 *N_A**-1*10e9
        k_deg = 0.1
        n = 12
        # print(R/(N_A**-1*10e9))

        return p["i_C"]*((v1*(C**n)/(K**n+C**n)) - k_deg*R)

    def step(self,t,dt,p,entity=None):

        r = odeint(self.f,p["R_il2"],[0,dt],args=(p,))
        N_A ** -1 * 10e9
        p["R_il2"] = r[1][0]
        return p
class RuleBasedSolver(InternalSolver):

    def __init__(self):
        self.transition = (False,0,cell_types.Tn)

    def step(self,t,dt,p,entity=None):

        if entity.cell_type.name == "Tn":
            if self.transition[0]:
                if self.transition[1] <= 0:
                    entity.set_cell_type(self.transition[2])
                else:
                    self.transition = (True,self.transition[1]-dt,self.transition[2])
            elif np.random.rand(1) > 0.5:  # chance for type change; uniform distribution
                draw = np.random.normal(1,0.2)  # time to type change; normal distribution

                if p["surf_c_il2"]*10e9 < 0.7 and p["surf_c_il6"]*10e9 > 0.35:  # il2neg and il6pos
                    self.transition = (True,draw,cell_types.Tfh)
                elif p["surf_c_infg"]*10e9 > 0.15:  #infg pos
                        self.transition = (True,draw,cell_types.Th1)

        return p
