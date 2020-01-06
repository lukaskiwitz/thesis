from copy import deepcopy
from CellType import CellType
from InternalSolver import InternalSolver
from scipy.integrate import odeint
from scipy.constants import N_A
import numpy as np

def make_cell_type(p,update,name):
    p_temp = deepcopy(p)
    p_temp.update(update)

    p_temp["type_name"]=name
    return CellType(p_temp,name)


x = [0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8,2, 2.2, 2.4, 2.6, 2.8]
y = [-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8]
z = y

# x = [-0.4, -0.2, 0, 0.2, 0.4]
# y = x
# z = y

# d_x = x[-1]+0.2
# d_y = y[-1]+0.2
# d_z = z[-1]+0.2

p_d = {"x":x,
       "y":y,
       "z":z,
       "margin":0.2
       }
p_sim = {#default values
        "R_il2":0,
        "q_il2":0,
        "R_il6":0,
        "q_il6":0,
        "R_infg":0,
        "q_infg":0,
        }
p = {
         "k_on":10e9*111.6/60**2,#111.6 per hour
         "rho": 0.05,#mu
         "D":(10**0.5*0.01)**2,#mu² per s
         "R_h":400*N_A**-1*10e9,
         "R_l":10*N_A**-1*10e9,
         "kd":0,#0.1/(60*2),
         "q_h":10*N_A**-1*10e9,
         "q_l":1*N_A**-1*10e9,
         "fraction":0.1
         }


Tn = make_cell_type(p,{
        "q_il2": p["q_h"],
        "R_il2": p["R_l"],
        "q_il6": 0,#p["q_l"],
        "R_il6": p["R_l"],
        "R_infg":p["R_l"],
        "q_infg":0,
        "type_int":1
} ,"Tn")

Tfh = make_cell_type(p,{
        "q_il2": p["q_h"],
        "R_il2": p["R_l"],
        "q_il6": p["q_h"],
        "R_il6": p["R_l"],
        "R_infg":0,
        "q_infg":0,
        "type_int":2
}, "Tfh")

Th1 = make_cell_type(p,{
        "q_il2": 0,#p["q_h"],
        "R_il2": 0,#p["R_l"],
        "q_il6": 0,#p["q_l"],
        "R_il6": 0,#p["R_l"],
        "q_infg": p["q_h"],
        "R_infg": p["R_l"],
        "type_int":3
}, "Th1")

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
        self.transition = (False, 0, Tn)

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
                    self.transition = (True, draw, Tfh)
                elif p["surf_c_infg"]*10e9 > 0.15:  #infg pos
                        self.transition = (True, draw, Th1)

        return p

Th1.internal_solver = RuleBasedSolver
Tfh.internal_solver = RuleBasedSolver
Tn.internal_solver = RuleBasedSolver