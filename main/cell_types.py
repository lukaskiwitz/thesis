from copy import deepcopy
from CellType import CellType
from InternalSolver import InternalSolver, RuleBasedSolver
from scipy.constants import N_A

def make_cell_type(p,update,solver,name):
    p_temp = deepcopy(p)
    p_temp.update(update)
    p_temp["type_name"]=name
    return CellType(p_temp,solver,name)

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

cell_type_1 = make_cell_type(p,{
        "q_il2": p["q_h"],
        "R_il2": p["R_l"],
        "q_il6": p["q_l"],
        "R_il6": p["R_l"],
        "R_infg":0,
        "q_infg":0,


    }, InternalSolver, "type_1")

Tn = make_cell_type(p,{
        "q_il2": p["q_h"],
        "R_il2": p["R_l"],
        "q_il6": 0,#p["q_l"],
        "R_il6": p["R_l"],
        "R_infg":p["R_l"],
        "q_infg":0,
        "type_int":1
}, RuleBasedSolver, "Tn")

Tfh = make_cell_type(p,{
        "q_il2": p["q_h"],
        "R_il2": p["R_l"],
        "q_il6": p["q_h"],
        "R_il6": p["R_l"],
        "R_infg":0,
        "q_infg":0,
        "type_int":2
}, RuleBasedSolver, "Tfh")

Th1 = make_cell_type(p,{
        "q_il2": 0,#p["q_h"],
        "R_il2": 0,#p["R_l"],
        "q_il6": 0,#p["q_l"],
        "R_il6": 0,#p["R_l"],
        "q_infg": p["q_h"],
        "R_infg": p["R_l"],
        "type_int":3
}, RuleBasedSolver, "Th1")