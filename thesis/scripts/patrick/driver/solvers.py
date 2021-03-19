import sys

sys.path.append("/home/brunner/thesis/thesis/main/")
# sys.path.append("/home/brunner/thesis/thesis/scenarios/")

import numpy as np
from scipy.constants import N_A
# from sympy import symbols, solve
from scipy.integrate import solve_ivp


from thesis.main.InternalSolver import InternalSolver

import mpi4py.MPI as MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

class kineticSolver(InternalSolver):

    """Controls the internal logic of cells. An instance is created for each cell.
    Persistent values should be stored in the cell ParameterSet, values attached to the InternalSolver
    object may get lost."""

    name = "kineticSolver"
    #
    # def __init__(self):
    #     self.il2_threshold = 0.0048
    def EC50_calculation(self, k_max,k_min,k,N,il2):
        return (k_max * k ** N + k_min * il2 ** N) / (k ** N + il2 ** N)

    def pSTAT5_response(self,c_0,EC50,N,eta,R_mean,k_max,k_min,il2,gamma,R_il2,t1,dt):
        typical_pSTAT5 = c_0 ** 3 / (c_0 ** 3 + EC50 ** 3)
        # k = solve((k_min * k_x ** N + k_max * typical_pSTAT5 ** N) / ((k_x ** N + typical_pSTAT5 ** N) * eta) - R_mean)[0]
        A = (typical_pSTAT5 ** N * (eta * R_mean - k_max))
        B = (k_min - eta * R_mean)
        k = (A / B) ** (1 / N)

        pSTAT5 = il2 ** 3 / (il2 ** 3 + EC50 ** 3)
        # print(pSTAT5)
        dict = {"N": N, "gamma": gamma, "eta": eta, "c_0": c_0, "pSTAT5": pSTAT5, "R_il2": R_il2, "R_mean": R_mean,
                "k_min": k_min, "k_max": k_max, "k": k}

        def func(t, y, p, dummy=0):
            R_il2 = y[0]
            dR_dt = ((p["k_min"] * p["k"] ** p["N"] + p["k_max"] * p["pSTAT5"] ** p["N"]) / (
                    p["k"] ** p["N"] + p["pSTAT5"] ** p["N"]) - p["eta"] * R_il2)
            return [dR_dt]

        y0 = [dict["R_il2"]]
        t_space = np.linspace(t1, t1 + dt, 100)
        result = solve_ivp(func, t_span=[t1, t1 + dt], t_eval=t_space, y0=y0, method="LSODA", args=(dict, 0))
        return pSTAT5, result.y[0][-1]

    def c_response(self,c_0,N,eta,R_mean,k_max,k_min,il2,gamma,R_il2,t1,dt):
        try:
            # k = solve((k_min * k_x ** N + k_max * c_0 ** N) / ((k_x ** N + c_0 ** N) * eta) - R_mean)[0]
            k = ((c_0 ** N * (eta * R_mean - k_max)) / (k_min - eta * R_mean)) ** (1 / N)
        except:  # if gamma ~ 1
            k = 1

        dict = {"N": N, "gamma": gamma, "eta": eta, "c_0": c_0, "il2": il2, "R_il2": R_il2, "R_mean": R_mean,
                "k_min": k_min, "k_max": k_max, "k": k}

        def func(t, y, p, dummy=0):
            R_il2 = y[0]
            dR_dt = ((p["k_min"] * p["k"] ** p["N"] + p["k_max"] * p["il2"] ** p["N"]) / (
                    p["k"] ** p["N"] + p["il2"] ** p["N"]) - p["eta"] * R_il2)
            return [dR_dt]

        y0 = [dict["R_il2"]]
        t_space = np.linspace(t1, t1 + dt, 100)
        result = solve_ivp(func, t_span=[t1, t1 + dt], t_eval=t_space, y0=y0, method="LSODA", args=(dict, 0))
        return result.y[0][-1]

    def step(self,t1,t2,dt,p,entity=None):
        if t1 > 0: # since pre_step calculates c_0 in the first step we can only start at the second step
            if p.get_physical_parameter("q", "IL-2").get_in_sim_unit() == 0:
                #get parameters from scan
                N = p.get_physical_parameter("hill_factor", "IL-2").get_in_post_unit()
                c_0 = p.get_physical_parameter("c0", "IL-2").get_in_post_unit() #1.90e-12 # M
                il2 = p.get_physical_parameter("surf_c", "IL-2").get_in_post_unit()*1e-9 #post is nM, neeeded in M

                gamma = p.get_physical_parameter("gamma", "IL-2").get_in_sim_unit()
                R_start = p.get_physical_parameter("R_start", "R_start").get_in_sim_unit()
                R_il2 = p.get_physical_parameter("R", "IL-2").get_in_sim_unit()
                eta = 1 / 72000  # 1/s

                R_mean = 1e4 * N_A ** -1 * 1e9
                k_min = R_mean / gamma * eta
                k_max = R_mean * gamma * eta

                if p.get_physical_parameter("pSTAT5_signal", "IL-2").get_in_post_unit() == True:
                    EC50 = self.EC50_calculation(k_max=4e-11, k_min=2e-12, k=1e4 * N_A ** -1 * 1e9, N=1.2, il2=il2)
                    pSTAT5,new_R = self.pSTAT5_response(c_0,EC50,N,eta,R_mean,k_max,k_min,il2,gamma,R_il2,t1,dt)
                    p.get_physical_parameter("EC50", "EC50").set_in_sim_unit(EC50)
                    p.get_physical_parameter("R", "IL-2").set_in_sim_unit(new_R)
                    p.get_physical_parameter("pSTAT5", "pSTAT5").set_in_sim_unit(pSTAT5)
                else:
                    new_R = self.c_response(c_0,N,eta,R_mean,k_max,k_min,il2,gamma,R_il2,t1,dt)
                    p.get_physical_parameter("R", "IL-2").set_in_sim_unit(new_R)
        return p