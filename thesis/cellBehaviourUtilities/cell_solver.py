import sys

sys.path.append("/home/brunner/thesis/thesis/main/")
# sys.path.append("/home/brunner/thesis/thesis/scenarios/")

import numpy as np
from scipy.constants import N_A
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

    def __init__(self):

        self.T = [0]
        super().__init__()

    def on_type_change(self, p, replicat_index, entity=None):
        pass

    def EC50_calculation(self, E_max,E_min,k,N,R):
        return (E_max * k ** N + E_min * R ** N) / (k ** N + R ** N)

    def c_response_delayed(self,p, Km,N,eta,k_max,k_min,il2,gamma,t1,dt, states, half_time, EC50, pSTAT5_response = True):
        if pSTAT5_response == True:
            pSTAT5 = il2 ** 3 / (EC50 ** 3 + il2 ** 3)
            try:
                p.get_physical_parameter("pSTAT5", "IL-2").set_in_post_unit(pSTAT5)
            except:
                pass  # in case this function is importet without sc or sm
        else:
            pSTAT5 = il2


        dict = {"N": N, "gamma": gamma, "eta": eta, "il2": il2, "pSTAT5": pSTAT5,
                "k_min": k_min, "k_max": k_max, "k": Km, "lambda": np.log(2)/(half_time*3600)}

        def func(t, y, dict, dummy=0):
            R0, R1, R2, R3, R4, R5 = y

            decay = dict["eta"]

            dR0_dt = ((dict["k_min"] * dict["k"] ** dict["N"] + dict["k_max"] * dict["pSTAT5"] ** dict["N"]) / (
                    dict["k"] ** dict["N"] + dict["pSTAT5"] ** dict["N"]) - decay * R0)

            dR1_dt = dict["lambda"] * R0 - dict["lambda"] * R1

            dR2_dt = dict["lambda"] * R1 - dict["lambda"] * R2

            dR3_dt = dict["lambda"] * R2 - dict["lambda"] * R3

            dR4_dt = dict["lambda"] * R3 - dict["lambda"] * R4

            dR5_dt = dict["lambda"] * R4 - dict["lambda"] * R5

            return [dR0_dt, dR1_dt, dR2_dt, dR3_dt, dR4_dt, dR5_dt]
        try:
            y0 = np.array(states.get_in_post_unit())/(N_A ** -1 * 1e9)
        except AttributeError:
            y0 = states

        t_space = np.linspace(t1, t1 + dt, 100)
        result = solve_ivp(func, t_span=[t1, t1 + dt], t_eval=t_space, y0=y0, method="LSODA", args=(dict, 0))
        return result.y[-1][-1], result.y[:,-1].tolist()

    def set_k(self, p, eta, gamma, initial_R):
        k_min = initial_R/(N_A ** -1 * 1e9) * eta
        k_max = k_min * gamma
        if gamma < 1: #neg
            Km = p.get_misc_parameter("Km_neg", "misc").get_in_post_unit()   # M
        elif gamma >= 1: #pos
            Km = p.get_misc_parameter("Km_pos", "misc").get_in_post_unit()   # M
        return k_min, k_max, Km

    def setup_initial(self, p):
        initial_R = p.get_physical_parameter("R", "IL-2").get_in_sim_unit()
        p.get_physical_parameter("R", "IL-2").set_in_sim_unit(initial_R)
        states = p.get_misc_parameter("states", "misc")
        states_values = states.get_in_post_unit()
        states_values = [initial_R for x in states_values]
        states.set_in_post_unit(states_values)
        return p

    def step(self,t1,t2,dt,p,entity=None):
        celltype = p.get_misc_parameter("name", "misc").get_in_post_unit()
        if celltype == "Tsec":
            try:
                start_t = p.get_misc_parameter("sec_start", "misc").get_in_post_unit()
                tmp_q = p.get_misc_parameter("tmp_q", "misc").get_in_post_unit()
                if t1 > start_t:
                    p.get_physical_parameter("q", "IL-2").set_in_post_unit(tmp_q)
                    # print("started")
                else:
                    pass
                    # p.get_physical_parameter("q", "IL-2").set_in_post_unit(0)
            except AttributeError:
                # print("tmp_q or sec_start not defined, skipping")
                pass
            # pass
        elif celltype != "Tsec" and celltype != "Tnaive":
            initial_R = p.get_physical_parameter("R_start", "R_start").get_in_sim_unit()
            try:
                gamma = p.get_physical_parameter("gamma", "IL-2").get_in_sim_unit()
            except Exception as e:
                try:
                    gamma = p.get_misc_parameter("gamma", "misc").get_in_sim_unit()
                except:
                    gamma = p.get_physical_parameter("gamma", "misc").get_in_sim_unit()
                if isinstance(gamma, str):
                    gamma = float(gamma)
            #get parameters from scan
            N = p.get_misc_parameter("hill_factor", "misc").get_in_post_unit()
            il2 = p.get_physical_parameter("surf_c", "IL-2").get_in_post_unit() * 1e-9 #post is nM, neeeded in units of Km --> M

            try:
                eta = p.get_physical_parameter("eta", "misc").get_in_post_unit()
            except Exception as e:
                eta = 1 / 72000

            R_il2 = p.get_physical_parameter("R", "IL-2").get_in_post_unit()

            k_min, k_max, Km = self.set_k(p, eta, gamma, initial_R)

            if self.T[0] == t1:  # set initial receptors
                p = self.setup_initial(p)
            elif gamma < 1 and celltype == "Treg":  # in negative fb Tregs = ILC = no feedback
                pass
            else:
                states = p.get_misc_parameter("states", "misc")
                EC50_k = p.get_misc_parameter("EC50_k", "misc").get_in_post_unit() if p.get_misc_parameter(
                    "EC50_k", "misc") is not None else 860
                EC50_N = p.get_misc_parameter("EC50_N", "misc").get_in_post_unit() if p.get_misc_parameter(
                    "EC50_N", "misc") is not None else 1.5

                EC50 = self.EC50_calculation(E_max=125e-12, E_min=0, k=EC50_k, N=EC50_N, R=R_il2)
                p.get_physical_parameter("pSTAT5", "misc").set_in_post_unit(il2 ** 3 / (EC50 ** 3 + il2 ** 3))
                p.get_physical_parameter("EC50", "EC50").set_in_post_unit(EC50)

                try:
                    half_time = p.get_misc_parameter("pos_half_time",
                                                     "misc").get_in_post_unit() if gamma >= 1 else p.get_misc_parameter(
                        "neg_half_time", "misc").get_in_post_unit()
                except:
                    print(p.get_as_dictionary())
                    print("halftime not defined")
                    exit()
                if celltype == "Treg":
                    k_max_factor = 1
                else:
                    k_max_factor = 1

                new_R, new_states = self.c_response_delayed(p, Km, N, eta, k_max * k_max_factor, k_min, il2, gamma,
                                                            t1, dt, states, half_time, EC50)
                states.set_in_post_unit([x * N_A ** -1 * 1e9 for x in new_states])

                p.get_physical_parameter("R", "IL-2").set_in_sim_unit(new_R * N_A ** -1 * 1e9)
        return p
