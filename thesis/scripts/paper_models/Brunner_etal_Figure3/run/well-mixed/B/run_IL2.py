import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
from copy import deepcopy
import os

def EC50_calculation(E_max, E_min, k, N, R):
    return (E_max * k ** N + E_min * R ** N) / (k ** N + R ** N)

def I2(R_resp, K_D, q, f_Tsec, f_Treg, k_endo, R_Tsec, R_Treg, N_cells, kd):
    a = q * N_cells * f_Tsec
    b = k_endo  * N_cells * (R_resp * (1 - (f_Tsec + f_Treg)) + R_Tsec * f_Tsec + R_Treg * f_Treg)
    return (np.sqrt(kd ** 2 * K_D ** 2 + 2 * kd * K_D * (a + b) + (a - b) ** 2) - kd * K_D + a - b) / (2 * kd)

def I3(R_resp, K_D, q, f_Tsec, f_Treg, k_endo, R_Tsec, R_Treg):
    a = q * f_Tsec
    b = k_endo * (R_resp * (1 - (f_Tsec + f_Treg)) + R_Tsec * f_Tsec + R_Treg * f_Treg)
    return -(a * K_D)/(a - b)

def I_lin(R_resp, q, f_Tsec, f_Treg, k_on, R_Tsec, R_Treg, N_cells, kd):
    A = q * f_Tsec * N_cells
    B = k_on / dr.molar_to_molecules(1) * R_resp * (1 - (f_Tsec + f_Treg)) * N_cells \
        + k_on / dr.molar_to_molecules(1) * R_Tsec * f_Tsec * N_cells \
        + k_on / dr.molar_to_molecules(1) * R_Treg * f_Treg * N_cells\
        + kd
    return A/B


from parameters import p
from driver import ODEdriver
from thesis.cellBehaviourUtilities.cell_solver import kineticSolver

p["satu"] = True

chain_length = 6

save_df = True
hdd = "/extra2" if os.path.exists("/extra2") else "/extra"
base_path = hdd + "/brunner/paper_models/ODE/saturated/kinetics/gamma_"


p["R_Th"] = 1.5e3
p["IL-2_sigma"] = 1
p["IL-2_Treg_fraction"] = 0
p["kd"] = 0.1/3600
p["half_time"] = 0.1

for gamma in [10]:
    save_df_path = base_path + str(gamma) + "/"
    dr = ODEdriver(p)
    p["IL-2_KD"] = dr.molar_to_molecules(7.423e-12)
    KD = p["IL-2_KD"]

    # time_range = dr.exp_time_range(dt=3600, max_T=250, length=100)
    time_range = dr.exp_time_range(dt=3600, max_T=4 * 24, length=60)

    repeats = 1

    # scan_names = ["R_Th"] # ["D", "R_Th", "eta", "q"]
    scan_names = ["IL-2_Tsec_fraction"]
    # a = np.linspace(0.01, 0.4, 18)
    # scan_values = np.around(a , 3)
    scan_values = [0.1]

    scan_indices = np.array([x for x in range(len(scan_names) * len(scan_values))])
    length = len(scan_indices)

    cell_df = pd.DataFrame()
    global_df = pd.DataFrame()

    cell_df_dict_list = []
    ########################################################################################################################
    base_p = deepcopy(p)
    # systemic_R = 0.1 * p["R_Tsec"] + 0.9 * p["R_Th"]

    for r,rep in enumerate(range(repeats)):
        print("repeat no. ", rep)
        result = np.empty(length)
        std = np.empty(length)
        act = np.empty(length)
        for sn, scan_name in enumerate(scan_names):
            for sv, scan_value in enumerate(scan_values):
                print(sv + 1, "/", len(scan_values))
                p = deepcopy(base_p)
                scan_index = sn * len(scan_values) + sv
                p[scan_name] = scan_value

                # Th_R =  (systemic_R - p["IL-2_Tsec_fraction"] * 100) / (1 - p["IL-2_Tsec_fraction"])
                varied_Rs = dr.set_R_lognorm(p["R_Th"], p["IL-2_sigma"] * p["R_Th"], length=int(p["N_cells"] * (1 - p["IL-2_Treg_fraction"] - p["IL-2_Tsec_fraction"])))
                states = dr.set_R_lognorm_states(chain_length, p["R_Th"], p["IL-2_sigma"] * p["R_Th"], length=int(p["N_cells"] * (1 - p["IL-2_Treg_fraction"] - p["IL-2_Tsec_fraction"])), write_draws = False)
                starting_R = deepcopy(states)
                q = p["q"] #* 0.1 / p["IL-2_Tsec_fraction"]

                k_mins = deepcopy(varied_Rs * p["eta"])
                k_maxes = k_mins * gamma
                if gamma < 1:
                    p["K_m"] = 0.5
                elif gamma > 1:
                    p["K_m"] = 0.5
                for t, time in enumerate(time_range):
                    print(t+1, "/", len(time_range))
                    dt = time - time_range[t-1]
                    dt = dt if dt > 0 else 0

                    if p["satu"] == False:
                        c = dr.molecules_to_molar(
                            I_lin(states[:,-1].mean(), q, p["IL-2_Tsec_fraction"], p["IL-2_Treg_fraction"], p["k_on"], p["R_Tsec"], p["R_Treg"], 1,
                            p["kd"])) * 1e12
                    else:
                        c = dr.molecules_to_molar(I2(states[:,-1].mean(), KD, q, p["IL-2_Tsec_fraction"], p["IL-2_Treg_fraction"], p["k_endo"], p["R_Tsec"], p["R_Treg"], p["N_cells"], p["kd"])) * 1e9
                        # c = dr.molecules_to_molar(I3(states[:,-1].mean(), KD, q, p["IL-2_Tsec_fraction"], p["IL-2_Treg_fraction"], p["k_endo"], p["R_Tsec"], p["R_Treg"])) * 1e9

                    assert c >= 0
                    for n in range(int(p["N_cells"] * (1 - p["IL-2_Treg_fraction"] - p["IL-2_Tsec_fraction"]))):
                        cell_df_dict_list.append({
                            "id": n,
                            "IL-2_Tsec_fraction": p["IL-2_Tsec_fraction"],
                            "IL-2_q": p["q"],
                            "IL-2_sigma": p["IL-2_sigma"],
                            "IL-2_R": states[n,-1],
                            "IL-2_k_endo": p["k_endo"],
                            "IL-2_KD": p["IL-2_KD"],
                            "IL-2_surf_c": c,
                            "IL-2_gamma": gamma,
                            "misc_pos_half_time": p["half_time"],
                            "type_name": "Th",
                            "time_index": t,
                            "time": time,
                            "replicat_index": r,
                            "scan_value": sv,
                            "scan_index": scan_index,
                            "scan_name": scan_name,
                            "R_start_R_start": starting_R[n][-1]})

                    if dt != 0:
                        EC50 = dr.EC50_calculation(E_max=125e-12, E_min=0, k=860, N=1.5, R=states[:,-1]) * 1e9
                        for n in range(len(states)):
                            states[n] = kineticSolver.c_response_delayed(kineticSolver,
                                                                      p,
                                                                      p["K_m"],
                                                                      p["N_R"],
                                                                      p["eta"],
                                                                      k_maxes[n],
                                                                      k_mins[n],
                                                                      c,
                                                                      gamma,
                                                                      time,
                                                                      dt,
                                                                      list(states[n]),
                                                                      p["half_time"],
                                                                      EC50=EC50[n],
                                                                      pSTAT5_response=True)[1]
                            assert np.isnan(states[n][0]) == False

    cell_df = pd.DataFrame(cell_df_dict_list)
    if save_df == True:
        try:
            # global_df.to_hdf(save_df_path + "global_df.h5", "w")
            cell_df.to_hdf(save_df_path + "cell_df.h5", "w")
        except OSError:
            os.mkdir(save_df_path)
            # global_df.to_hdf(save_df_path + "global_df.h5", "w")
            cell_df.to_hdf(save_df_path + "cell_df.h5", "w")

        cell_df.to_hdf(save_df_path + "cell_df.h5", "w")
