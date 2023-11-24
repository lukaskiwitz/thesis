import numpy as np
import pandas as pd
import os
from parameters import p
from driver import ODEdriver

def I2(R_resp, K_D, q, f_Tsec, f_Treg, k_endo, R_Tsec, R_Treg, N_cells, kd):
    a = q * N_cells * f_Tsec
    b = k_endo  * N_cells * (R_resp * (1 - (f_Tsec + f_Treg)) + R_Tsec * f_Tsec + R_Treg * f_Treg)
    return (np.sqrt(kd ** 2 * K_D ** 2 + 2 * kd * K_D * (a + b) + (a - b) ** 2) - kd * K_D + a - b) / (2 * kd)

def I_lin(R_resp, q, f_Tsec, f_Treg, k_on, R_Tsec, R_Treg, N_cells, kd):
    A = q * f_Tsec * N_cells
    B = k_on / dr.molar_to_molecules(1) * R_resp * (1 - (f_Tsec + f_Treg)) * N_cells \
        + k_on / dr.molar_to_molecules(1) * R_Tsec * f_Tsec * N_cells \
        + k_on / dr.molar_to_molecules(1) * R_Treg * f_Treg * N_cells\
        + kd
    return A/B

p["satu"] = True
Tregs = False

save_df = True
hdd = "/extra2" if os.path.exists("/extra2") else "/extra"
bc = "saturated" if p["satu"] else "linear"
save_df_path = hdd + f"/brunner/paper_models/ODE/{bc}/activation_q_10_R_1e4_Tsec_scan_5/"

p["f_Treg"] = 0.25 if Tregs == True else 0.0

repeats = 1
dr = ODEdriver(p)
fractions = np.linspace(0.5, 55, 60)*0.01
length = len(fractions)


systemic_R = p["f_Tsec"] * p["R_Tsec"] + (1 - p["f_Tsec"]) * p["R_Th"]

result = np.empty(length)
std = np.empty(length)
act = np.empty(length)

KD = p["K_D"]

rep_result = []
rep_std = []
rep_act = []
cell_list = []
for r,rep in enumerate(range(repeats)):
    print("repeat no. ", rep)
    result = np.empty(length)
    std = np.empty(length)
    act = np.empty(length)
    for f,fraction in enumerate(fractions):

        Th_R =  p["R_Th"]
        q = p["q"] #* p["f_Tsec"] / fraction
        if p["satu"] == False:
            c = dr.molecules_to_molar(
                I_lin(Th_R, q, fraction, p["f_Treg"], p["k_on"], p["R_Tsec"], p["R_Treg"], p["N_cells"],
                p["kd"])) * 1e12
        else:
            c = dr.molecules_to_molar(I2(Th_R, KD, q, fraction, p["f_Treg"], p["k_endo"], p["R_Tsec"], p["R_Treg"], p["N_cells"], p["kd"])) * 1e12

        for n in range(int(p["N_cells"] * (1 - p["f_Treg"] - fraction))):
            cell_list.append({"id":n,
                              "IL-2_Tsec_fraction": fraction,
                              "IL-2_R": Th_R,
                              "IL-2_surf_c": c * 1e-3,
                              "type_name": "Th",
                              "time": 1,
                              "replicat_index": r})

        for n in range(int(p["N_cells"] * fraction)):
            cell_list.append({"id":n + int(p["N_cells"] * (1 - p["f_Treg"] - fraction)),
                              "IL-2_Tsec_fraction": fraction,
                              "IL-2_R": p["R_Tsec"],
                              "IL-2_surf_c": c * 1e-3,
                              "type_name": "Tsec",
                              "time": 1,
                              "replicat_index": r})
        result[f] = np.mean(c)
        std[f] = c.std()
    rep_act.append(act)
    rep_result.append(result)
    rep_std.append(std)

mean_act = np.mean(rep_act, axis=0)
mean_result = np.mean(rep_result, axis=0)
mean_std = np.mean(rep_std, axis=0)

global_df = pd.DataFrame(columns=["Th_R", "IL-2_Tsec_fraction", "surf_c"])
global_df["IL-2_Tsec_fraction"] = fractions
global_df["surf_c"] = mean_result * 1e-3
global_df["surf_c_std"] = mean_std * 1e-3
global_df["Th_R"] = p["R_Th"]
global_df["type_name"] = "Th"
global_df["time"] = 1
global_df["replicat_index"] = 0
global_df["IL-2_sigma"] = p["sigma"]

cell_df = pd.DataFrame(cell_list)
if save_df == True:
    try:
        global_df.to_hdf(save_df_path + "global_df.h5", "w")
    except OSError:
        os.mkdir(save_df_path)
        global_df.to_hdf(save_df_path + "global_df.h5", "w")

    cell_df.to_hdf(save_df_path + "cell_df.h5", "w")