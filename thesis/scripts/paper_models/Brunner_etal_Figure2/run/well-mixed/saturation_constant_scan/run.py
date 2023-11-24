import numpy as np
import pandas as pd
import os
from parameters import p
from driver import ODEdriver

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

p["satu"] = True
Tregs = False


save_df = True
hdd = "/extra2" if os.path.exists("/extra2") else "/extra"
save_df_path = hdd + "/brunner/paper_models/ODE/saturated/KD_scan_fixed/"


p["f_Treg"] = 0.25 if Tregs == True else 0.0

repeats = 20

dr = ODEdriver(p)

systemic_R = p["f_Tsec"] * p["R_Tsec"] + (1 - p["f_Tsec"]) * p["R_Th"]
Th_R = (systemic_R - p["f_Tsec"] * p["R_Tsec"]) / (1 - p["f_Tsec"])

KDs = np.logspace(-1, 1, 31) * 7.437
KDs = dr.molar_to_molecules(np.sort(np.concatenate([KDs, [7.437 * 1/5, 7.437 * 5]])) * 1e-12)

length = len(KDs)

rep_result = []
rep_std = []
cell_list = []
for r,rep in enumerate(range(repeats)):
    print("repeat no. ", rep)
    result = np.empty(length)
    std = np.empty(length)
    act = np.empty(length)
    for i, KD in enumerate(KDs):
        varied_Rs = dr.set_R_lognorm(Th_R, p["sigma"] * Th_R, length=int(p["N_cells"] * (1 - p["f_Treg"] - p["f_Tsec"])))
        if p["satu"] == False:
            c = dr.molecules_to_molar(
                I_lin(varied_Rs.mean(), p["q"], p["f_Tsec"], p["f_Treg"], p["k_on"], p["R_Tsec"], p["R_Treg"], p["N_cells"],
                p["kd"])) * 1e12
        else:
            c = dr.molecules_to_molar(I2(varied_Rs.mean(), KD, p["q"], p["f_Tsec"], p["f_Treg"], p["k_endo"], p["R_Tsec"],
                                         p["R_Treg"], p["N_cells"], p["kd"])) * 1e12

        for n in range(int(p["N_cells"] * (1 - p["f_Treg"] - p["f_Tsec"]))):
            cell_list.append({"id":n,
                              "IL-2_Tsec_fraction": p["f_Tsec"],
                              "IL-2_sigma": p["sigma"],
                              "IL-2_KD": dr.molecules_to_molar(KD) * 1e12,
                              "IL-2_R": varied_Rs[n],
                              "IL-2_surf_c": c * 1e-3,
                              "type_name": "Th",
                              "time": 1,
                              "replicat_index": r,
                              "scan_name": "KD"})

        for n in range(int(p["N_cells"] * p["f_Tsec"])):
            cell_list.append({"id":n + int(p["N_cells"] * (1 - p["f_Treg"])),
                              "IL-2_Tsec_fraction": p["f_Tsec"],
                              "IL-2_sigma": p["sigma"],
                              "IL-2_KD": dr.molecules_to_molar(KD) * 1e12,
                              "IL-2_R": p["R_Tsec"],
                              "IL-2_surf_c": c * 1e-3,
                              "type_name": "Tsec",
                              "time": 1,
                              "replicat_index": r,
                              "scan_name": "KD"})

        result[i] = np.mean(c)
        std[i] = c.std()
    rep_result.append(result)
    rep_std.append(std)

mean_result = np.mean(rep_result, axis=0)
mean_std = np.mean(rep_std, axis=0)

global_df = pd.DataFrame(columns=["Th_R", "IL-2_KD", "surf_c", "surf_c_std"])
global_df["IL-2_KD"] = dr.molecules_to_molar(KDs) * 1e12
global_df["IL-2_sigma"] = p["sigma"]
global_df["surf_c"] = mean_result * 1e-3
global_df["surf_c_std"] = mean_std * 1e-3
global_df["Th_R"] = 5000
global_df["type_name"] = "Th"
global_df["time"] = 1
global_df["replicat_index"] = 0
global_df["IL-2_Tsec_fraction"] = p["f_Tsec"]
global_df["scan_name"] = "KD"

cell_df = pd.DataFrame(cell_list)
if save_df == True:
    try:
        global_df.to_hdf(save_df_path + "global_df.h5", "w")
    except OSError:
        os.mkdir(save_df_path)
        global_df.to_hdf(save_df_path + "global_df.h5", "w")

    cell_df.to_hdf(save_df_path + "cell_df.h5", "w")