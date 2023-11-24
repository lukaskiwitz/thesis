def get_path(bc, hdd, model, scan_variable):
    if scan_variable == "IL-2_Tsec_fraction":
        if model == "loc. q and R":
            path = f"/extra2/brunner/paper_models/statics/saturated/Figure_2C/dataframes_Tsec_0.05_timeseries/" if bc == "saturated" else f"/extra2/brunner/paper_models/statics/linear/Tsec_scan_4/"
        if model == "well_mixed":
            path = "/extra2/brunner/paper_models/ODE/saturated/Tsec_scan_fixed/" if bc == "saturated" else "/extra2/brunner/paper_models/ODE/{bc}/Tsec_scan_2/".format(bc=bc)

    elif scan_variable == "IL-2_sigma":
        if model == "loc. q and R":
            path = f"/extra2/brunner/paper_models/statics/saturated/Figure_2C/dataframes_R_Tsec_0.05_steady_state/" if bc == "saturated" else "/extra2/brunner/paper_models/statics/linear/R_hetero_log10/"
        if model == "well_mixed":
            path = f"/extra2/brunner/paper_models/ODE/saturated/R_hetero_0.05/" if bc == "saturated" else "/extra2/brunner/paper_models/ODE/linear/R_hetero_log10/"

    elif scan_variable == "IL-2_KD":
        if model == "loc. q and R":
            path = "/extra2/brunner/paper_models/statics/saturated/Figure_2C/dataframes_KD_Tsec_0.05_steady_state/"
        if model == "well_mixed":
            path = "/extra2/brunner/paper_models/ODE/saturated/KD_scan_fixed/"
        if model == "k_on_linear":
            path = "/extra2/brunner/paper_models/statics/linear/k_on_scan/"
    return path

def scale_dataframes(c_df, g_df, model, scan_variable, sv, standards):
    g_df["surf_c_std_norm"] = g_df["surf_c_std"] / g_df["surf_c"]
    g_df["surf_c"] *= 1e3
    g_df["surf_c_std"] *= 1e3
    c_df["IL-2_surf_c"] *= 1e3
    if scan_variable == "IL-2_KD":
        if model != "well_mixed":
            try:
                if g_df.scan_name_scan_name.unique()[0][:2] == "KD":
                    g_df["IL-2_KD"] = g_df.scan_value * 1e3
                    c_df["IL-2_KD"] = c_df.scan_value * 1e3
            except AttributeError:
                if g_df.scan_name.unique()[0][:2] == "KD":
                    g_df["IL-2_KD"] = g_df.scan_value * 1e3
                    c_df["IL-2_KD"] = c_df.scan_value * 1e3
        if model == "loc. q" or model == "well_mixed":
            from thesis.scripts.patrick.ODE.driver import ODEdriver
            from thesis.scripts.patrick.ODE.parameters import p
            try:
                p["N_cells"] = len(c_df["id"].unique())
            except KeyError:
                p["N_cells"] = len(c_df["id_id"].unique())
            ODEdr = ODEdriver(p)

        c_df["IL-2_KD"] = 1 * c_df["IL-2_KD"]
        g_df["IL-2_KD"] = 1 * g_df["IL-2_KD"]

    g_df["fold_change"] = g_df[scan_variable] / standards[sv]
    c_df["fold_change"] = c_df[scan_variable] / standards[sv]

    g_df["scale"] = g_df["fold_change"] / (g_df["fold_change"].max())
    try:
        c_df["id"]
    except KeyError:
        c_df["id"] = c_df["id_id"]
    return c_df, g_df