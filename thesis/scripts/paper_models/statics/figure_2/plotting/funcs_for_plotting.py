def get_path(bc, hdd, model, scan_variable):
    if scan_variable == "IL-2_Tsec_fraction":
        if model == "loc. q and R":
            # path = hdd + "/brunner/paper_models/statics/standard/Tsec_scan/"
            path = hdd + "/brunner/paper_models/statics/{bc}/Tsec_scan_0.05_standard_3/".format(bc=bc)
        elif model == "loc. q":
            # path = hdd + "/brunner/paper_models/yukawa/standard/Tsec_scan/"
            path = hdd + "/brunner/paper_models/yukawa/{bc}/Tsec_scan_low_R_log2_1/".format(bc=bc)
        if model == "well_mixed":
            # path = hdd + "/brunner/paper_models/ODE/standard/Tsec_scan/"
            path = hdd + "/brunner/paper_models/ODE/{bc}/Tsec_scan_0.05_standard_2/".format(bc=bc)

    elif scan_variable == "IL-2_sigma":
        if model == "loc. q and R":
            # path = hdd + "/brunner/paper_models/statics/standard/R_lognorm/"
            path = hdd + "/brunner/paper_models/statics/{bc}/R_lognorm_0.05_standard_3/".format(bc=bc)
        elif model == "loc. q":
            # path = hdd + "/brunner/paper_models/yukawa/standard/R_lognorm/"
            path = hdd + "/brunner/paper_models/yukawa/{bc}/R_lognorm_low_R_log2_1/".format(bc=bc)
        if model == "well_mixed":
            # path = hdd + "/brunner/paper_models/ODE/standard/R_lognorm/"
            path = hdd + "/brunner/paper_models/ODE/{bc}/R_lognorm_0.05_standard_2/".format(bc=bc)

    elif scan_variable == "IL-2_KD":
        if model == "loc. q and R":
            path = hdd + "/brunner/paper_models/statics/saturated/KD_scan_0.05_standard_3/"
        elif model == "loc. q":
            path = hdd + "/brunner/paper_models/yukawa/saturated/KD_scan_log2_5/"
        if model == "well_mixed":
            path = hdd + "/brunner/paper_models/ODE/saturated/KD_scan_0.05_standard_2/"
    return path

def scale_dataframes(c_df, g_df, model, scan_variable, sv, min_value, max_value, standards):
    g_df["surf_c_std_norm"] = g_df["surf_c_std"] / g_df["surf_c"]
    g_df["surf_c"] *= 1e3
    g_df["surf_c_std"] *= 1e3
    c_df["IL-2_surf_c"] *= 1e3
    if scan_variable == "IL-2_KD":
        if model != "well_mixed":
            try:
                if g_df.scan_name_scan_name.unique()[0] == "KD":
                    g_df["IL-2_KD"] = g_df.scan_value * 1e3
                    c_df["IL-2_KD"] = c_df.scan_value * 1e3
            except AttributeError:
                if g_df.scan_name.unique()[0] == "KD":
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
            g_df["IL-2_KD"] = ODEdr.molecules_to_molar(g_df["IL-2_KD"].values) * 1e12
            c_df["IL-2_KD"] = ODEdr.molecules_to_molar(c_df["IL-2_KD"].values) * 1e12

        c_df["IL-2_KD"] = 1 / c_df["IL-2_KD"]
        g_df["IL-2_KD"] = 1 / g_df["IL-2_KD"]

    # if model == "loc. q":
    #     g_df = g_df.loc[g_df["scan_index"] < g_df["scan_index"].max()]
    #     c_df = c_df.loc[c_df["scan_index"] < c_df["scan_index"].max()]

    c_df = c_df.loc[c_df[scan_variable] <= max_value * standards[sv]]
    g_df = g_df.loc[g_df[scan_variable] <= max_value * standards[sv]]

    c_df = c_df.loc[c_df[scan_variable] >= min_value * standards[sv]]
    g_df = g_df.loc[g_df[scan_variable] >= min_value * standards[sv]]

    g_df["fold_change"] = g_df[scan_variable] / standards[sv]
    c_df["fold_change"] = c_df[scan_variable] / standards[sv]

    g_df["scale"] = g_df["fold_change"] / (g_df["fold_change"].max())
    try:
        c_df["id"]
    except KeyError:
        c_df["id"] = c_df["id_id"]
    return c_df, g_df