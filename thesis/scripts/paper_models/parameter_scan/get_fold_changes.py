from thesis.main.MyPlotter import Plotter
import numpy as np
from thesis.main.my_debug import message
from parameters import path

path = path
plotter = Plotter(path)
message("dataset loaded")


def get_fold_change(df, aggregator, group_by_columns):
    def f(df):

        df = df.groupby("scan_value").apply(aggregator)
        df = df.drop(columns=["scan_value"]).reset_index()

        low = 0.1
        high = 10

        df = df.loc[(df.scan_value == low) | (df.scan_value == high)]

        return df

    gb = df.groupby(group_by_columns)
    df = gb.apply(f)
    for k in group_by_columns:
        if k in df.columns:
            df = df.drop(columns=[k])
    df = df.reset_index()
    return df


""" here std could be intepreted as the sem of average cytokine concentration"""
global_mean = get_fold_change(plotter.global_df, np.mean, ["numeric_linear", plotter.scan_name_key])
global_std = get_fold_change(plotter.global_df, np.std, ["numeric_linear", plotter.scan_name_key])

"""mean and standard deviation over cells pooled from all replicats"""
cell_aggregrate_mean = get_fold_change(plotter.cell_df, np.mean, ["numeric_linear", plotter.scan_name_key, "type_name"])
cell_aggregate_std = get_fold_change(plotter.cell_df, np.std, ["numeric_linear", plotter.scan_name_key, "type_name"])

"""
mean and standard error for replicats
"""
cell_temp = get_fold_change(plotter.cell_df, np.mean,
                            ["numeric_linear", plotter.scan_name_key, "type_name", plotter.time_index_key])
cell_sem = cell_temp.groupby(["numeric_linear", plotter.scan_name_key, "type_name", "scan_value"]).std().reset_index()
cell_replicats_mean = cell_temp.groupby(
    ["numeric_linear", plotter.scan_name_key, "type_name", "scan_value"]).mean().reset_index()
