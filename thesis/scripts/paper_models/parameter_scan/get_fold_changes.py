import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from parameters import path, IMGPATH
from thesis.main.MyPlotter import Plotter
from thesis.main.my_debug import message

path = path
plotter = Plotter(path)
message("dataset loaded")


def get_fold_change(df, aggregator, group_by_columns):
    def f(df):

        scan_name = df[plotter.scan_name_key].iloc[0]
        # df = df.groupby("scan_value").apply(aggregator)
        # df = df.drop(columns=["scan_value"]).reset_index()

        low = np.min(df.scan_value)
        high = np.max(df.scan_value)
        ds = df.loc[(df.scan_value == 1)]

        df = df.loc[(df.scan_value == low) | (df.scan_value == high) ]
        props = ["Concentration","SD","CV","surf_c","surf_c_std","surf_c_cv","MeshVolume"]
        for p in props:
            df[p] = df[p] - ds[p].mean()

        df["fold_change"] = high
        df["explode"] = np.abs(df["Concentration"].mean()) > 1

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



# global_std = get_fold_change(plotter.global_df, np.std, ["numeric_linear", plotter.scan_name_key])\
#              - get_standard_change(plotter.global_df, np.std, ["numeric_linear", plotter.scan_name_key])
#

# """mean and standard deviation over cells pooled from all replicats"""
# cell_aggregrate_mean = get_fold_change(plotter.cell_df, np.mean, ["numeric_linear", plotter.scan_name_key, "type_name"])
# cell_aggregate_std = get_fold_change(plotter.cell_df, np.std, ["numeric_linear", plotter.scan_name_key, "type_name"])
#
# """
# mean and standard error for replicats
# """
# cell_temp = get_fold_change(plotter.cell_df, np.mean,
#                             ["numeric_linear", plotter.scan_name_key, "type_name", plotter.time_index_key])
# cell_sem = cell_temp.groupby(["numeric_linear", plotter.scan_name_key, "type_name", "scan_value"]).std().reset_index()
# cell_replicats_mean = cell_temp.groupby(
#     ["numeric_linear", plotter.scan_name_key, "type_name", "scan_value"]).mean().reset_index()


df = global_mean.loc[
    # (global_mean["explode"] == False) &
    (global_mean["numeric_linear"] == False)
    ]


a = 8.3/6
b = np.sqrt(2) * a

sns.catplot(x = plotter.scan_name_key,
            y = "CV",
            data=df,
            kind = "bar",
            hue = "scan_value",
            height = b,
            legend=None,
            aspect=2**0.5,
            linewidth=1, facecolor=(1, 1, 1, 0),
            errcolor=".2", edgecolor=".1",
            sharey=True,
            sharex=False,
            dodge = True,
            order = ["D","ratio","abs_R","sec_q","distance","Koff","kendo"]
            )


# plt.ylim([-10,10])
plt.tight_layout()
plt.savefig(IMGPATH+"fold_changes.pdf")
plt.show()
