import logging

import KDEpy
import numpy as np
import pandas as pd

module_logger = logging.getLogger(__name__)


def evalutate_kernel_on_grid(kernel, grid_points):
    grid, points = kernel.evaluate(grid_points)
    x, y = np.unique(grid[:, 0]), np.unique(grid[:, 1])
    v = points.reshape(grid_points, grid_points).T

    return x, y, v


def evaluate_kernel_at_points(kernel, points):
    pass


def get_kde_from_df(df, kernel_type, bw, visual=True):
    data = []

    if df.shape[0] == 0:
        pass
    elif df.shape[0] == 1:
        data = np.expand_dims(np.array([df["x"].iloc[0], df["y"].iloc[0], df["z"].iloc[0]]), axis=1).T
    else:
        data = np.array([df["x"], df["y"], df["z"]]).T

    kernel = KDEpy.TreeKDE(kernel_type, bw=bw).fit(data)

    return kernel


def get_cell_df(cell_list):
    df = pd.DataFrame()

    for cell in cell_list:
        d = cell.p.get_as_dictionary()
        d["x"] = cell.center[0]
        d["y"] = cell.center[1]
        d["z"] = cell.center[2]

        df = df.append(d, ignore_index=True)
    return df

# def compute_kde(cell_df,kernel_type = "tri", bw = 20,lag = 0):
#
#     result = []
#     vis_kernel_list = []
#     tree_kernel_list = []
#
#     for time_index,df in cell_df.groupby("time_index"):
#
#         if lag == "init":
#             lag = time_index
#
#         t = time_index-lag if time_index >= lag else 0
#
#         tree_kernels,vis_kernels = get_kde_estimators(cell_df,t-lag,cell_df["type_name"].unique(),kernel_type,bw)
#         # axes_list.append(axes)
#         vis_kernel_list.append(vis_kernels)
#         tree_kernel_list.append(tree_kernels)
#
#
#         for type_name, kernel in tree_kernels.items():
#             score_name = type_name + "_score"
#             values = kernel(np.array([df.loc[:, "x"], df.loc[:, "y"], df.loc[:, "z"]]).T)
#
#             values = pd.Series(values)
#             values.index = df.index
#             values = values.div(values.mean())
#             df[score_name] = values
#
#
#         result.append(df)
#         result = pd.concat(result)
#
#     return result,tree_kernel_list, vis_kernel_list
