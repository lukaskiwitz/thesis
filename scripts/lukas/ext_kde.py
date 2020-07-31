import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import getpass
import KDEpy
import os
import multiprocessing as mp
from skimage.measure import label

# user = getpass.getuser()
# model_name = "example_min"
# name = "test"
# path = "/extra/{u}/{mn}/{n}/".format(u=user, n=name, mn=model_name)
# IMGPATH = path
# os.makedirs(IMGPATH,exist_ok = True)
#
# keep_columns = ["id","x","y","z","type_name","IL-2_surf_c","time_index"]
#
# cell_df: pd.DataFrame = pd.read_hdf(path + "cell_df.h5", mode="r")
# cell_df = cell_df.drop(columns=cell_df.columns.drop(keep_columns))

def evalutate_kernel_on_grid(kernel,grid_points):


    grid, points = kernel.evaluate(grid_points)
    x, y = np.unique(grid[:, 0]), np.unique(grid[:, 1])
    v = points.reshape(grid_points, grid_points).T

    return x,y,v


def get_kde_estimators(ts, time_index, type_names,kernel_type, bw):
    kernels = {}
    vis_kernels = {}

    fig,ax = plt.subplots(
        3,3,figsize = (10,10)
    )
    ax = np.ravel(ax)

    for i,type_name in enumerate(type_names):
        inital_cells = ts.loc[(ts["time_index"] == time_index) & (ts["type_name"] == type_name)]
        if inital_cells.shape[0] == 0:
            break
        elif inital_cells.shape[0] == 1:
            data = np.array([inital_cells["x"].iloc[0], inital_cells["y"].iloc[0], inital_cells["z"].iloc[0]])

        else:
            data = np.array([inital_cells["x"], inital_cells["y"], inital_cells["z"]]).T

        kernel = KDEpy.TreeKDE(kernel_type, bw=bw).fit(data)
        kernels[type_name] = kernel

        data = data.T[0:2].T
        kernel = KDEpy.FFTKDE(kernel_type,bw=bw).fit(data)
        vis_kernels[type_name] = kernel

        x,y,v = evalutate_kernel_on_grid(kernel,200)

        ax[i].set_title(type_name)
        ax[i].contourf(x,y,v,100)
        sns.scatterplot(x="x",y="y",data=inital_cells,s=0.1,ax = ax[i])

    return kernels, vis_kernels, ax

def get_t_df (df,t):
    return df.loc[df["time_index"] == t]

def compute_kde(cell_df,kernel_type = "tri", bw = 20,lag = 1, dim = 3):

    result = []
    axes_list =  []
    vis_kernel_list = []
    for time_index,df in cell_df.groupby("time_index"):

        if lag == "init":
            lag = time_index

        t = time_index-lag if time_index >= lag else 0

        kernels,vis_kernels, axes = get_kde_estimators(cell_df,t-lag,cell_df["type_name"].unique(),kernel_type,bw)
        axes_list.append(axes)
        vis_kernel_list.append(vis_kernels)


        for type_name, kernel in kernels.items():
            score_name = type_name + "_score"
            values = kernel(np.array([df.loc[:, "x"], df.loc[:, "y"], df.loc[:, "z"]]).T)

            values = pd.Series(values)
            values.index = df.index
            values = values.div(values.mean())
            df[score_name] = values


        result.append(df)



        sns.boxplot(x="type_name", y="Th1_score", data=df,ax = axes[3])
        sns.boxplot(x="type_name", y="Th_score", data=df, ax=axes[4])
        sns.boxplot(x="type_name", y="Tfh_score", data=df, ax=axes[5])


    return pd.concat(result),axes_list, vis_kernel_list

def get_connected_regions(kernel_list ,ax = None):
    result_outer = []
    for k_dict in kernel_list:
        result = {}
        for i,(type_name, kernel) in enumerate(k_dict.items()):

            x,y,v = evalutate_kernel_on_grid(kernel,1000)

            m = np.mean(v)
            v[v < m] = 0
            v[v >= m] = 1

            labels = label(v,connectivity=2)
            # plt.colorbar(
            #     plt.contourf(x,y,v,10)
            # )

            # plt.show()

            if not ax is None:

                ax[i].set_title(type_name)

                ax[i].contourf(x,y,labels,100, cmap = "Dark2")

            result[type_name] = len(np.unique(labels))
        result_outer.append(result)
    return result_outer


def init(df):

    global cell_df_g
    cell_df_g = df

def target(tupel):

    index = tupel[0]
    bw = tupel[1]

    cell_df,axes_list,vis_kernel_list = compute_kde(cell_df_g,kernel_type="tri", bw=bw,lag = "init")
    bw_series = pd.Series([bw]*len(cell_df))
    bw_series.index = cell_df.index
    cell_df["bw"] = bw_series

    labels = get_connected_regions(vis_kernel_list,ax = axes_list[0][-3:])[0]

    for type_name in cell_df["type_name"].unique():
        n = len(cell_df.loc[cell_df["type_name"] == type_name])
        labels[type_name] = labels[type_name]/n

    labels_df = pd.DataFrame([labels])



    bw_series = pd.Series([bw] * len(labels_df))
    bw_series.index = labels_df.index
    labels_df["bw"] = bw_series

    plt.tight_layout()
    img = IMGPATH + "images/"
    os.makedirs(img, exist_ok=True)
    plt.gcf().suptitle(r'bw $= {bw}$'.format(bw=bw))

    plt.tight_layout(pad = 2, rect=[0, 0.03, 1, 0.95])
    plt.savefig(img + "{i}.png".format(i = index), dpi=600)

    return [cell_df, labels_df]


def time_plot(type_name,data, ax_l,):
    data = data.loc[data["type_name"] == type_name]
    sns.lineplot(x="time_index", y=type_name+"_score", data=data, hue="bw", ci=None, ax=ax_l)

def bw_plot(type_name,data,ax_l):
    # data = data.loc[data["type_name"] == type_name]
    sns.lineplot(x="bw", y=type_name+"_score", data=data, hue="type_name", ax=ax_l)

import glob

user = getpass.getuser()
path = "/home/{u}/hauser_ext_kde/".format(u=user)
IMGPATH = path
os.makedirs(IMGPATH,exist_ok = True)

files = []
for file_path in glob.glob(path+"*.csv"):
    with open(file_path) as f:

        df = pd.read_csv(f,names=["x","y"])

        type_name = file_path.split("/")[-1].split(".")[0].split("_")[0]

        series = pd.Series([type_name]*len(df))
        series.index = df.index

        df["type_name"] = series

        series = pd.Series([0] * len(df))
        series.index = df.index

        df["time_index"] = series

        series = pd.Series([0] * len(df))
        series.index = df.index

        df["z"] = series

        files.append(df)


cells = pd.concat(files)
steps = 1000
n = 80

with mp.Pool(n,initializer=init,initargs=[cells]) as p:
    result = p.map(target,enumerate(np.linspace(1,50,steps)))

result = np.transpose(result)
labels = pd.concat(result[1])
result = pd.concat(result[0])

plt.figure()
sns.lineplot(x = "bw", y = "Th1", data = labels)
sns.lineplot(x = "bw", y = "Tfh", data = labels)
sns.lineplot(x = "bw", y = "Th", data = labels)
plt.legend(["Th1","Tfh","Th"])
plt.ylabel("number of connected regions")
plt.savefig(IMGPATH+"connected.pdf")
plt.xscale("log")
plt.savefig(IMGPATH+"connected_log.pdf")
plt.show()

fig,ax = plt.subplots(3,2,figsize = (10,10))
ax = np.ravel(ax)

ax_l = ax[0]
bw_plot("Th1",result,ax_l)
# ax_l = ax[1]
# time_plot("abs",result,ax_l)

ax_l = ax[2]
bw_plot("Tfh",result,ax_l)
# ax_l = ax[3]
# time_plot("sec",result,ax_l)

ax_l = ax[4]
bw_plot("Th",result,ax_l)
# ax_l = ax[5]
# time_plot("default",result,ax_l)


plt.savefig(IMGPATH+"ext_kde.pdf")
plt.show()

# for i,data in cells.groupby(["type_name"]):
#     sns.scatterplot(x = "x", y = "y", data = data)
#     plt.show()


# # r,l,kernel_list = compute_kde(cells, bw = 50, lag=0)
# global cell_df_g
# cell_df_g = cells
# target(10)
# plt.figure()

