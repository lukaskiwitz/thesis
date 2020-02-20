import getpass

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from my_debug import message

user = getpass.getuser()
path = "/extra/{u}/scan_example_small/".format(u=user)

sns.set_context("paper",font_scale=1,rc={
            "lines.markersize":0,
            "lines.linewidth":2
            }
            )

global_df = pd.read_hdf(path+"global_df.h5", mode="r")
cell_df = pd.read_hdf(path+"cell_df.h5",mode="r")
stats_df = pd.read_hdf(path+"cell_stats_df.h5",mode="r")
#
#
#
means = stats_df.loc[:, (slice(None), "mean")]
means.columns = means.columns.droplevel(1)
means.reset_index(inplace=True)

std = stats_df.loc[:, (slice(None), "std")]
std.columns = std.columns.droplevel(1)
std.reset_index(inplace=True)

counts = stats_df.loc[:, (slice(None), "count")]
counts.columns = counts.columns.droplevel(1)
counts["n"] = counts["n"].apply(lambda x:int(x))
counts.reset_index(inplace=True)


# fig,ax = plt.subplots(2,2,sharex=True,sharey=False)
# sns.lineplot(x="t", y="surf_c_il2", data=means,hue="type_name",ax=ax[0][0],legend="full")
# sns.lineplot(x="t", y="surf_c_il6", data=means,hue="type_name",ax=ax[0][1],legend=False)
# sns.lineplot(x="t", y="surf_c_infg", data=means,hue="type_name",ax=ax[1][0],legend=False)
#
# plt.show()
# #
# fig,ax = plt.subplots(2,2,sharex=True,sharey=False)
# sns.lineplot(x="timeIndex", y="concentration", data=global_df,hue="field_name",ax=ax[0][0],legend=False,palette=sns.color_palette(["red","seagreen","navy"]))
# sns.lineplot(x="timeIndex", y="gradient", data=global_df,hue="field_name",ax=ax[0][1],palette=sns.color_palette(["red","seagreen","navy"]))
# sns.lineplot(x="timeIndex", y="sd", data=global_df,hue="field_name",ax=ax[1][0],palette=sns.color_palette(["red","seagreen","navy"]))
# sns.lineplot(x="t", y="n", data=counts,hue="type_name",ax=ax[1][1],legend="full")
# plt.show()
# # #
# #
# fig,ax = plt.subplots(2,2,sharex=True,sharey=False)
# sns.lineplot(x="t",y="Th1_score_norm",hue="type_name",data=cell_df,ax=ax[0][0])
# sns.lineplot(x="t",y="Tfh_score_norm",hue="type_name",data=cell_df,ax=ax[0][1])
# sns.lineplot(x="t",y="Tn_score_norm",hue="type_name",data=cell_df,ax=ax[1][0])
# plt.show()



fig,ax = plt.subplots(2,2,sharex=False,sharey=False)
sns.distplot(cell_df.loc[(cell_df["type_name"] == "Tn") & (cell_df["t"] == 1)]["surf_c_il2"],ax=ax[0][0])
sns.distplot(cell_df.loc[(cell_df["type_name"] == "Tn") & (cell_df["t"] == 1)]["surf_c_il6"],ax=ax[0][1])
sns.distplot(cell_df.loc[(cell_df["type_name"] == "Tn") & (cell_df["t"] == 1)]["surf_c_infg"],ax=ax[1][0])
plt.show()

# fig,ax = plt.subplots(2,2,sharex=False,sharey=False)
# sns.distplot(cell_df.loc[(cell_df["type_name"] == "Tn") & (cell_df["t"] == 50)]["surf_c_il2"],ax=ax[0][0])
# sns.distplot(cell_df.loc[(cell_df["type_name"] == "Tn") & (cell_df["t"] == 50)]["surf_c_il6"],ax=ax[0][1])
# sns.distplot(cell_df.loc[(cell_df["type_name"] == "Tn") & (cell_df["t"] == 50)]["surf_c_infg"],ax=ax[1][0])
# plt.show()


tfh = cell_df.loc[
    (cell_df["surf_c_il6"] > 0.047)&
    (cell_df["surf_c_il2"] < 0.065)&
    (cell_df["type_name"] == "Tn")
            ].groupby(["t","scan_index","type_name"]).count()["n"].reset_index()
th1 = cell_df.loc[
    (cell_df["surf_c_infg"] > 0.055)&
    (cell_df["type_name"] == "Tn")
            ].groupby(["t","scan_index","type_name"]).count()["n"].reset_index()

fig,ax = plt.subplots(2,2,sharex=True,sharey=False)
sns.lineplot(x="t",y="n",data=tfh,ax=ax[0][0])
sns.lineplot(x="t",y="n",data=th1,ax=ax[0][1])
plt.show()
