import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

get_dataframes = []#[[]]*3
saving_dataframe = pd.DataFrame()

from parameters_q_fraction import path
# get_dataframes.append([])
# # path = "/extra/brunner/10x10x10/R_lognorm/run" + str(j) + "/"
# # path = "/extra/brunner/10x10x10/q_fraction_exact/run" + str(j) + "/"
# user = getpass.getuser()
path = "/extra/brunner/thesis/static/q_fraction_medium/"
# # ext_cache="/extra/brunner/para_handling/kinetic/R_lognorm_ext_cache/"

T = np.arange(0, 200, 1)
dt = 3600

global_df = pd.read_hdf(path + 'global_df.h5', mode="r")
cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")
# print(global_df["Concentration"]*1e3)


# subplot_style
a_x = 3
a_y = 2
x_size = 4
y_size = 4

# sns.set_context("talk")
# sns.set_context("talk", font_scale=1,rc={"lines.linewidth": 1})
sns.set(rc={'figure.figsize':(11,8.27)})
fig = plt.figure()
plt.subplots_adjust(wspace=.3)

myhue = None
fig.add_subplot(a_x,a_y, 1)
#sns.lineplot(x="fraction", y="mean_surf_c_il2", data=saving_dataframe,hue=myhue)
ax1  = sns.lineplot(x="IL-2_fraction", y="surf_c", data=global_df,hue=myhue)
#ax1.errorbar(np.linspace(0.00005,1.0,20), global_df["surf_c"], yerr=global_df["sd"], fmt='-o')
ax1.set(xlabel="IL-2_fraction", ylabel="mean surface c. in nM")#, ylim=(0.14125,0.142))
# fig.savefig("plots/q_fraction_1.png", bbox_inches='tight')


# fig = plt.figure()
fig.add_subplot(a_x,a_y, 2)
ax3 = sns.lineplot(x="IL-2_fraction", y="Concentration", data=global_df,hue=myhue, legend=False)
ax3.set(xlabel="IL-2_fraction", ylabel="concentration in nM")#, ylim=(0.159,0.164))
# fig.savefig("plots/q_fraction_2.png", bbox_inches='tight')


# fig = plt.figure()
fig.add_subplot(a_x,a_y, 3)
global_df["Gradient"] = global_df["Gradient"]*1e5
ax2 = sns.lineplot(x="IL-2_fraction", y="Gradient", data=global_df,hue=myhue, legend=False)#
ax2.set(xlabel="IL-2_fraction", ylabel="mean gradient in nM/um")
# fig.savefig("plots/q_fraction_3.png", bbox_inches='tight')

# fig = plt.figure()
fig.add_subplot(a_x,a_y, 4)
global_df["surf_c_std_norm"] = global_df["surf_c_std"]/global_df["surf_c"]
ax4 = sns.lineplot(x="IL-2_fraction", y="surf_c_std_norm", data=global_df,hue=myhue, legend=False)
ax4.set(xlabel="IL-2_fraction", ylabel="surface c. std/mean ")
# fig.savefig("plots/q_fraction_4.png", bbox_inches='tight')

# fig = plt.figure()
fig.add_subplot(a_x,a_y, 5)
global_df["std_norm"] = global_df["SD"]*1e9/global_df["Concentration"]
ax4 = sns.lineplot(x="IL-2_fraction", y="std_norm", data=global_df,hue=myhue, legend=False)
ax4.set(xlabel="IL-2_fraction", ylabel="std/mean")

plt.show()