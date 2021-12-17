import getpass
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os
from thesis.main.my_debug import message

def EC50_calculation(E_max, E_min, k, N, R):
    return (E_max * k ** N + E_min * R ** N) / (k ** N + R ** N)

save_plot = True

hdd = "extra2" if os.path.exists("/extra2") else "extra"
user = getpass.getuser()
model_name = "clustering"
subname = "test"
saving_string = "/{extra}/{u}/paper_models/kinetics/{mn}/{sn}/".format(u=user, mn=model_name, sn = subname, extra = hdd)

yscale = "linear"
xscale = "linear"

# define data start and end point
startingPoint = None
stoppingPoint = None
# plotting limits
xlim = (None, None)
ylim = (None, None)
# time and fraction max limit
max_time = 100
max_frac = -1
# colours, first = positive, second = negative


sns.set(rc={'figure.figsize': (7, 5.5)})
sns.set_style("ticks")
sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 2.5})
fig = plt.figure()

pos_path = "/{extra}/{u}/paper_models/kinetics/{mn}/{sn}/".format(u=user, sn=subname, mn=model_name, extra = hdd)
# neg_path = "/{extra}/{u}/paper_models/kinetics/{mn}/{sn}/".format(u=user, sn=subname, mn=model_name, extra = hdd)

pos_activation_df = pd.read_hdf(pos_path + 'activation_df_0.h5', mode="r")
amount_of_replicats = 1
# neg_activation_df = load_runs(neg_path, 1)

# pos_activation_df = pos_activation_df[pos_activation_df.index < max_time]
# neg_activation_df = neg_activation_df[neg_activation_df.index < max_time]

pos_activation_df["time"] /= 3600

for d,df in enumerate([pos_activation_df]):
    print(df)
    ste = []
    mean = []
    # sliced_df = df.loc[df.index ==  df.index.max()]
    for s, scan in enumerate(df["scan_index"].unique()):
        #     ste.append(df[scan].std() / np.sqrt(amount_of_replicats))
        #     mean.append(df[scan].mean())
        # plt.fill_between(np.sort(df.columns[1:max_frac].astype("float")), np.clip(np.array(mean) - np.array(ste), 0, 1),
        #                  np.clip(np.array(mean) + np.array(ste), 0, 1), alpha=0.15, color=colours[d])

        ax_1 = sns.scatterplot(x=np.sort(df.time.unique().astype("float")), y=df.loc[df["scan_index"] == scan, "activation"].values, s=30)
        ax_2 = sns.lineplot(x=np.sort(df.time.unique().astype("float")), y=df.loc[df["scan_index"] == scan, "activation"].values, linewidth = "2.5", label=scan)


# import matplotlib.lines as mlines
#
# white_line = mlines.Line2D([], [], color='white', alpha=1)
# red_line = mlines.Line2D([], [], color='red')
# blue_line = mlines.Line2D([], [], color='blue')
# # red_line = mlines.Line2D([], [], color='red')
#
# plt.legend([white_line, red_line, blue_line],
#            ["Spatial", "positive", "negative"], loc='best',
#            prop={'size': 17})
handles, labels = ax_2.get_legend_handles_labels()
# labels = ["no clustering", "medium", "strong"]
plt.legend(handles, labels)

plt.ylabel("fraction of activated cells")
plt.xlabel("time (h)")
plt.yscale(yscale)
plt.xscale(xscale)
plt.ylim(ylim)
plt.xlim(xlim)

plt.tight_layout()
if save_plot == True:
    fig.savefig(saving_string + "test_activation.svg", bbox_inches='tight')
plt.show()