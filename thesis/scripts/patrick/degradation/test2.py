import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from scipy.stats import lognorm
from scipy.constants import N_A

sns.set(rc={'figure.figsize':(7,7)})
sns.set_style("ticks")
sns.set_context("talk", font_scale=1, rc={"lines.linewidth": 5})
fig = plt.figure()

p = {
    "k_on": 15 * 10**7,  #1e9 *  111.6 / 60 ** 2,  # 111.6 per hour # now 540 per hour per nM
    "rho": 0.05,  # mu
    "D": 0.001,  # mu² per s
    "R_h": 4000 * N_A ** -1 * 1e9,
    "R_l": 20000 * N_A ** -1 * 1e9,
    "kd": 0,#0.1/(60*2),
    "q_h": 10 * N_A ** -1 * 1e9,
    "q_l": 60 * N_A ** -1 * 1e9,
    "fraction": 1,
    "sigma" : 7764,
    "gamma": 10
}

def func(R):
	if R < 0:
		R = 0
	c = p["q_l"]/p["k_on"] * (4*np.pi * p["D"] * 0.01 * 0.02 + p["k_on"] * 0.01 * 999 * p["R_l"])/ \
	(p["k_on"] * 0.01 * 999 * p["R_l"] * R + 4*np.pi * p["D"] * 0.01 * 0.02*(R + 999*p["R_l"]))
	return c
sigmas = np.linspace(1,20000,50) * N_A ** -1 * 1e9
results = []
for sigma in sigmas:
	temp = np.zeros(10)
	for i in range(len(temp)):
		func_values = np.zeros(1000)
		for ii in range(len(func_values)):
			draw = np.random.normal(p["R_l"],sigma)
			func_values[ii] = func(draw)
		temp[i] = np.mean(func_values)
	results.append(np.mean(temp)*1e9)
ax = sns.lineplot(sigmas/p["R_l"], results, label="normal")
ax.set(xlabel = "sigma/E", ylabel = "concentration in pM")
# print(results)
# plt.show()

p = {
    "k_on": 15 * 10**7,  #1e9 *  111.6 / 60 ** 2,  # 111.6 per hour # now 540 per hour per nM
    "rho": 0.05,  # mu
    "D": 0.001,  # mu² per s
    "R_h": 4000 * N_A ** -1 * 1e9,
    "R_l": 20000 * N_A ** -1 * 1e9,
    "kd": 0,#0.1/(60*2),
    "q_h": 10 * N_A ** -1 * 1e9,
    "q_l": 60 * N_A ** -1 * 1e9,
    "fraction": 1,
    "sigma" : 7764,
    "gamma": 10
}

#high cell density
def func(R):
    c = p["q_l"]/p["k_on"] * (4*np.pi * p["D"] * 0.01 * 0.02 + p["k_on"] * 0.01 * 999 * p["R_l"])/ \
        (p["k_on"] * 0.01 * 999 * p["R_l"] * R + 4*np.pi * p["D"] * 0.01 * 0.02*(R + 999*p["R_l"]))
    return c

vars1 = np.linspace(1,20000,50) * N_A ** -1 * 1e9
results = []
for var in vars1:
	sigma = np.sqrt(np.log(var**2/p["R_l"]**2 + 1))
	mean = np.log(p["R_l"]) - 1/2 * sigma**2
	temp = []
	for i in range(100):
		temp.append(np.mean(func(abs(np.random.lognormal(mean,sigma,1000)))))
	results.append(np.mean(temp)*1e9)
ax1 = sns.lineplot(vars1/p["R_l"], results, label="lognormal", color="green")
ax1.set(xlabel = "sigma/E", ylabel = "concentration in pM")
resultsLow = []

#low cell density
def funcLow(R):
    c = p["q_l"]/(p["k_on"]*R + 4*np.pi*p["D"]*p["rho"])
    return c
for var in vars1:
	sigma = np.sqrt(np.log(var**2/p["R_l"]**2 + 1))
	mean = np.log(p["R_l"]) - 1/2 * sigma**2
	temp = []
	for i in range(100):
		temp.append(np.mean(funcLow(abs(np.random.lognormal(mean,sigma,1000)))))
	resultsLow.append(np.mean(temp)*1e9)
# ax2 = sns.lineplot(vars1/p["R_l"], resultsLow, label="lognormLow")
# ax2.set(xlabel = "sigma/E", ylabel = "concentration in pM")
# print(results)
# plt.legend()
# plt.show()


import pandas as pd

group_variables = ["IL-2_sigma", None]
# subplot_style
D_to_iter_over = ['$10.0$']  # ['$1.0$', '$10.0$', '$100.0$']

common_y_axis = False
save_plots = False

# sns.set(rc={'figure.figsize':(5,5)})
# sns.set_style("ticks")
# sns.set_context("talk", font_scale=1, rc={"lines.linewidth": 5, 'figure.figsize':(6.4, 4.8)})
xscale = "linear"
yscale = "log"
if "IL-2_fraction" in group_variables:
    data_start = 5
    data_end = -2
else:
    data_start = 1
    data_end = None

plot_variables = ["std", "concentration"]  # ["surf_c_std"]#, "std", "gradient", "concentration", "surf_c"]
# plot_variables = ["surf_c", "surf_c_std"]  # ["surf_c_std"]#, "std", "gradient", "concentration", "surf_c"]
# plot_variables = ["gradient"]

get_dataframes = []  # [[]]*3
saving_dataframe = pd.DataFrame()
for j in range(14):
    get_dataframes.append([])
    print("reading in run", j)
    if group_variables[0] == "IL-2_fraction":
        # path = "/extra/brunner/10x10x10/q_fraction_exact/run" + str(j) + "/"
        path = "/extra/brunner/para_handling/static/q_fraction_large/run" + str(j) + "/"
        ext_cache = "/extra/brunner/10x10x10/q_fraction/run_ext_cache/"
    elif group_variables[0] == "IL-2_sigma":
        # path = "/extra/brunner/10x10x10/R_lognorm/run" + str(j) + "/"
        path = "/extra/brunner/para_handling/static/R_lognorm/run" + str(j) + "/"
        ext_cache = "/extra/brunner/10x10x10/R_lognorm/run_ext_cache/"
    else:
        print("Unknown grouping variable")
        exit()

    T = range(1)
    dt = 1

    global_df = pd.read_hdf(path + 'global_df.h5', mode="r")
    cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")
    # cell_stats_df = pd.read_hdf(path + 'cell_stats_dataframe_' + str(j) + '.h5', mode="r")


    get_dataframes[j] = [global_df, cell_df]

concat_df = pd.concat((get_dataframes[x][0] for x in range(len(get_dataframes))))

if "IL-2_sigma" in group_variables:
    plotting_df = concat_df.groupby(["IL-2_sigma"], as_index=False).mean()
    plotting_df["IL-2_D"] = "$10.0$"
else:
    plotting_df = concat_df.groupby(group_variables, as_index=False).mean()
plotting_df.reset_index(inplace=True)
# print(plotting_df["Concentration"].to_numpy())
plotting_df["sigma/E"] = plotting_df["IL-2_sigma"]/20000
plotting_df["Concentration"] = plotting_df["Concentration"]#*1e3
ax2 = sns.lineplot("sigma/E", "Concentration", data=plotting_df, label="modelled surface concentration", color="darkorange")
ax2.set(xlabel = "Heterogeneity", ylabel = "concentration in nM", xticks=[0,0.5,1])
# fig.savefig("/home/brunner/Documents/Current work/05062020/" + "interesting_behaviour" + ".png", bbox_inches='tight')
# plt.show()

fig = plt.figure()
ax1 = sns.distplot(np.abs(np.random.normal(20000,10000,10000)), label="normal")
var = 10000
sigma = np.sqrt(np.log(var ** 2 / 20000 ** 2 + 1))
mean = np.log(20000) - 1 / 2 * sigma ** 2
ax2 = sns.distplot(np.random.lognormal(mean,sigma,10000), label = "lognormal")
ax2.set(xlabel = "Receptors", ylabel = "probability")
plt.legend()
plt.show()
