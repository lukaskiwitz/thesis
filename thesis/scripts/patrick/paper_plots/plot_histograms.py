import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import random
from my_debug import message
from scipy import spatial
from scipy.optimize import curve_fit

#load the cell_df
# from parameters_q_fraction import path_kinetic
# path = path_kinetic
# path = "/extra/brunner/thesis/kinetic/q_fraction/"
# path = "/extra/brunner/thesis/kinetic/q_fraction_gamma_scan/"
# path = "/extra/brunner/thesis/kinetic/q_fraction_wave_multi/q_fraction_wave_multi_0/"
# path = "/extra/brunner/thesis/kinetic/q_fraction_large_gamma_scan_new_paras/"
# path = "/extra/brunner/thesis/kinetic/q_fraction_test/"
# path = "/extra/brunner/thesis/kinetic/q_fraction_small_gamma_scan_g_0.1/"
# path = "/extra/brunner/thesis/kinetic/q_fraction_small_g_0.1/"


verbose = False
save_avg_hist_df = True

no_of_runs = 1
dt = 1

cell_cell_distance = 20

# base_path = "/extra/brunner/thesis/kinetic/q_fraction_medium_pos_g_scan_multi/"
base_path = "/extra/brunner/thesis/kinetic/kin_large_new_hill_4/"

runs_hist_df = pd.read_hdf(base_path + str(no_of_runs) + '_runs_0.25_f_hist_df.h5', mode="r")
avg_hist_df = pd.read_hdf(base_path + str(no_of_runs) + '_runs_0.25_f_avg_hist_df.h5', mode="r")
mean_sec_dist = float(np.loadtxt(base_path + str(no_of_runs) + "_runs_0.25_f_avg_sec_dists.txt"))/cell_cell_distance

global_df = pd.read_hdf(base_path + "global_df.h5", mode="r")
cell_df = pd.read_hdf(base_path + "cell_df.h5", mode="r")
fractioned_cell_df = cell_df.loc[(cell_df["IL-2_Tsec_fraction"] == 0.25)]

message("loaded runs_hist_df from " + base_path)

timepoints = list(runs_hist_df.index[:-1])
timepoints = [timepoints[x] for x in [5]]
# timepoints = np.array([0])*dt
# gammas = [runs_hist_df.columns[x][1] for x in range(len(runs_hist_df.columns))]
gammas = [1/10]

#############################################################################################################################################

def AC_plotting(avg_hist_df, gamma, time, color):
    hist = avg_hist_df[gamma][time]
    auto_corr = []
    for r_dash in range(len(hist) - 1):
        R = 0
        for r in range(len(hist)):
            try:
                R += hist[r + r_dash] * hist[r]
            except IndexError:
                pass
        auto_corr.append(R)
    AC_niche = []
    def func(x, a, b, c):
        return a * b * np.exp(-x / b) + c

    popt, pcov = curve_fit(func, np.arange(len(auto_corr)), auto_corr, maxfev=10000)
    AC_niche.append(popt[1])

    # ax1 = sns.lineplot(np.arange(0,len(avg_hist_df[gamma][time])), avg_hist_df[gamma][time], label=variable + "=" + str(time*10) + "h", legend=myLegend, color=palette[i], marker="o")
    sns.lineplot(np.arange(len(auto_corr)), func(np.arange(len(auto_corr)), *popt), color=color)
    return AC_niche

def histogram_plotting(avg_hist_df, gamma, time, color):
    plt.plot(avg_hist_df[gamma][time], label=str(round(time/3600,1)) + "h", color=color, marker="o")


# Niche
# Fitting

# AC_niche = []
#
# def func(x, a, b, c):
#     return a * b * np.exp(-x/b) + c
# from scipy.optimize import curve_fit
# for i in range(len(auto_corr_data_array)):
#     popt, pcov = curve_fit(func, np.arange(len(auto_corr_data_array[i])), auto_corr_data_array[i], maxfev=10000)
#     AC_niche.append(popt[1])


# Thresholding

# niche = []
#
# for j, hist in enumerate(hist_data_array):
#     if j == 0:
#         threshold = 200
#     for i,c in enumerate(hist):
#         if c < threshold:
#             niche.append(i)
#             break
#         if i == len(hist)-1:
#             niche.append(i)

# ax13 = sns.lineplot(points, niche, label="hist_niche", legend=False)
# ax13 = sns.lineplot(points[:], AC_niche[:], label="AC_niche")
# ax13.set(xlabel="time", ylabel="cell-cell-distance", yscale="linear", xscale="linear", title="Niche size at threshold " + str(threshold))
# plt.show()

#############################################################################################################################################

# Plotting


fig = plt.figure(figsize=(12,12))
sns.set_style("ticks")
sns.set_context("talk", font_scale=2, rc={"lines.linewidth": 3, 'lines.markersize': 5})
plt.subplots_adjust(wspace=.4, hspace=0.6)
from math import ceil
if len(gammas) > 1:
    # a_x = ceil(len(slope_thresholds)/3)
    # a_y = ceil(len(slope_thresholds)/a_x)
    a_x = ceil(np.sqrt(len(gammas)))
    a_y = ceil(len(gammas)/a_x)
else:
    a_x = 1
    a_y = 1


xscale = "linear"
yscale = "linear"

if np.isnan(mean_sec_dist) == True:
    mean_sec_dist = 0

variable = "t"

for i, gamma in enumerate(gammas):
    if gamma == 0.1:
        gamma_str = "str.neg."
        myLegend = None
    elif gamma == 0.5:
        gamma_str = "neg."
        myLegend = "brief"
    elif gamma == 2:
        gamma_str = "pos."
        myLegend = None
    elif gamma == 10:
        gamma_str = "str.pos."
        myLegend = None
    else:
        myLegend = None
        gamma_str = ""

    fig.add_subplot(a_x, a_y, i+1)

    palette = sns.color_palette("bwr", len(timepoints)*2)
    palette = palette[len(timepoints):]
    palette = ["red"]

    for i,time in enumerate(timepoints):
        # AC_plotting(avg_hist_df, gamma, time, color=palette[i])
        histogram_plotting(avg_hist_df,gamma,time,palette[i])

    # ax2 = sns.lineplot(center[:len(hist_data_array[1])], hist_data_array[1], label=str(float(time_list[1])*10)+"h", legend=False)
    # ax3 = sns.lineplot(center[:len(hist_data_array[2])], hist_data_array[2], label=str(float(time_list[2])*10)+"h", legend=False)
    plt.axvline(mean_sec_dist, color="black", linestyle="dashed")
    if gamma < 1:
        plt.title("g = 1/" + str(int(1/gamma)))
    else:
        plt.title("g = " + str(gamma))
    # ax1.set(xlabel="cell-cell-distance", ylabel="Avg. c. in pM", yscale=yscale, xscale=xscale, xticks=[0,5,10], title="g = " + str(gamma)+" , "+gamma_str)
    plt.xticks([0,2,4,6,8,10,12])
    plt.xlabel("c-c-d")
    plt.ylabel("Concentration (pM)")
    plt.legend()
plt.show()

# from scipy.constants import N_A
# p = {
#     "kon": 100 * 1e9 / 60 ** 2,  #1e9 *  111.6 / 60 ** 2,  # 111.6 per hour # now 540 per hour per nM
#     "D": 10,  # mu² per s -> m
#     "R": 1e4 * N_A ** -1 * 1e9,
#     "R_Treg": 1e4 * N_A ** -1 * 1e9,
#     "kd": 0.1/(60**2), # 1/s
#     "q": 10 * N_A ** -1 * 1e9,
#     "gamma": 0.01,
#     "L": 20e-6, #m
#     "rho": 5e-6, #m
#     "N": 1000,
# }
#
#
# def low_density(r, q, rho, kon,R,D):
#     return q*rho/(r*(kon*R + 4*np.pi*D*rho))*1e12
#
# def high_density(r, q, rho, kon, R, D, L, N, Rresp):
#     c = q*rho/(kon*r) * (4*np.pi*D*r*(L+rho) + kon*(L-r+rho)*N*Rresp)/(kon*L*R*N*Rresp + 4*np.pi*D*rho*(L+rho)*(R+N*Rresp))
#     c = 4*np.pi*D*r*(L+rho)
#     return c*1e12
#
# def test(r, q, rho, kon, R, D, L, N, Rresp):
#     c = 4*np.pi*D*r*(L+rho) + kon*(L-r+rho)*N*Rresp
#     return c*1e12
#
# myRange = np.linspace(p["rho"],16*p["L"] + p["rho"], 16)
# sns.set()
# # plt.plot(low_density(myRange, p["q"], p["rho"], p["kon"], p["R"], p["D"]), label="low density")
# # plt.plot(high_density(myRange, p["q"], p["rho"], p["kon"], p["R"], p["D"], p["L"], p["N"], p["R_Treg"]), label="high_density")
# plt.plot(test(myRange, p["q"], p["rho"], p["kon"], p["R"], p["D"], p["L"], p["N"], p["R_Treg"]), label="high_density")
# # plt.legend()
# plt.show()
# fig.savefig("plots/" + "averaged_histograms" + ".svg", bbox_inches='tight')
# fig.savefig("plots/" + "averaged_histograms" + ".png", bbox_inches='tight')
exit()




























fig = plt.figure(figsize=(9,8))
sns.set_context("talk", font_scale=1.7, rc={"lines.linewidth": 5})

# plt.figure(figsize=(12,10))
# sns.set_context("talk", font_scale=1.1, rc={"lines.linewidth": 2.5})
# sns.set(rc={'figure.figsize':(8,8)})
sns.set_style("ticks")
sns.set_context("talk", font_scale=1.7, rc={"lines.linewidth": 5, 'lines.markersize': 10})


time = 20
palette = sns.color_palette("bwr", 4)

from scipy.optimize import curve_fit

def func(x, a,b,c):
    return a*np.exp(-x/b) + c


sns.lineplot(np.arange(0,len(avg_hist_df[0.1][time])), avg_hist_df[0.1][time], label="strong neg.", color=palette[3], marker="o")
popt = curve_fit(func, np.arange(0,len(avg_hist_df[0.1][time])), avg_hist_df[0.1][time], maxfev=10000)[0]
sns.lineplot(np.arange(0,len(avg_hist_df[0.1][time])), func(np.arange(0,len(avg_hist_df[0.1][time])), *popt), color=palette[3], linewidth=2)

sns.lineplot(np.arange(0,len(avg_hist_df[0.5][time])), avg_hist_df[0.5][time], label="neg.", color=palette[2], marker="o")
popt = curve_fit(func, np.arange(0,len(avg_hist_df[0.5][time])), avg_hist_df[0.5][time], maxfev=10000)[0]
sns.lineplot(np.arange(0,len(avg_hist_df[0.5][time])), func(np.arange(0,len(avg_hist_df[0.5][time])), *popt), color=palette[2], linewidth=2)


sns.set_context("talk", font_scale=1.7, rc={"lines.linewidth": 3, 'lines.markersize': 5})
sns.lineplot(np.arange(0,len(avg_hist_df[0.5][0])), avg_hist_df[0.5][0], label="no", color="black", marker="o")
sns.set_context("talk", font_scale=1.7, rc={"lines.linewidth": 5, 'lines.markersize': 10})


sns.lineplot(np.arange(0,len(avg_hist_df[2][time])), avg_hist_df[2][time], label="pos.", color=palette[1], marker="o")
popt = curve_fit(func, np.arange(0,len(avg_hist_df[2][time])), avg_hist_df[2][time], maxfev=10000)[0]
sns.lineplot(np.arange(0,len(avg_hist_df[2][time])), func(np.arange(0,len(avg_hist_df[2][time])), *popt), color=palette[1], linewidth=2)

sns.lineplot(np.arange(0,len(avg_hist_df[10][time])), avg_hist_df[10][time], label="strong pos.", color=palette[0], marker="o")
popt = curve_fit(func, np.arange(0,len(avg_hist_df[10][time])), avg_hist_df[10][time], maxfev=10000)[0]
sns.lineplot(np.arange(0,len(avg_hist_df[10][time])), func(np.arange(0,len(avg_hist_df[10][time])), *popt), color=palette[0], linewidth=2)



plt.xlabel("cell-cell-distance")
plt.ylabel("Avg. c. in pM")
plt.xticks([0,2,4,6,8,10])

# ax2 = sns.lineplot(center[:len(hist_data_array[1])], hist_data_array[1], label=str(float(time_list[1])*10)+"h", legend=False)
# ax3 = sns.lineplot(center[:len(hist_data_array[2])], hist_data_array[2], label=str(float(time_list[2])*10)+"h", legend=False)
import matplotlib.lines as mlines

zero_line = mlines.Line2D([], [], color=palette[0], label='lightcoral line')
one_line = mlines.Line2D([], [], color=palette[1], label='red line')
black_line = mlines.Line2D([], [], color='black', label='Black line')
two_line = mlines.Line2D([], [], color=palette[2], label='lightsteelblue line')
three_line = mlines.Line2D([], [], color=palette[3], label='blue line')

grey_line = mlines.Line2D([], [], color="grey", label='grey line')
grey_dashed_line = mlines.Line2D([], [], color="grey", label='grey dashed line', linestyle="dashed")
white_line = mlines.Line2D([], [], color="white", label='white line')


new_handles = [three_line, black_line, zero_line]
new_labels = ["neg. feedback", "no feedback", "pos. feedback"]
# new_handles = [three_line, two_line, black_line, one_line, zero_line]
# new_labels = ["strong negative feedback", "negative feedback", "no feedback", "positive feedback", "strong positive feedback"]
plt.legend(handles=new_handles, labels=new_labels,loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True)
# fig.savefig("plots/" + "averaged_histograms_t_200" + ".svg", bbox_inches='tight')
plt.show()
# plt.plot(niche, label="c_niche")
# plt.plot(AC_niche, label="AC_niche")
# plt.legend()
# plt.show()

# # fig.add_subplot(a_x, a_y, 2)
# # for i in range(len(c_data_array)):
# #     ax4 = sns.lineplot(center[:len(c_data_array[i])], c_data_array[i], label=variable + "=" + str(float(points[i])), legend=False)
# # # ax5 = sns.lineplot(center[:len(c_data_array[1])], c_data_array[1], label=str(float(time_list[1])*10)+"h", legend=False)
# # # ax6 = sns.lineplot(center[:len(c_data_array[2])], c_data_array[2], label=str(float(time_list[2])*10)+"h", legend=False)
# # plt.axvline(mean_sec_dist, label="sec. avg. dist.", color="black")
# # ax4.set(xlabel="cell-cell-distance", yscale=yscale, xscale=xscale, title="c_0*c_r", xticks=[0,5,10])
#
# # fig.add_subplot(a_x,a_y,3)
# # for i in range(len(conv_data_array)):
# #     ax7 = sns.lineplot(center[:len(conv_data_array[i])], conv_data_array[i], label=variable + "=" + str(float(points[i])), legend=False)
# # # ax8 = sns.lineplot(center[:len(conv_data_array[1])], conv_data_array[1], label=str(float(time_list[1])*10)+"h", legend=False)
# # # ax9 = sns.lineplot(center[:len(conv_data_array[2])], conv_data_array[2], label=str(float(time_list[2])*10)+"h", legend=False)
# # plt.axvline(mean_sec_dist, label="sec. avg. dist.", color="black")
# # ax7.set(xlabel="cell-cell-distance", ylabel="Convolution", yscale=yscale, xscale=xscale, title="conv(c_r * c_r)", xticks=[0,5,10])
# #
# fig.add_subplot(a_x,a_y,2)
# for i in range(len(auto_corr_data_array)):
#     ax10 = sns.lineplot(center[:len(auto_corr_data_array[i])], auto_corr_data_array[i], label=variable + "=" + str(float(points[i])*10) + "h", legend=False)
# # ax11 = sns.lineplot(center[:len(auto_corr_data_array[1])], auto_corr_data_array[1], label=str(float(time_list[1])*10)+"h", legend=False)
# # ax12 = sns.lineplot(center[:len(auto_corr_data_array[-1])], auto_corr_data_array[-1], label=str(float(time_list[-1])*10)+"h", legend=False)
#
# ax10.set(xlabel="cell-cell-distance", ylabel="G(r)", yscale=yscale, xscale=xscale, title="AC with r")
# # plt.axvline(mean_sec_dist, label="sec. avg. dist.", color="black")
# #
# # fig.add_subplot(a_x, a_y, 3)
# # for i in range(len(sec_cell_corr_array)):
# #     ax14 = sns.lineplot(center[:len(sec_cell_corr_array[i])], sec_cell_corr_array[i], label=variable + "=" + str(float(points[i])), legend=False)
# # # ax14.set(xlabel="cell-cell-distance", ylabel="G(r)", yscale=yscale, xscale=xscale, title="iterating AC", xticks=[0,5,10])
# #
# # plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
# #
# #
# fig.add_subplot(a_x,a_y,3)
# ax13 = sns.lineplot(np.array(points)*10, niche, label="hist_niche", legend=False)
# ax13 = sns.lineplot(np.array(points[:])*10, AC_niche[:], label="AC_niche")
# ax13.set(xlabel="time in h", ylabel="cell-cell-distance", yscale=yscale, xscale=xscale, title="Niche size at threshold " + str(threshold))
#
#
# # fig.savefig("/home/brunner/Documents/Current work/28082020/" + "niche_" + "f=" + str(cell_df["IL-2_fraction"].unique()[0]) + "_t=" + str(threshold) + ".pdf", bbox_inches='tight')
# # fig.savefig("/home/brunner/Documents/Current work/28082020/" + "niche_" + "f=" + str(cell_df["IL-2_fraction"].unique()[0]) + "_t=" + str(threshold) + ".png", bbox_inches='tight')
#
# plt.show()