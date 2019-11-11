#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 11:19:45 2019

@author: kiwitz

"""
import os
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from typing import cast, List, Tuple

# import matplotlib as mpl
# mpl.use('Agg')

sns.set_context("paper", font_scale=1.5, rc={
    "lines.linewidth": 1,
    "lines.markersize": 20,
    'xtick.labelsize': 'small',
    'ytick.labelsize': 'small',
    "xtick.major.width": 0.8,
    "xtick.minor.width": 0.8,
    "xtick.major.size": 3.5,
    "xtick.minor.size": 3.5})  #


def group_count(data_frame: pd.DataFrame, label: str) -> pd.DataFrame:
    data_frame_g = data_frame.groupby(["x", "time"], as_index=False)
    series: pd.DataFrame = pd.DataFrame([data_frame_g.size().to_numpy()], ["n"]).T
    data_frame = data_frame_g.mean().join(series)
    data_frame = data_frame.assign(l=label)
    return data_frame


def make_simple_fractions(key: str, t_hold: float, label: str, data_in: pd.DataFrame) -> pd.DataFrame:
    data_frame = data_in.loc[(data_in[key] > t_hold)]
    return group_count(data_frame, label)


def make_conditional_fraction(keys: List[str], t_holds: List[float], label: str, bool_tuple: Tuple[bool, bool], data_in: pd.DataFrame):
    if not bool_tuple[0] and not bool_tuple[1]:
        # noinspection PyUnusedLocal
        data_frame: pd.DataFrame = data_in.loc[
            (data_in[keys[0]] < t_holds[0]) & (data_in[keys[1]] < t_holds[1])
            ]
    if bool_tuple[0] and bool_tuple[1]:
        # noinspection PyUnusedLocal
        data_frame: pd.DataFrame = data_in.loc[
            (data_in[keys[0]] > t_holds[0]) & (data_in[keys[1]] > t_holds[1])
            ]
    if not bool_tuple[0] and bool_tuple[1]:
        data_frame: pd.DataFrame = data_in.loc[
            (data_in[keys[0]] < t_holds[0]) & (data_in[keys[1]] > t_holds[1])
            ]
    else:
        data_frame: pd.DataFrame = data_in.loc[
            (data_in[keys[0]] > t_holds[0]) & (data_in[keys[1]] < t_holds[1])
            ]
    return group_count(data_frame, label)


def make_average(data_in: pd.DataFrame) -> pd.DataFrame:
    return data_in.groupby(["x", "time"], as_index=False).mean()


def plot_averages(img_path: str, x_ticks: List[float], x_scale: List[float], data_frame: pd.DataFrame) -> None:
    average = make_average(data_frame)

    color = "tab:red"
    fig, ax1 = plt.subplots()
    del fig
    sns.lineplot(x="x", y="il2", ci="sd",
                 data=average,
                 ax=ax1,
                 color=color,
                 markers=["."],
                 legend=False,
                 style=True,
                 err_style="bars")

    ax1.set_ylim(0)
    ax1.tick_params(axis="y", labelcolor=color)
    ax1.set_ylabel(r'IL-2 on cell surface (nM)', color=color)
    ax2 = ax1.twinx()
    color = "tab:blue"
    sns.lineplot(x="x", y="il6", ci="sd",
                 data=average,
                 ax=ax2,
                 color=color,
                 markers=["."],
                 legend=False,
                 style=True,
                 err_style="bars")

    ax2.tick_params(axis="y", labelcolor=color)
    ax2.set_ylabel(r'IL-21 on cell surface (nM)', color=color)
    ax2.set_ylim(0)
    ax1.set_xlabel(r'distance from left boundary $[\mu m]$')
    ax1.set_xlabel(r'cell distance from boundary')
    plt.xticks(x_ticks, x_scale)
    plt.tight_layout()
    plt.savefig(img_path + "/" + "averageConcentrations.pdf", dpi=1200)
    plt.close()


def plot_simple_fractions(img_path: str, x_ticks: List[float], x_scale: List[float], data_frame: pd.DataFrame, t_holds: List[float]):
    il2p = make_simple_fractions("il2", t_holds[0], "IL-2", data_frame)
    il6p = make_simple_fractions("il6", t_holds[1], "IL-6", data_frame)
    simple_fractions = il2p.append(il6p)

    #    color = {"IL-2": "tab:red", "IL-6": "tab:blue"}
    fig, ax1 = plt.subplots()
    del fig
    sns.lineplot(x="x", y="n", ci="sd", data=simple_fractions,
                 markers={"IL-2": ".", "IL-6": "."},
                 style="l",
                 hue="l",
                 dashes={"IL-2": [1, 0], "IL-6": [1, 0]})

    ax1.set_ylim(0)
    ax1.tick_params(axis="y")
    ax1.set_ylabel(r'Fraction of positive cells (%)')
    ax1.set_xlabel(r'cell distance from boundary')
    plt.xticks(x_ticks, x_scale)
    plt.tight_layout()
    plt.savefig(img_path + "/" + "simple_fractions.pdf", dpi=1200)
    plt.close()


def plot_conditional_fractions(img_pat_holds, x_ticks, x_scale, data_frame, t_holds) -> None:
    il2pil6n = make_conditional_fraction(
        ["il2", "il6"], t_holds, "IL2$^+$ IL21$^-$", (True, False), data_frame)
    il2nil6p = make_conditional_fraction(
        ["il2", "il6"], t_holds, "IL2$^-$ IL21$^+$", (False, True), data_frame)
    il2pil6p = make_conditional_fraction(
        ["il2", "il6"], t_holds, "IL2$^+$ IL21$^+$", (True, True), data_frame)
    il2nil6n = make_conditional_fraction(
        ["il2", "il6"], t_holds, "IL2$^-$ IL21$^-$", (False, False), data_frame)

    conditional_fractions = il2pil6n
    for i in [il2nil6p, il2pil6p, il2nil6n]:
        conditional_fractions = conditional_fractions.append(i)

    sns.lineplot(x="x", y="n", ci="sd", data=conditional_fractions,
                 hue="l", err_style="bars", marker=".")
    plt.xlabel(r'distance from left boundary $[\mu m]$')
    plt.xlabel(r'cell distance from boundary')
    plt.ylabel(r'fraction of positive cells (%)')
    plt.ylim(-5, 105)
    plt.yticks(np.arange(0, 120, 20), np.arange(0, 120, 20))
    plt.xticks(x_ticks, x_scale)
    plt.tight_layout()
    plt.savefig(img_pat_holds + "/" + "conditional_fractions.pdf", dpi=600)
    plt.close()


def single_scan_plots(data: pd.DataFrame, img_path: str) -> None:
    scan_groups = data.groupby(by=["scanIndex"])
    t_holds = [0.1, 2]
    keys = list(scan_groups.groups.keys())
    data_frame = scan_groups.get_group(keys[0])

    x_ticks = list(np.unique(data_frame["x"].to_numpy()))
    x_scale = [0, "", 2, "", 5, "", 6, "", 8, "", 10]

    for key in keys:
        scan_index = data_frame.get("scanIndex").values[0]
        img_path = img_path + str(scan_index)
        os.makedirs(img_path, exist_ok=True)
        data_frame = scan_groups.get_group(key)

        # plotting
        plt.figure(1)
        plot_averages(img_path, x_ticks, x_scale, data_frame)

        plt.figure(2)
        plot_simple_fractions(img_path, x_ticks, x_scale, data_frame, t_holds)

        plt.figure(3)
        plot_conditional_fractions(img_path, x_ticks, x_scale, data_frame, t_holds)


def parameter_plot(data: pd.DataFrame, img_path: str) -> None:
    # x = np.unique(data.filter("D").to_numpy())
    # data = data.div(x[3],axis="D")
    plt.figure(1)
    sns.lineplot(x="scanIndex", y="concentration", ci="sd", data=data,
                 hue="l", err_style="bars", marker=".")
    # plt.xlabel(r'$D$')
    # plt.ylabel(r'concentration (nM)')
    plt.ylabel(r'concentrations')
    # plt.legend(["IL-2", "IL-6"], title="field")
    plt.tight_layout()
    # plt.savefig(img_path+"/"+"global_concentration.pdf", dpi=600)
    plt.ylim([0,20])
    plt.show()
    # plt.close()


PATH = "/extra/kiwitz/results_parameter_scan_Diffusion/"
# PATH = "/home/chrx/results_parameter_scan_Diffusion/"

IMG = "/home/kiwitz/postProcessResult_img/"

# tree = ET.parse(pat_holds+"postProcess.xml")
# res = _prepData(tree)

PATH_LIST = [
    {"path":"/extra/kiwitz/results_parameter_scan_Diffusion/","key":"D"},
    {"path":"/extra/kiwitz/results_parameter_scan_fraction/","key":"fraction"},
    {"path":"/extra/kiwitz/results_parameter_scan_kd/","key":"decay"},
    {"path":"/extra/kiwitz/results_parameter_scan_kON/","key":"k_{on}"},
    {"path":"/extra/kiwitz/results_parameter_scan_q_s/","key":"q Secretors"},
    {"path":"/extra/kiwitz/results_parameter_scan_R_il2_f/","key":"R il2 f"},
    {"path":"/extra/kiwitz/results_parameter_scan_R_il2_s/","key":"R il2 s"}
             ]

res_global: pd.DataFrame = pd.DataFrame()
for e in PATH_LIST:
    path = e["path"]
    frame: pd.DataFrame = cast(pd.DataFrame, pd.read_hdf(path + 'global_dataframe.h5', "data"))
    frame = frame.assign(l = e["key"])
    frame = frame.loc[frame["fieldName"] == "il2"]
    res_global= res_global.append(frame)

o = parameter_plot(res_global, IMG)

# RES: pd.DataFrame = cast(pd.DataFrame, pd.read_hdf(PATH + 'dataframe.h5', "data"))
# single_scan_plots(RES, IMG)

