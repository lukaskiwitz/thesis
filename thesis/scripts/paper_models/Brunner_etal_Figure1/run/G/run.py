import json
import os
import sys
from copy import deepcopy

import numpy as np
import sympy

from build_scans import compute_samples_points, post_process, get_new_samples_points, build_scan_function,\
    get_existing_path_list,extract_scan_axis_from_dictionary,ScanAxis
from multiprocessing import Process
from parameters import cytokines, cell_types_dict, geometry, numeric
from parameters import ext_cache
from parameters import path as old_path
from thesis.main.MyScenario import MyScenario
from thesis.scenarios.box_grid import setup

ext_cache = "../" + ext_cache

T = 100
dt = 3

start = 0
stop = 2
grad_max = 1e5

dt = 3600  # 1h
ramp_length = 8
ramp_max_T = dt * 4

T_length = 28
max_T = dt * 48

T_ramp = np.linspace(0,ramp_max_T, ramp_length)
T = np.linspace(np.max(T_ramp), max_T, T_length)
# T = np.linspace(0, max_T, T_length)

# time_range = np.concatenate([T_ramp[:-1], T])
time_range = [0,1]

scenario: MyScenario = setup(deepcopy(cytokines), deepcopy(cell_types_dict), [], deepcopy(geometry),
                             deepcopy(numeric), deepcopy(old_path))

# scan_kwargs = {ScanAxis.treg_cs_strength: 1, ScanAxis.ftreg: 0, ScanAxis.sigma: 1,ScanAxis.sec_q:[1,30,6]}

if __name__ == "__main__":
    for a in sys.argv:
        print(a)

    if len(sys.argv)  == 2:

        with open(str(sys.argv[1]),"r") as fp:
            scan_kwargs = json.load(fp)


scan_value_name,scan_scale,scan_kwargs = extract_scan_axis_from_dictionary(scan_kwargs)

scan_samples, scan_name = build_scan_function(scenario, scan_value_name,scan_scale, **scan_kwargs)
old_path = old_path[:-1] + "_" + scan_name
path_list,start = get_existing_path_list(old_path, start=start)

for i in range(start,stop+1):

    if len(path_list) > 0:
        print("calculating sample based on step {i}".format(i = i-1))
        scan_scale = get_new_samples_points([item for sublist in path_list for item in sublist],i-1, grad_max=grad_max)
        scan_samples, scan_name = build_scan_function(scenario, scan_value_name, scan_scale, **scan_kwargs)

    if i == stop: break
    if len(scan_samples) == 0: break
    if len(scan_samples) > 64: break

    path_list.append([os.path.join(old_path, "refine_{i}_{o}".format(i=i, o=o)) for o, s in enumerate(scan_samples)])

    processes = []
    for o, sample in enumerate(scan_samples):
        p = Process(target=compute_samples_points,
                    args=([sample], scenario, path_list[-1][o], time_range,),
                    kwargs={"ext_cache": ext_cache}
                    )
        p.start()
        processes.append(p)
    for p in processes:
        p.join()
        p.close()

    # processes = []
    # n = int(8 / len(path_list[-1]))
    # n = n if n > 0 else 1
    #
    # for o, path in enumerate(path_list[-1]):
    #     p = Process(target=post_process,
    #                 args=(path, [0,len(time_range) - 2]),
    #                 kwargs={"n_procs": n}
    #                 )
    #     p.start()
    #     processes.append(p)
    # for p in processes:
    #     p.join()
    #     p.close()