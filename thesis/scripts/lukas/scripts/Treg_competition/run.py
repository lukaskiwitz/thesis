import os
import sys
from copy import deepcopy

import numpy as np

from build_scans import compute_samples_points, post_process, get_new_samples_points
from parameters import cytokines, cell_types_dict, geometry, numeric
from parameters import ext_cache
from parameters import path as old_path
from thesis.main.MyScenario import MyScenario
from thesis.scenarios.box_grid import setup
import os
import sys
from copy import deepcopy

import numpy as np

from build_scans import compute_samples_points, post_process, get_new_samples_points
from parameters import cytokines, cell_types_dict, geometry, numeric
from parameters import ext_cache
from parameters import path as old_path
from thesis.main.MyScenario import MyScenario
from thesis.scenarios.box_grid import setup

ext_cache = "../" + ext_cache

T = 200
dt = 3
time_range = np.linspace(0, T, int(np.ceil(T / dt))) * 3600

scenario: MyScenario = setup(deepcopy(cytokines), deepcopy(cell_types_dict), [], deepcopy(geometry),
                             deepcopy(numeric), deepcopy(old_path))


class ParameterError(Exception): pass


if __name__ == "__main__":

    cs = float(sys.argv[1])
    if not (cs >= 0 and cs <= 1):
        raise ParameterError

    ec50 = float(sys.argv[2])
    if ec50 < 0:
        raise ParameterError

    k = float(sys.argv[3])
    if k < 0:
        raise ParameterError

    treg_half_time = float(sys.argv[4])
    if treg_half_time <= 0:
        raise ParameterError

    scan_kwargs = {"cs_strength": cs, "pSTAT_ec50": ec50, "treg_k": k, "treg_half_time": treg_half_time}

else:
    scan_kwargs = {"cs_strength": 1, "pSTAT_ec50": 1, "treg_k": 860}

from build_scans import scan_2
from multiprocessing import Process

scan_samples, scan_name = scan_2(scenario, np.linspace(0.01, 0.5, 16), **scan_kwargs)
old_path = old_path[:-1] + "_" + scan_name

path_list = []
for i in [0, 1]:

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

    processes = []
    n = int(8 / len(path_list[-1]))
    n = n if n > 0 else 1

    for o, path in enumerate(path_list[-1]):
        p = Process(target=post_process,
                    args=(path, [len(time_range) - 2]),
                    kwargs={"n_procs": n}
                    )
        p.start()
        processes.append(p)

    for p in processes:
        p.join()
        p.close()

    scan_scale = get_new_samples_points([item for sublist in path_list for item in sublist], grad_max=0.01)
    scan_samples, scan_name = scan_2(scenario, scan_scale, **scan_kwargs)
