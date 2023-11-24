import os
from multiprocessing import Pool
import json
from typing import List

import numpy as np

from parameters import path
from build_scans import get_scan_name,extract_scan_axis_from_dictionary, ScanAxis, chunks


def run(parameter_dict):
    scan_value_name,scan_scale,scan_kwargs = extract_scan_axis_from_dictionary(parameter_dict)

    join_path = "/".join(path.split("/")[0:-2]) + "/"
    os.makedirs(join_path,exist_ok=True)

    parameter_file_path = join_path + get_scan_name({scan_value_name:None}, **scan_kwargs) + ".json"
    with open(parameter_file_path, "w+") as fp:
        json.dump({k.value:v for k,v in parameter_dict.items()}, fp)

    os.system("python3 run.py "+parameter_file_path)


secreting_cells = np.array([2.5, 5, 10, 15, 20]) / 100
conditions = []
for sec in secreting_cells:
    conditions.append(
        {ScanAxis.gamma: [np.log10(1.1), np.log10(100), 15], ScanAxis.ftreg: 0, ScanAxis.fth: 0.999, ScanAxis.fsec: sec})


a = ["_".join([str(k)+"_"+str(v) if not isinstance(v,List) else str(k) for k,v in i.items()]) for i in conditions]
assert len(a)== len(np.unique(a))

n = 3
for sublst in list(chunks(conditions, n)):

    if n == 1:
        run(sublst[0])
    else:
        with Pool(len(sublst), maxtasksperchild=1) as pool:
            pool.map(run, sublst)