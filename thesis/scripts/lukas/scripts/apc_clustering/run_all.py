import json
import os
from multiprocessing import Pool
from typing import List

import numpy as np

from build_scans import get_scan_name, extract_scan_axis_from_dictionary, ScanAxis, chunks
from parameters import path


def run(parameter_dict):
    scan_value_name, scan_scale, scan_kwargs = extract_scan_axis_from_dictionary(parameter_dict)

    join_path = "/".join(path.split("/")[0:-2]) + "/"
    os.makedirs(join_path, exist_ok=True)

    parameter_file_path = join_path + get_scan_name({scan_value_name: None}, **scan_kwargs) + ".json"
    with open(parameter_file_path, "w+") as fp:
        json.dump({k.value: v for k, v in parameter_dict.items()}, fp)

    os.system("python run.py " + parameter_file_path)


NUMBER_OF_SAMPLES = 8
SAMPLES_IN_PARALLEL = 1

conditions = [
    {ScanAxis.treg_cs_strength: [0, 1, NUMBER_OF_SAMPLES], ScanAxis.ftreg: 0.2, ScanAxis.sigma: 1, ScanAxis.fsec: 0.1},
    {ScanAxis.treg_cs_strength: [0, 1, NUMBER_OF_SAMPLES], ScanAxis.ftreg: 0.025, ScanAxis.sigma: 1,
     ScanAxis.fsec: 0.1},
    {ScanAxis.treg_cs_strength: [0, 1, NUMBER_OF_SAMPLES], ScanAxis.ftreg: 0.05, ScanAxis.sigma: 1, ScanAxis.fsec: 0.1},
    {ScanAxis.treg_cs_strength: [0, 1, NUMBER_OF_SAMPLES], ScanAxis.ftreg: 0, ScanAxis.sigma: 1, ScanAxis.fsec: 0.1},
    {ScanAxis.sec_cs_strength: [0, 1, NUMBER_OF_SAMPLES], ScanAxis.ftreg: 0, ScanAxis.sigma: 1, ScanAxis.fsec: 0.01},
    {ScanAxis.sec_cs_strength: [0, 1, NUMBER_OF_SAMPLES], ScanAxis.ftreg: 0, ScanAxis.sigma: 1, ScanAxis.fsec: 0.05},
    {ScanAxis.sec_cs_strength: [0, 1, NUMBER_OF_SAMPLES], ScanAxis.ftreg: 0, ScanAxis.sigma: 1, ScanAxis.fsec: 0.1},
    {ScanAxis.sec_cs_strength: [0, 1, NUMBER_OF_SAMPLES], ScanAxis.ftreg: 0, ScanAxis.sigma: 1, ScanAxis.fsec: 0.2},
]

a = ["_".join([str(k) + "_" + str(v) if not isinstance(v, List) else str(k) for k, v in i.items()]) for i in conditions]
assert len(a) == len(np.unique(a))

for sublst in list(chunks(conditions, SAMPLES_IN_PARALLEL)):

    if SAMPLES_IN_PARALLEL == 1:
        run(sublst[0])
    else:
        with Pool(len(sublst), maxtasksperchild=1) as pool:
            pool.map(run, sublst)
