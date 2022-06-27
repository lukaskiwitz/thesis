import os
from copy import deepcopy

from build_scans import compute_samples_points, build_scan_function, ScanAxis, extract_scan_axis_from_dictionary
from parameters import cytokines, cell_types_dict, geometry, numeric
from parameters import ext_cache
from parameters import path as old_path
from thesis.main.MyScenario import MyScenario
from thesis.scenarios.box_grid import setup

scenario: MyScenario = setup(deepcopy(cytokines), deepcopy(cell_types_dict), [], deepcopy(geometry),
                             deepcopy(numeric), deepcopy(old_path))

scan_kwargs = {ScanAxis.treg_cs_strength: 0.9, ScanAxis.ftreg: [0.1, 0.1, 1], ScanAxis.sigma: 1, ScanAxis.fsec: 0.1}

scan_value_name, scan_scale, scan_kwargs = extract_scan_axis_from_dictionary(scan_kwargs)
scan_samples, scan_name = build_scan_function(scenario, scan_value_name, scan_scale, **scan_kwargs)

old_path = os.path.join("/".join(old_path.split("/")[:-2]), "meshing_run")
compute_samples_points(scan_samples[0:1], scenario, old_path, [0, 1], ext_cache=ext_cache, number_of_replicats=1)
