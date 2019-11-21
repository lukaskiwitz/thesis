#!/home/kiwitz/anaconda3/envs/fenicsproject/bin/python

import multiprocessing as mp
import subprocess as sp
import os
from glob import glob
from parameter_scan_test import path_prefix
import importlib.util

def f(script: str) -> None:
    # print("mpirun -n 4 {c}/{s} run".format(s=script, c=os.getcwd()))
    # sp.run("mpirun -n 4 {c}/{s} run".format(s=script,c=os.getcwd()),shell=True)
    spec = importlib.util.spec_from_file_location(script,"{c}/{s}".format(s=script,c=os.getcwd()))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    # module.name
    print("{c}/post_process.py {p}{s}/".format(s=module.name, p=path_prefix, c=os.getcwd()))
    sp.run("{c}/post_process.py {p}{s}/".format(s=module.name,p=path_prefix,c=os.getcwd()),shell=True)


process_list = []

for file in glob("parameter_scans/*.py"):
    # process_list.append(mp.Process(target=f,args=(file,)))
    f(file)
#process_list = process_list[0:2]
# for i in process_list:
#     i.start()