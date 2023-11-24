import glob
import os
from multiprocessing import Pool
from build_scans import get_existing_path_list, chunks
import tqdm
import gc
import re

from parameters import path, scan_name
print(path)

def to_rgb(h):
    return [int(h[2*i:2*i+2],16)/255 for i in range(3)]

def run(path):

    if path is not None:
        if TIME_INDEX is not None:
            os.system("python post_process.py {p} {n} {ti}".format(p=path + "/", n = processes_per_sample, ti = int(TIME_INDEX)))
        else:
            os.system("python post_process.py {p} {n}".format(p=path + "/", n=processes_per_sample))


all_refine_paths = []
path_list = glob.glob(os.path.join("/".join(path.split("/")[:-2]),scan_name+"*"))
# path_list = [
#     [re.search(r".*",p) for p in path_list]
# ]

for p in path_list:
    all_refine_paths = all_refine_paths + get_existing_path_list(p, start = None)[0]
all_refine_paths = [item for sublist in all_refine_paths for item in sublist]

new_refine_paths = []
for entry in all_refine_paths:
    if not ("cell_df.h5" in os.listdir(entry)) or ("cell_df.pkl" in os.listdir(entry)):
        new_refine_paths.append(entry)

print(len(all_refine_paths))
print(len(new_refine_paths))
if len(new_refine_paths) == 0:
    new_refine_paths = all_refine_paths
    print("length was zero, set to full redo")

samples_in_parallel = 3
processes_per_sample = 15
TIME_INDEX = 57
if TIME_INDEX is not None:
	print(f"time_index set to {TIME_INDEX}, only steadystate is calculated!")

# chks = list(chunks(all_refine_paths, samples_in_parallel))
chks = list(chunks(new_refine_paths, samples_in_parallel))
for sublst in tqdm.tqdm(chks):

    if samples_in_parallel == 1:
        run(sublst[0])
    else:
        with Pool(samples_in_parallel, maxtasksperchild=1) as pool:
            pool.map(run,sublst)
    gc.collect()