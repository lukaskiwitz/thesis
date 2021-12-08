import glob
import os
import re
import sys
from multiprocessing import Pool

import tqdm

from parameters import path, scan_name

# logging.basicConfig(
#     filename=os.path.join(".", "post_process.log"),
#     level=logging.INFO,
#     filemode="w",
#     format='%(levelname)s::%(asctime)s %(message)s',
#     datefmt='%I:%M:%S')

"""number of threads can be passed as first cli argument"""
if len(sys.argv) > 1:
    n_processes = int(sys.argv[1])
else:
    n_processes = 4

"""
setting filepath to look for sim results. This is setup so that it works on the itb computers.
"""


def to_rgb(h):
    return [int(h[2 * i:2 * i + 2], 16) / 255 for i in range(3)]


all_refine_paths = []

path_list = glob.glob(os.path.join(
    "/".join(path.split("/")[:-2])
    , scan_name + "*"))

for path in path_list:
    refine_paths = glob.glob(os.path.join(path, "refine_*"))
    refine_paths = [re.search(r"refine_(?P<n>\d)_.*", rp) for rp in refine_paths]
    refine_paths = [rp.string if (int(rp.groupdict()["n"]) in [0, 1, 2, 3, 4, 5]) else None for rp in refine_paths]
    all_refine_paths = all_refine_paths + refine_paths


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def run(path):
    if path is not None:
        os.system("python post_process.py {n} {p}".format(p=path + "/", n=n_processes))


n = 1
chks = list(chunks(all_refine_paths, n))

for sublst in tqdm.tqdm(chks):

    if n == 1:
        run(sublst[0])
    else:
        with Pool(n, maxtasksperchild=1) as pool:
            pool.map(run, sublst)
