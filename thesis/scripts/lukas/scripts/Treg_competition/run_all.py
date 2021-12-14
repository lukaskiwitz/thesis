import os
from multiprocessing import Pool


def run(command):
    os.system(command)


runs = [

    "python run.py 1 1 860 1",  # cs ec50 treg_k treg_halftime
    "python run.py 1 3 860 1",  # cs ec50 treg_k treg_halftime
    "python run.py 1 1 215 1",  # cs ec50 treg_k treg_halftime
    "python run.py 1 1 860 0.25",  # cs ec50 treg_k treg_halftime
    "python run.py 1 3 215 0.25",  # cs ec50 treg_k treg_halftime

    "python run.py 0.5 1 860 1",  # cs ec50 treg_k treg_halftime
    "python run.py 0.5 3 860 1",  # cs ec50 treg_k treg_halftime
    "python run.py 0.5 1 215 1",  # cs ec50 treg_k treg_halftime
    "python run.py 0.5 1 860 0.25",  # cs ec50 treg_k treg_halftime
    "python run.py 0.5 3 215 0.25",  # cs ec50 treg_k treg_halftime

    "python run.py 0 1 860 1",  # cs ec50 treg_k treg_halftime
    "python run.py 0 3 860 1",  # cs ec50 treg_k treg_halftime
    "python run.py 0 1 215 1",  # cs ec50 treg_k treg_halftime
    "python run.py 0 1 860 0.25",  # cs ec50 treg_k treg_halftime
    "python run.py 0 3 215 0.25",  # cs ec50 treg_k treg_halftime
]

with Pool(4, maxtasksperchild=1) as pool:
    pool.map(run, runs)
