import logging
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import poisson

from parameters import cytokines, cell_types_dict, geometry, numeric, boundary
from thesis.main.InternalSolver import InternalSolver
from thesis.main.PostProcess import PostProcessor
from thesis.main.StateManager import StateManager
from thesis.scenarios.box_grid import setup, assign_fractions

"""Define paths for result folder structure"""
solution_path = "/extra/kiwitz/cell_behavior_example/test_1/"
os.makedirs(solution_path, exist_ok=True)


class SimpleThresholdSolver(InternalSolver):
    name = "SimpleThresholdSolver"

    def on_type_change(self, p, replicat_index, entity=None):
        pass

    def step(self, t1, t2, dt, p, entity=None, **kwargs):

        il2_threshold = p.get_physical_parameter("ths", "IL-2").get_in_post_unit()
        il2 = p.get_physical_parameter("surf_c", "IL-2").get_in_post_unit()

        if entity.type_name == "default" and np.random.uniform(0, 1) > poisson.cdf(k=1, mu=0.25):
            if 2 * il2_threshold < il2:
                entity.change_type = "abs"
            elif 0.5 * il2_threshold > il2:
                entity.change_type = "sec"
        return p


"""Simulation"""
scenario = setup(
    cytokines,
    cell_types_dict,
    boundary,
    geometry,
    numeric)

stMan = StateManager(solution_path)
scenario.internal_solvers.append(SimpleThresholdSolver)
stMan.scenario = scenario
t_unit = 3600
stMan.T = np.arange(0, 20 * t_unit, t_unit / 4)


def pre_scan(state_manager, scan_index):
    assign_fractions(state_manager.sim_container, 0)


stMan.pre_scan = pre_scan
stMan.run()

"""Post Processing"""
pp = PostProcessor(solution_path)
pp.unit_length_exponent = -6
pp.run_post_process(4)

"""Plotting"""
global_df = pd.read_hdf(os.path.join(solution_path, "global_df.h5"))
cell_df = pd.read_hdf(os.path.join(solution_path, "cell_df.h5"))
global_df.time = global_df.time / 3600
cell_df.time = cell_df.time / 3600

fig, ax = plt.subplots(2, 1)

sns.lineplot(x="time", y="Concentration", hue="model_name", data=global_df, ax=ax[0])
sns.lineplot(x="time", y="id", hue="model_name", style="type_name",
             data=cell_df.groupby(
                 ["time", "time_index", "model_index", "model_name", "type_name"]).count().reset_index(),
             ax=ax[1]
             )

ax[0].set_xlabel("time (h)")
ax[0].set_ylabel("Mean IL-2")
ax[1].set_xlabel("time (h)")
ax[1].set_ylabel("Number of cells")
plt.tight_layout()
plt.savefig(os.path.join(solution_path, "plot.pdf"))
plt.show()
