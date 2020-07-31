import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

get_dataframes = []#[[]]*3
saving_dataframe = pd.DataFrame()

from parameters_q_fraction import path
# get_dataframes.append([])
# # path = "/extra/brunner/10x10x10/R_lognorm/run" + str(j) + "/"
# # path = "/extra/brunner/10x10x10/q_fraction_exact/run" + str(j) + "/"
# user = getpass.getuser()
# path = "/extra/brunner/thesis/kinetic/q_fraction_wave_gamma_scan/"
# # ext_cache="/extra/brunner/para_handling/kinetic/R_lognorm_ext_cache/"

T = np.arange(0, 200, 1)
dt = 3600

global_df = pd.read_hdf(path + 'global_df.h5', mode="r")
cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")
print(global_df["Concentration"]*1e3)