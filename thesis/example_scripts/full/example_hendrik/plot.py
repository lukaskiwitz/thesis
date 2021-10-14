import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from parameters import path

"""
post_process saves the data in pandas dataframes. These can then be used to visualise the results.
seaborn is very convenient for plotting directly from dataframes.
global_df contains "global" results. These are calculated from cell_df and typically contain averages and statistics.
"""
f = lambda df: df.loc[df.model_name == "pde_model"]

global_df = pd.read_hdf(path + 'global_df.h5', mode="r")
global_df = f(global_df)
cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")
cell_df = f(cell_df)

sns.set(rc={'figure.figsize': (14, 6)})
sns.set_context("talk", font_scale=1, rc={"lines.linewidth": 2.5})
sns.set_style("ticks")
fig = plt.figure()
plt.subplots_adjust(wspace=.3)

"""
IL-2_surf_c is the IL-2 surface concentration on each cell. Here seaborn averaged over all the cells.
"""
fig.add_subplot(1,2, 1)
ax1 = sns.lineplot(x="fractions_sec", y="IL-2_surf_c", data=cell_df)
ax1.set(xlabel="Fraction of secreting cells", ylabel="Surface concentration [nM]", title="cell_df")

"""
Same can be done with the globel_df's "surf_c", just without the standard deviation:
"""
fig.add_subplot(1,2, 2)
ax2 = sns.lineplot(x="fractions_sec", y="surf_c", data=global_df)
ax2.set(xlabel="Fraction of secreting cells", ylabel="", title="global_df")
plt.tight_layout()
plt.show()

"""
seaborn can also show every single cells surface concentration
"""
sns.set(rc={'figure.figsize': (7, 6)})
sns.set_context("talk", font_scale=1, rc={"lines.linewidth": 0.1})
sns.set_style("ticks")
fig = plt.figure()
plt.subplots_adjust(wspace=.3)
fig = plt.figure()
ax3 = sns.lineplot(x="fractions_sec", y="IL-2_surf_c", data=cell_df, estimator=None, units="id_id")
ax3.set(xlabel="Fraction of secreting cells", ylabel="Each cells surface concentration [nM]")
plt.tight_layout()
plt.show()

