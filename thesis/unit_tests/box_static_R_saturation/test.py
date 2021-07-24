from  thesis.main.MyPlotter import Plotter
from parameters import path
import numpy as np

plotter = Plotter(path)
gdf = plotter.global_df
cdf = plotter.cell_df
from parameters import geometry

perfect_volume = (geometry["x_grid"]**3)- (4/3* np.pi * geometry["rho"]**3 * len(cdf["id"].unique()))


def run_successful(gdf):
    return len(np.where(gdf.success == False)[0]) == 0

def concentration_bounded(gdf):
    return len(np.where(gdf.Concentration > 1)[0]) == 0

def surface_concentration_bounded(cdf):
    return len(np.where(cdf["IL-2_surf_c"] > 1)[0]) == 0

def mesh_volume(gdf):
    return len(np.where(np.abs(gdf.MeshVolume - perfect_volume) > perfect_volume*0.1)[0]) == 0

assert run_successful(gdf)
assert concentration_bounded(gdf)
assert surface_concentration_bounded(cdf)
assert mesh_volume(gdf)



