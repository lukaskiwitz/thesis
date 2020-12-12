try:
    import fenics as fcs
except RuntimeError:
    import os
    os.environ['PATH'] = '/home/brunner/anaconda3/envs/Lukas2/bin:/home/brunner/.local/bin:/home/brunner/anaconda3/condabin:/usr/local/bin:/usr/bin:/bin:/usr/local/games:/usr/games:/opt/puppetlabs/bin'
    import fenics as fcs

import dolfin as dlf
import numpy as np
from scipy import spatial
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from mpi4py import MPI
from scipy.integrate import solve_ivp
from scipy.constants import N_A
comm = MPI.COMM_WORLD
rank = comm.rank

save_path = r"/home/brunner/thesis/thesis/scripts/patrick/yukawa_analytical/test/"

########################################################################################################################
ext_path = r"/extra/brunner/thesis/kinetic/standard/3D_ext_cache/"
mesh_path = ext_path + "mesh.xdmf"
markers_path = ext_path + "boundary_markers.h5"

df_path = "/extra/brunner/thesis/kinetic/standard/8_big_scan/"

print("replicating setup from", df_path)
global_df = pd.read_hdf(df_path + 'global_df.h5', mode="r")
cell_df = pd.read_hdf(df_path + 'cell_df.h5', mode="r")
########################################################################################################################

starting_surf_c = cell_df.loc[(cell_df["time_index"] == 1) & (cell_df["type_name"] == "Tsec"), "IL-2_surf_c"].mean()
time_range = np.sort(cell_df["time"].unique())
dt = 1
frac = 0.25
gammas = np.sort(cell_df["IL-2_gamma"].unique())

########################################################################################################################

cell_df_split = cell_df.loc[(cell_df["time"] == 0.0) & (cell_df["IL-2_gamma"] == 1/100) & (cell_df["IL-2_Tsec_fraction"] == frac)].copy()
cell_points = np.array([cell_df_split["x"].values, cell_df_split["y"].values, cell_df_split["z"].values]).T
cell_points_id_map = {}
for e, entry in enumerate(cell_df_split["id"].values):
    cell_points_id_map[entry] = e


mesh = dlf.Mesh()
with dlf.XDMFFile(dlf.MPI.comm_world, mesh_path) as f:
    f.read(mesh)
print("mesh loaded from", mesh_path)

# get cells boundary vertices
# read in boundary markers
boundary_markers = fcs.MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
with fcs.HDF5File(fcs.MPI.comm_world, markers_path, "r") as f:
    f.read(boundary_markers, "/boundaries")
# determine which boundary marker id correspondes to which vertices
patch_indices = np.unique(boundary_markers.array())
pi_to_verts = {}
log = fcs.get_log_level()
fcs.set_log_level(50)
for patch_index in patch_indices:
    facets = np.array([i for i in fcs.SubsetIterator(boundary_markers, patch_index)])
    verts = np.unique(np.array(list(map(lambda x: x.entities(0), facets))))
    if len(verts) < 100:
        pi_to_verts[patch_index] = verts
fcs.set_log_level(log)
# determine which boundary marker id corresponds to which cell id
cell_id_vertices = {}
tree = spatial.KDTree(cell_points)
for key,values in pi_to_verts.items():
    mean_values = np.mean(mesh.coordinates()[values], axis=0)
    cell_id_vertices[tree.query(mean_values)[1]] = values

print("boundary markers loaded")


########################################################################################################################
def hill(c, k, N, min, max):
    return (min * k ** N + max * c ** N) / (k ** N + c ** N)



base_path = '/home/brunner/thesis/thesis/scripts/patrick/yukawa_analytical/test/'
path = base_path
cell_df = pd.read_hdf(path + 'cell_df.h5', mode="r")
for id in cell_df_split.loc[cell_df_split["type_name"] == "Tsec", "id"]:
    cell_df.loc[cell_df["id"] == id, "type_name"] = "Tsec"


cell_df = pd.DataFrame()
global_df = pd.DataFrame()

for g,gamma in enumerate(gammas):
    xdmffile = fcs.XDMFFile(save_path + "test_" + str(gamma) + ".xdmf")
    V = fcs.FunctionSpace(mesh, "P", 1)
    v3 = fcs.Function(V)

    avg_c = starting_surf_c
    N = 3
    rho = 5

    min_x = cell_points[0][0] - rho
    max_x = cell_points[-1][0] + rho
    width = max_x - min_x
    volume = width**3
    Treg_density = len(cell_df_split.loc[cell_df_split["type_name"] == "Treg"])/volume

    # for Receptor feedback
    D = 100 #um^2/s
    kon = 100 # 1/s
    R = 1e4 #Receptors
    u = avg_c
    Kc = 1/(kon*1e3)
    Amax = 1 # molecule/s
    eta = 1/72000 # 1/s

    # for stimulus scaling
    niche_base = np.sqrt(D/(kon*R*avg_c**N/(Kc**N+avg_c**N)))
    K = avg_c

    Tsec_surf_c_list = np.zeros((len(time_range),len(cell_df_split.loc[cell_df_split["type_name"] == "Tsec", "id"])))
    surf_c_list = np.zeros((len(time_range),len(cell_df_split["id"].values)))

    niche_list = []
    mean_c_list = []
    R_list = [R]

    for t,time in enumerate(time_range[1:]):
        dt = time - time_range[t]
        # niche = np.sqrt(D * Kc / ((u/(Kc + u))*Treg_density) * 1e-9 * N_A / 1e15) #np.sqrt(D*u*(u/(Kc+u))/(kon*R*Kc*Treg_density) * 1e-9*N_A/1e15) # last factor: nmol/l -> molecules/um^3
        niche = np.sqrt(D / (kon * R * Kc * (u / (Kc + u)) * Treg_density))
        niche_list.append(niche)
        # print(niche)

        stimuli_list = []
        added_stimuli = None
        for id in cell_df_split.loc[cell_df_split["type_name"] == "Tsec", "id"]:
            stimuli_list.append('(%f * %f)/std::sqrt(pow(x[0] - %f, 2) + pow(x[1] - %f, 2) + pow(x[2] - %f, 2)) * exp(-(std::sqrt(pow(x[0] - %f, 2) + pow(x[1] - %f, 2) + pow(x[2] - %f, 2)) - %f)/%f)' \
                                 %(K, rho, cell_points[cell_points_id_map[id]][0], cell_points[cell_points_id_map[id]][1], cell_points[cell_points_id_map[id]][2], cell_points[cell_points_id_map[id]][0],cell_points[cell_points_id_map[id]][1],cell_points[cell_points_id_map[id]][2], rho, niche))
            # if added_stimuli == None:
            #     added_stimuli = '(%f * %f)/std::sqrt(pow(x[0] - %f, 2) + pow(x[1] - %f, 2) + pow(x[2] - %f, 2)) * exp(-(std::sqrt(pow(x[0] - %f, 2) + pow(x[1] - %f, 2) + pow(x[2] - %f, 2)) - %f)/%f)' \
            #                     %(K, rho, cell_points[cell_points_id_map[id]][0], cell_points[cell_points_id_map[id]][1], cell_points[cell_points_id_map[id]][2], cell_points[cell_points_id_map[id]][0],cell_points[cell_points_id_map[id]][1],cell_points[cell_points_id_map[id]][2], rho, niche)
            # else:
            #    added_stimuli = added_stimuli + " + " + '(%f * %f)/std::sqrt(pow(x[0] - %f, 2) + pow(x[1] - %f, 2) + pow(x[2] - %f, 2)) * exp(-(std::sqrt(pow(x[0] - %f, 2) + pow(x[1] - %f, 2) + pow(x[2] - %f, 2)) - %f)/%f)' \
            #                    %(K, rho, cell_points[cell_points_id_map[id]][0], cell_points[cell_points_id_map[id]][1], cell_points[cell_points_id_map[id]][2], cell_points[cell_points_id_map[id]][0],cell_points[cell_points_id_map[id]][1],cell_points[cell_points_id_map[id]][2], rho, niche)
        added_stimuli = "+".join(stimuli_list)
        exit()
        kwargs = {'element' : V.ufl_element()}
        exp = fcs.Expression(added_stimuli, degree=5)
        v3.interpolate(exp)
        test = fcs.project(exp,V)

        for id in cell_df_split["id"].values:
            id = int(id)
            tmp_list = []
            for entry in mesh.coordinates()[cell_id_vertices[id]]:
                tmp_list.append(v3(entry))
            mean_surf_c = np.mean(tmp_list)
            surf_c_list[t+1][id] = mean_surf_c

        for i,id in enumerate(cell_df_split.loc[cell_df_split["type_name"] == "Tsec", "id"]):
            id = int(id)
            tmp_list = []
            for entry in mesh.coordinates()[cell_id_vertices[id]]:
                tmp_list.append(v3(entry))
            Tsec_surf_c_list[t+1][i] = np.mean(tmp_list)

        mean_c = np.mean(v3.compute_vertex_values())
        mean_c_list.append(mean_c)

        ########## Receptor feedback #############
        R_mean = 1e4
        kmin = R_mean / gamma * eta
        kmax = R_mean * gamma * eta

        try:
            # k = solve((kmin * k_x ** N + kmax * c_0 ** N) / ((k_x ** N + c_0 ** N) * eta) - R_mean)[0]
            k = ((avg_c ** N * (eta * R_mean - kmax)) / (kmin - eta * R_mean)) ** (1 / N)
        except:  # if gamma ~ 1
            k = 1

        dict = {"N": N, "gamma": gamma, "eta": eta, "c_0": avg_c, "il2": mean_c, "R_il2": R, "R_mean": R_mean,
                "kmin": kmin, "kmax": kmax, "k": k}


        def func(t, y, p, dummy=0):
            R_il2 = y[0]
            dR_dt = ((p["kmin"] * p["k"] ** p["N"] + p["kmax"] * p["il2"] ** p["N"]) / (
                    p["k"] ** p["N"] + p["il2"] ** p["N"]) - p["eta"] * R_il2)
            return [dR_dt]


        y0 = [dict["R_il2"]]
        t_space = np.linspace(time_range[t], time_range[t] + dt, 100)
        result = solve_ivp(func, t_span=[time_range[t], time_range[t] + dt], t_eval=t_space, y0=y0, method="LSODA", args=(dict, 0))
        R = result.y[0][-1]
        R_list.append(R)

        print("time =", time)
        xdmffile.write(v3,t)
    xdmffile.close()

    ########################################################################################################################
    if g == 0:
        ids = np.concatenate([cell_df_split["id"].values for x in range(len(time_range))])
        times = np.array([[x] * len(cell_df_split["id"].values) for x in time_range]).flatten()
        time_indices = np.array([[x] * len(cell_df_split["id"].values) for x in range(len(time_range))]).flatten()
    cells_R = np.array([[x] * len(cell_df_split["id"].values) for x in R_list]).flatten()
    xs = np.concatenate([cell_df_split["x"].values for x in range(len(time_range))])
    ys = np.concatenate([cell_df_split["y"].values for x in range(len(time_range))])
    zs = np.concatenate([cell_df_split["z"].values for x in range(len(time_range))])

    std_list = []
    for t,time in enumerate(time_range[1:]):
        t += 1
        std_list.append(surf_c_list[t].std())

    #dict with which names correnspond to which list
    d_cell_df = {
        "IL-2_surf_c": surf_c_list.flatten(),
        "id": ids,
        "time" : times,
        "time_index": time_indices,
        "IL-2_R": cells_R,
        "x": xs,
        "y": ys,
        "z": zs,
        "type_name": "Treg",
        "IL-2_Tsec_fraction": frac,
        "IL-2_gamma": gamma,
    }

    d_global_df = {
        "time": time_range[1:],
        "surf_c": np.mean(surf_c_list.flatten()),
        "Concentration": mean_c_list,
        "no_of_vertices": len(v3.vector()),
        "IL-2_R": R_list[1:],
        "niche_size": niche_list,
        "IL-2_gamma": gamma,
        "surf_c_std": std_list
    }

    #fill df's
    tmp_cell_df = pd.DataFrame()
    for entry in d_cell_df.items():
        if g == 0:
            cell_df[entry[0]] = entry[1]
        else:
            tmp_cell_df[entry[0]] = entry[1]
    cell_df = cell_df.append(tmp_cell_df, ignore_index=True)

    tmp_global_df = pd.DataFrame()
    for entry in d_global_df.items():
        if g == 0:
            global_df[entry[0]] = entry[1]
        else:
            tmp_global_df[entry[0]] = entry[1]
    global_df = global_df.append(tmp_global_df, ignore_index=True)


global_df.to_hdf(save_path + "global_df.h5", "w")
cell_df.to_hdf(save_path + "cell_df.h5", "w")

global_df["time"] /= 3600


sns.set(rc={'figure.figsize':(11,8.27)})
fig = plt.figure()
plt.subplots_adjust(wspace=.3)

# subplot_style
a_x = 2
a_y = 2


my_hue = "IL-2_gamma" #"D"
x_axis = "time"
averaging_values = np.sort(global_df[x_axis].unique())

plot_D = 0
stoppingPoint = None
startingPoint = None

yscale = "linear"
xscale = "linear"

yscaleR = "log"

# ylim = (1.0, 1.55) #(0.01, 0.045)
ylim = (None, None)
ylimR = (None, None) #(1, 23000)

# hue_order = [1/10]
# hue_order = ["$%s$" % x for x in hue_order]
hue_order = None


fig.add_subplot(a_x,a_y, 1)
global_df["surf_c"] = global_df["surf_c"]
#sns.lineplot(x="fraction", y="mean_surf_c_il2", data=saving_dataframe,hue=group_variables[1])
ax3  = sns.lineplot(x=x_axis, y="surf_c", data=global_df.sort_values(by="time")[startingPoint:stoppingPoint],legend=False, color="black")
# ax3  = sns.lineplot(x="t", y="surf_c_il2", data=cell_df.loc[cell_df["id"] == 1][:stoppingPoint],hue=my_hue, hue_order=hue_order, color="black")
#ax1.errorbar(np.linspace(0.00005,1.0,20), plotting_df[5:]["surf_c"], yerr=plotting_df[5:]["sd"], fmt='-o')
ax3.set(xlabel="time (h)", ylabel="mean surf.c. in nM", title="surface c.", yscale=yscale, xscale=xscale, ylim=ylim)#, xticklabels=[]) #, ylim=(0.275,0.45))
# ax3.set_xticklabels([])

fig.add_subplot(a_x,a_y, 2)
ax6 = sns.lineplot(x=x_axis, y="surf_c_std", data=global_df.sort_values(by="time")[startingPoint:stoppingPoint], label="surf_c_std", legend=False, color="black")
ax6.set(xlabel="time (h)", ylabel="surf.c. std. in nM", title="surf. std", xticklabels=[], yscale=yscale, xscale=xscale, ylim=ylim)
# ax6.set_xticklabels([])
#
#
fig.add_subplot(a_x,a_y, 3)
global_df["surf_c_std_norm"] = global_df["surf_c_std"]/global_df["surf_c"]
ax7 = sns.lineplot(x=x_axis, y="surf_c_std_norm", data=global_df.sort_values(by="time")[startingPoint:stoppingPoint], label="surf_c_std", legend=False, color="black")
ax7.set(xlabel="time (h)", ylabel="surf.c. std./mean surf.c.", title="normed surf. std", yscale=yscale, xscale=xscale, ylim=(None, None))
# ax7.set_xticklabels([])


fig.add_subplot(a_x, a_y, 4)
ax8 = sns.lineplot(x=x_axis, y="IL-2_R", data=global_df.sort_values(by="time")[startingPoint:stoppingPoint], label="R_mean", legend=False, color="black")
ax8.set(xlabel="time (h)", ylabel="Receptors", title="R mean", xticklabels=[], xscale="linear", yscale=yscaleR, ylim=ylimR)
# ax8.set_xticklabels([])

plt.show()

