#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 12:56:59 2019

@author: Lukas Kiwitz
"""
import multiprocessing as mp
import os
import random
from typing import List, Dict

import KDEpy
import dolfin as dlf
import fenics as fcs
import lxml.etree as et
import matplotlib.pyplot as plt
import mpi4py.MPI as MPI
import numpy as np
import pandas as pd

import thesis.main.StateManager as st
from thesis.main.ParameterSet import ParameterSet
from thesis.main.myDictSorting import groupByKey
from thesis.main.my_debug import message


class ComputeSettings:
    """
    Class to store computation input and output

    :param file_path: file path to volume file
    :type file_path: str

    :param field: field name in h5 file
    :type field: str

    :param cell_data: cell data
    :type cell_data: List[Dict]

    :param dynamic: dynamic parameters
    :type dynamic: Dict

    :param scan_index: scan index
    :type scan_index: int

    :param time_index: time index
    :type time_index: float
    """

    # def __init__(self, file_path: str, field: str, mesh: fcs.Mesh, u: fcs.Function, cell_data: List[Dict],
    #              boundary_markers: fcs.MeshFunction, dynamic: Dict, scan_index: int, time_index: float) -> None:
    def __init__(self) -> None:
        self.file_path: str = ""
        """file path to volume file
        """
        self.field: str = ""
        self.dynamic: Dict = {}
        self.scan_index: int = 0
        self.time_index: float = 0
        self.tmp_path = ""

    def set_mesh(self, mesh: fcs.Mesh):
        self._mesh = mesh

    def set_u(self, u: fcs.Function):
        self._u: fcs.Function = u

    # noinspection PyPep8Naming,PyPep8Naming
    def set_V(self, V: fcs.FunctionSpace):
        self._V: fcs.FunctionSpace = V

    def set_boundary_markers(self, bm: fcs.MeshFunction):
        self._boundary_markers: fcs.MeshFunction = bm

    def get_mesh(self) -> fcs.Mesh:
        return self._mesh

    def get_u(self) -> fcs.Function:
        return self._u

    # noinspection PyPep8Naming
    def get_V(self) -> fcs.FunctionSpace:
        return self._V

    def get_boundary_markers(self) -> fcs.MeshFunction:
        return self._boundary_markers


class PostProcessor:

    def __init__(self, path: str) -> None:
        self.cellDump = []
        self.out_tree_path = path + "postProcess.xml"
        self.path = path
        self.cell_dataframe: pd.DataFrame = pd.DataFrame()
        self.global_dataframe: pd.DataFrame = pd.DataFrame()
        self.cell_stats: pd.DataFrame = pd.DataFrame()
        self.unit_length_exponent: int = 1
        self.computations = [c(), grad(), SD()]

        self.image_settings = {
            "cell_colors": "Dark2",
            "cell_color_key": "type_name",
            "legend_title": "",
            "dpi": 350,
        }

    def write_post_process_xml(self, threads, debug=False, make_images=False):
        """
        runs compute-function for all scans an writes result to xml file

        """

        assert type(threads) == int

        # initializes state manager from scan log
        self.stateManager = st.StateManager(self.path)
        self.stateManager.load_xml()

        tmp_path = self.path + "tmp/"
        os.makedirs(tmp_path, exist_ok=True)

        scatter_list: List[ComputeSettings] = []

        # loads timestep logs
        for scan_index, scan_sample in enumerate(self.stateManager.element_tree.findall("ScanContainer/ScanSample")):
            parameters = scan_sample.find("Parameters")

            for field_step in scan_sample.findall("TimeSeries/Step/Fields/Field"):
                compute_settings: ComputeSettings = ComputeSettings()
                compute_settings.path = self.path
                compute_settings.file_path = field_step.get("dist_plot_path")
                compute_settings.field = field_step.get("field_name")
                compute_settings.mesh_path = field_step.get("mesh_path")

                compute_settings.parameters = et.tostring(parameters)
                compute_settings.scan_index = int(scan_index)
                compute_settings.time_index = int(field_step.getparent().getparent().get("time_index"))
                compute_settings.time = field_step.getparent().getparent().get("time")
                compute_settings.make_images = make_images
                compute_settings.tmp_path = tmp_path
                compute_settings.unit_length_exponent = self.unit_length_exponent

                for k, v in self.image_settings.items():
                    compute_settings.__setattr__(k, v)

                compute_settings.computations = self.computations

                scatter_list.append(compute_settings)

        mesh_path = self.path + "/" + scatter_list[0].mesh_path
        threads = threads if threads <= os.cpu_count() else os.cpu_count()
        message("distributing to {threads} threads".format(threads=threads))

        with mp.Pool(processes=threads, initializer=initialise(mesh_path, self.cell_dataframe)) as p:
            result_list = p.map(compute, scatter_list)

        element_list = []
        for file in result_list:
            with open(file, "rb") as f:
                element_list.append(et.fromstring(f.read()))

        indexed_list = [
            {"scan_index": i.get("scan_index"),
             "time_index": i.get("time_index"),
             "time": i.get("time"),
             "entry": i}
            for i in element_list]

        indexed_list = groupByKey(indexed_list, ["scan_index"])
        for i, e in enumerate(indexed_list):
            indexed_list[i] = groupByKey(e, ["time_index"])

        post_process_result = et.Element("PostProcess")
        tree = et.ElementTree(element=post_process_result)
        for s in indexed_list:
            scan = et.SubElement(post_process_result, "Scan")
            #            message(s[0])
            scan.set("i", str(s[0][0]["scan_index"]))
            for t in s:
                for i in t:
                    time = et.SubElement(scan, "Step")
                    time.set("i", i["time_index"])
                    time.set("time", i["time"])
                    time.append(i["entry"])
        message("writing post process output to {p}".format(p=self.out_tree_path))
        tree.write(self.out_tree_path, pretty_print=True)

    def get_global_dataframe(self) -> pd.DataFrame:
        message("collecting global dataframe from post_process.xml")
        in_tree = et.parse(self.out_tree_path)

        result = []
        if self.cell_dataframe.empty:
            raise MyError.DataframeEmptyError("Post Processor Cell Dataframe")

        for scan in in_tree.findall("Scan"):

            for step in scan.findall("Step"):

                for file in step.findall("File"):
                    g_values = file.find("./GlobalResults")
                    scan_index = int(file.get("scan_index"))
                    time_index = int(file.get("time_index"))
                    time = float(file.get("time"))
                    field_name = file.get("field_name")
                    filter = lambda x: (x["time_index"] == time_index) & \
                                       (x["scan_index"] == scan_index)

                    d = {
                        "scan_index": scan_index,
                        "time_index": time_index,
                        "time": time,
                        "field_name": field_name,
                        "surf_c": self.cell_dataframe.loc[filter(self.cell_dataframe)][
                            "{field}_surf_c".format(field=field_name)].mean(),
                        "surf_c_std":self.cell_dataframe.loc[filter(self.cell_dataframe)][
                            "{field}_surf_c".format(field=field_name)].std()
                    }
                    parameter_set = ParameterSet("dummy_set", [])
                    parameter_set.deserialize_from_xml(file.find("Parameters/ParameterSet[@name='dynamic']"))
                    d.update(parameter_set.get_as_dictionary())

                    for v in g_values:
                        d.update({v.tag: float(v.text)})
                    result.append(d)

        return pd.DataFrame(result)

    def get_kde_estimators(self, n, ts, time_index, type_names):

        kernels = {}
        message("computing kde for time series: {n} and timestep {t}".format(n=n, t=time_index))
        for type_name in type_names:
            inital_cells = ts.loc[(ts["time_index"] == time_index) & (ts["type_name"] == type_name)]
            if inital_cells.shape[0] == 0:
                break
            elif inital_cells.shape[0] == 1:
                data = np.array([inital_cells["x"].iloc[0], inital_cells["y"].iloc[0], inital_cells["z"].iloc[0]])
                kernel = KDEpy.NaiveKDE("tri", bw=10).fit(data)
            else:
                data = np.array([inital_cells["x"], inital_cells["y"], inital_cells["z"]]).T
                kernel = KDEpy.TreeKDE("tri", bw=10).fit(data)

            # kernel = KDEpy.TreeKDE(bw='ISJ').fit(data)

            # grid_points = 100
            # grid, points = kernel.evaluate(grid_points)
            # x, y, z = np.unique(grid[:, 0]), np.unique(grid[:, 1]), np.unique(grid[:, 2])
            # v = points.reshape(grid_points, grid_points, grid_points).T
            #
            # plt.title(type_name)
            # plt.contour(x,y,v[:,:,0])
            # sns.scatterplot(x="x",y="y",data=inital_cells)
            # plt.show()

            kernels[type_name] = kernel

        return kernels

    def get_cell_dataframe(self, kde=False):

        self.stateManager = st.StateManager(self.path)
        self.stateManager.load_xml()

        result: pd.DataFrame = self.stateManager.get_cell_ts_data_frame()

        result = result.groupby("id").apply(lambda x: x.ffill().bfill()).drop_duplicates()

        result["time_index"] = result["time_index"].apply(lambda x: int(x))
        result["scan_index"] = result["scan_index"].apply(lambda x: int(x))
        result["time"] = result["time"].apply(lambda x: float(x))


        """---------------------------"""

        if kde:
            message("running kernel density estimation")
            r_grouped = result.groupby(["scan_index"], as_index=False)
            kde_result = pd.DataFrame()
            for scan_index, ts in r_grouped:
                kernels = []
                for time_index, step in ts.groupby(["time_index"]):
                    kde_time_index = time_index - 1 if time_index > 0 else 0
                    kernels.append(
                        self.get_kde_estimators(scan_index, ts, kde_time_index, result["type_name"].unique()))

                for time_index, step in ts.groupby(["time_index"]):

                    for type_name, kernel in kernels[time_index].items():
                        positions = np.array([step["x"], step["y"], step["z"]]).T

                        scores = kernel.evaluate(positions).T
                        scores = pd.Series(scores)
                        scores.index = step.index

                        step.insert(step.shape[1], "{type_name}_score".format(type_name=type_name), scores)

                    for type_name, kernel in kernels[0].items():
                        positions = np.array([step["x"], step["y"], step["z"]]).T

                        scores = kernel.evaluate(positions).T
                        scores = pd.Series(scores)
                        scores.index = step.index

                        step.insert(step.shape[1], "{type_name}_score_init".format(type_name=type_name), scores)

                    kde_result = kde_result.append(step)
            result = self._normalize_cell_score(kde_result)

        return result.drop_duplicates()

    def _normalize_cell_score(self, x):

        x = x.groupby(["time_index", "scan_index"], as_index=False)

        result = pd.DataFrame()
        for group in x:
            i = group[0]
            group = group[1]

            # ids = x.loc[
            #     (x["type_name"] != group["type_name"][0]) &
            #     (x["time_index"] == 0)
            #     ]["id"].unique()

            no_init = group  # group.loc[~group["id"].isin(ids)]

            for old in pd.Series(group.columns).str.extract("(.*_score)").dropna()[0].unique():
                new = "{old}_norm".format(old=old)
                group.insert(group.shape[1], new, group[old] / float(no_init.mean()[old]))

            for old in pd.Series(group.columns).str.extract("(.*_score_init)").dropna()[0].unique():
                new = "{old}_norm".format(old=old)
                group.insert(group.shape[1], new, group[old] / float(no_init.mean()[old]))

            result = result.append(group)

        return result

    def get_stats(self):

        # temp_cell_df = self.cell_dataframe.drop(columns = ["x","y","z","id"])
        grouped = self.cell_dataframe.groupby(["type_name", "time_index", "scan_index"], as_index=False)

        des = grouped.describe(include='all')

        # des = des.reset_index()

        return des

    def run_post_process(self, threads, kde=False, make_images=False):
        self.cell_dataframe = self.get_cell_dataframe(kde=kde)
        self.write_post_process_xml(threads, make_images=make_images)
        self.global_dataframe = self.get_global_dataframe()


def get_concentration_conversion(unit_length_exponent: int):
    assert isinstance(unit_length_exponent, int)
    return 10 ** (-1 * (unit_length_exponent * 3 + 3))


def get_gradient_conversion(unit_length_exponent: int):
    assert isinstance(unit_length_exponent, int)
    exp = (-1 * (unit_length_exponent * 3 + 3)) - 1
    exp -= (6 + unit_length_exponent)  # to nM/um

    return 10 ** (exp)


class PostProcessComputation():
    add_xml_result = True

    def __init__(self):
        self.name = "PostProcessComputation"

    def __call__(self, u, grad, c_conv, grad_conv, mesh_volume, **kwargs):
        return 0


class SD(PostProcessComputation):

    def __init__(self):
        self.name = "SD"

    def __call__(self, u, grad, c_conv, grad_conv, mesh_volume, **kwargs):
        sd: float = np.std(np.array(u.vector())) * c_conv

        return sd


class c(PostProcessComputation):

    def __init__(self):
        self.name = "Concentration"

    def __call__(self, u, grad, c_conv, grad_conv, mesh_volume, **kwargs):
        nodal_values = np.array(u.vector())

        mean_c = np.sum(nodal_values) * (1 / nodal_values.shape[0])
        return mean_c * c_conv


class grad(PostProcessComputation):

    def __init__(self):
        self.name = "Gradient"

    def __call__(self, u, grad, c_conv, grad_conv, mesh_volume, **kwargs):
        g = np.reshape(grad.vector().vec().array, (u.vector().vec().size, 3))
        g = np.transpose(g)

        my_grad = np.sqrt(np.power(g[0], 2) + np.power(g[1], 2) + np.power(g[2], 2))
        my_grad = np.mean(my_grad) * grad_conv

        return my_grad


def get_mesh_volume(mesh):
    sum = 0
    for cell in fcs.cells(mesh):
        sum += cell.volume()
    return sum


def compute(compute_settings: ComputeSettings) -> str:
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    """
    performs computations according to compute_settings and return result as element tree string

    :param compute_settings:
    :return:
    """

    u: fcs.Function = fcs.Function(function_space)
    with fcs.HDF5File(comm, compute_settings.path + "/" + compute_settings.file_path, "r") as f:
        f.read(u, "/" + compute_settings.field)

    result: et.Element = et.Element("File")
    global_results: et.Element = et.SubElement(result, "GlobalResults")

    result.set("field_name", str(compute_settings.field))
    result.set("dist_plot_path", str(compute_settings.file_path))

    result.append(et.fromstring(compute_settings.parameters))

    scan_index = str(compute_settings.scan_index)
    result.set("scan_index", scan_index)
    time_index = str(compute_settings.time_index)
    result.set("time_index", time_index)
    result.set("time", str(compute_settings.time))

    message("running global computations for step {t} in scan {s}".format(t=time_index, s=scan_index))

    mesh_volume_element = et.SubElement(global_results, "MeshVolume")
    mesh_volume_element.text = str(mesh_volume)

    p = ParameterSet("dummy", [])
    p.deserialize_from_xml(et.fromstring(compute_settings.parameters))

    V_vec: fcs.VectorFunctionSpace = fcs.VectorFunctionSpace(mesh, "P", 1)

    grad: fcs.Function = fcs.project(fcs.grad(u), V_vec, solver_type="gmres")

    c_conv = get_concentration_conversion(compute_settings.unit_length_exponent)
    g_conv = get_gradient_conversion(compute_settings.unit_length_exponent)

    if compute_settings.make_images:
        try:
            make_images(u, c_conv, compute_settings)
        except RuntimeError as e:
            message("could not make images for scanindex {si} at t = {ti}".format(si = scan_index, ti = time_index))

    for comp in compute_settings.computations:

        assert isinstance(comp, PostProcessComputation)

        if comp.add_xml_result:
            result_element: et.Element = et.SubElement(global_results, comp.name)
            result_element.text = str(comp(
                u,
                grad,
                c_conv,
                g_conv,
                mesh_volume,
                V=function_space,
                V_vec=V_vec)
            )
        else:
            comp(
                u,
                grad,
                c_conv,
                g_conv,
                mesh_volume,
                V=function_space,
                V_vec=V_vec)
    tmp_path = compute_settings.tmp_path

    filename = tmp_path + "post_{r}.txt".format(r=str(random.randint(0, 2 ** 32)))
    while filename in os.listdir(tmp_path):
        filename = tmp_path + "post_{r}.txt".format(r=str(random.randint(0, 2 ** 32)))

    with open(filename, 'wb') as f:
        f.write(et.tostring(result))

    return filename


def initialise(mesh_path, _cell_dataframe):
    global mesh
    mesh = fcs.Mesh()
    with fcs.XDMFFile(mesh_path) as f:
        f.read(mesh)

    global mesh_volume
    mesh_volume = get_mesh_volume(mesh)

    global function_space
    function_space = fcs.FunctionSpace(mesh, "P", 1)

    global cell_dataframe
    cell_dataframe = _cell_dataframe


def get_color_dictionary(cell_df, cell_color_key, cell_colors, round_legend_labels=3):
    cell_color_key = cell_color_key

    if isinstance(cell_colors, Dict):
        color_dict = cell_colors
    else:
        try:
            cmap = plt.get_cmap(cell_colors)
        except:
            message("could not inteperet cell_colors as colormap")
            cmap = plt.get_cmap("Dark2")

        from numbers import Number

        if not isinstance(cell_df[cell_color_key].iloc[0], Number):
            categorical = True
            cell_types = cell_df[cell_color_key].unique()
        else:
            categorical = False
            cell_types = cell_df[cell_color_key].round(round_legend_labels).unique()

        color_dict = {}
        if cmap.N <= 8:
            for i, t in enumerate(cell_types):
                color_dict[t] = cmap(i)
        else:

            for i, t in enumerate(cell_types):
                i *= cmap.N / len(cell_types)
                color_dict[t] = cmap(int(i))

    legend_items = []

    if len(color_dict) > 8:
        step = int(len(color_dict) / 8)
        items = list(color_dict.items())[0::step]
        legend_color_dict = dict(items)
    else:
        legend_color_dict = color_dict

    for type, color in legend_color_dict.items():
        legend_items.append(
            plt.Line2D([0], [0], marker='o', color='w', label=type,
                       markerfacecolor=color, markersize=15, markeredgewidth=0, linewidth=0)
        )

    return color_dict, legend_items, categorical


def make_images(u, conv_factor, compute_settings):
    scan_index = compute_settings.scan_index
    time_index = compute_settings.time_index

    cell_color_key = compute_settings.cell_color_key
    legend_title = compute_settings.legend_title
    cell_colors = compute_settings.cell_colors
    round_legend_labels = compute_settings.round_legend_labels

    color_dict, legend_items, categorical = get_color_dictionary(cell_dataframe, cell_color_key, cell_colors,
                                                                 round_legend_labels=round_legend_labels)

    cell_df = cell_dataframe.loc[
        (cell_dataframe["scan_index"] == scan_index) &
        (cell_dataframe["time_index"] == time_index)
        ]
    u.set_allow_extrapolation(True)
    fig = plt.figure()
    for index, cell in cell_df.iterrows():
        p = [
            cell["x"],
            cell["y"]
        ]

        key = cell[cell_color_key] if categorical else round(cell[cell_color_key], round_legend_labels)
        if not key in color_dict:
            color = "black"
            message("no color provided for key {k}".format(k=key))
        else:
            color = color_dict[key]

        c = plt.Circle(p, 5, color=color)
        fig.gca().add_artist(c)

    coords = u.function_space().mesh().coordinates()
    coords_t = np.transpose(coords)
    x_limits = (
        np.min(coords_t[0]), np.max(coords_t[0])
    )
    y_limits = (
        np.min(coords_t[1]), np.max(coords_t[1])
    )

    rec_mesh = fcs.RectangleMesh(fcs.Point(x_limits[0], y_limits[0]), fcs.Point(x_limits[1], y_limits[1]), 200, 200)
    rev_V = fcs.FunctionSpace(rec_mesh, "P", 1)
    u_slice = fcs.interpolate(u, rev_V)

    mesh = u_slice.function_space().mesh()
    slice_v = u_slice.compute_vertex_values(mesh) * conv_factor

    triang = dlf.common.plotting.mesh2triang(mesh)

    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize

    # plt.tricontourf(triang, slice_v, levels=np.linspace(0, 0.1, 100))


    vmin = min(slice_v)
    vmax = max(slice_v)

    norm_cytokine = Normalize(vmin=vmin, vmax=vmax)

    mappable = ScalarMappable(norm=norm_cytokine, cmap="viridis")
    if hasattr(compute_settings,"colorbar_range"):
        tpl = compute_settings.colorbar_range
        mappable.set_clim(tpl[0], tpl[1])

    fig.gca().tricontourf(triang, slice_v, levels=100, norm = norm_cytokine)
    cb_cytokine = fig.colorbar(mappable)

    cb_cytokine.set_label("cytokine (nM)")

    fig.gca().set_xlabel(r"x ($\mu m $)")
    fig.gca().set_ylabel(r"y ($\mu m $)")

    if not categorical:
        vmin = min(color_dict.keys())
        vmax = max(color_dict.keys())

        norm = Normalize(vmin=vmin, vmax=vmax)
        cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cell_colors))
        cb.set_label(cell_color_key)
    else:
        plt.legend(handles=legend_items, title=legend_title)

    img_path = "images/scan_{scan_index}/{field}/".format(
        field=compute_settings.field,
        scan_index=scan_index
    )
    file_name = "{time_index}".format(
        time_index=time_index
    )
    os.makedirs(compute_settings.path + img_path, exist_ok=True)

    fig.savefig(compute_settings.path + img_path + file_name + ".png", dpi=compute_settings.dpi)
