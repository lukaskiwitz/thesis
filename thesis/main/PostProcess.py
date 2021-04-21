#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 12:56:59 2019

@author: Lukas Kiwitz
"""
import multiprocessing as mp
import os
import random
import json
from typing import List, Dict

import KDEpy
import dolfin as dlf
import fenics as fcs
import lxml.etree as et
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import mpi4py.MPI as MPI
import numpy as np
import pandas as pd

import thesis.main.StateManager as st
from thesis.main.ParameterSet import ParameterSet
from thesis.main.myDictSorting import groupByKey
from thesis.main.my_debug import message, warning
import thesis.main.MyError as MyError

from copy import deepcopy

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
        self.figure_width = 8.3/4
        self.pad_l = 0.3
        self.pad_r = 0.1
        self.render_paraview = False
        self.unit_name = "nM"
        self.marker_lookup = {}

        self.paraview_settings = {

            "cell_type_title_font_size": 20,
            "cell_type_label_font_size": 20,
            "slice_zoom": 1.2,
            "slice_origin": [0, 0, 40],
            "slice_normal": [0, 0, 1],
            "layer_distance": 20,
            "axis_title_font_size": 15,
            "axis_label_font_size": 15,
            "number_format": '%2.1g',
            "axis_ticks": True,
            "axis_edges": True,
            "field_color_preset": "Cool to Warm",
            "volume_camera_pos": [850,70,-110],
            "volume_raytracing": True,
            "volume_raytracing_progressive_passes": 0,
            "volume_raytracing_samples": 2,
            "volume_raytracing_ambient_samples": 2,
            "volume_raytracing_light_scale": 1,
            "marker_view_uniform_opacity":True,
            "lookup":{}

        }


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
        self.timing_dataframe: pd.DataFrame = pd.DataFrame()

        self.cell_stats: pd.DataFrame = pd.DataFrame()
        self.unit_length_exponent: int = 1
        self.computations = [c(), grad(), SD(), CV()]
        self.rc = {}

        self.image_settings = {
            "cell_colors": "Dark2",
            "cell_color_key": "type_name",
            "legend_title": "",
            "dpi": 350,
        }
        self.image_settings_fields = None
        self.visual_conversion = 1 #conversion factor from nano molar


    def write_post_process_xml(self, n_processes, debug=False, render_paraview = False, make_images=False, time_indices = None):
        """
        runs compute-function for all scans an writes result to xml file

        """

        plt.rcParams.update(self.rc)
        assert type(n_processes) == int

        # initializes state manager from scan log
        self.stateManager = st.StateManager(self.path)
        self.stateManager.load_xml()


        tmp_path = self.path + "tmp/"
        os.makedirs(tmp_path, exist_ok=True)

        scatter_list: List[ComputeSettings] = []


        # loads timestep logs

        for scan_index, scan_sample in enumerate(self.stateManager.element_tree.findall("ScanContainer/ScanSample")):
            parameters = scan_sample.find("Parameters")
            self.stateManager.rebuild_tree(scan_sample)
            for field_step in scan_sample.findall("TimeSeries/Step/Fields/Field"):

                markers = field_step.findall("../../Marker")
                compute_settings: ComputeSettings = ComputeSettings()
                compute_settings.additional_conversion = self.visual_conversion
                self.render_paraview = render_paraview
                compute_settings.markers = {e.get("marker_key"): e.get("path") for e in markers}

                marker_lookup = field_step.findall("../../MarkerLookupTable/MarkerLookup")
                compute_settings.marker_lookup = {e.get("key"): e.get("value") for e in marker_lookup}

                compute_settings.path = self.path
                compute_settings.file_path = field_step.get("dist_plot_path")
                compute_settings.solution_path = field_step.get("solution_path")

                compute_settings.field = field_step.get("field_name")
                compute_settings.mesh_path = field_step.get("mesh_path")
                compute_settings.dynamic_mesh = True if field_step.get("dynamic_mesh") == 'True' else False

                compute_settings.parameters = et.tostring(parameters)
                compute_settings.scan_index = int(scan_index)
                compute_settings.time_index = int(field_step.getparent().getparent().get("time_index"))
                if (not time_indices is None) and not (compute_settings.time_index in time_indices):
                    continue
                compute_settings.time = field_step.getparent().getparent().get("time")
                compute_settings.make_images = make_images
                compute_settings.render_paraview = render_paraview
                compute_settings.tmp_path = tmp_path
                compute_settings.unit_length_exponent = self.unit_length_exponent
                compute_settings.cell_df = self.cell_dataframe.loc[
                    (self.cell_dataframe["time_index"] == compute_settings.time_index)&
                    (self.cell_dataframe["scan_index"] == compute_settings.scan_index)
                ]


                image_settings = deepcopy(self.image_settings)

                if hasattr(self,"image_settings_fields") and isinstance(self.image_settings_fields,Dict):
                    for k,v in self.image_settings_fields.items():
                        if k == compute_settings.field:
                            if isinstance(v,Dict):
                                for kk,vv in v.items():
                                    image_settings[kk].update(vv)
                            else:
                                image_settings[k] = v


                for k, v in image_settings.items():
                    if isinstance(v,Dict):
                        if hasattr(compute_settings,k) and isinstance(compute_settings.__getattribute__(k),Dict):
                            compute_settings.__getattribute__(k).update(v)
                    else:
                        compute_settings.__setattr__(k, v)

                compute_settings.computations = self.computations

                scatter_list.append(compute_settings)

        mesh_prefix = self.path
        n_processes = n_processes if n_processes <= os.cpu_count() else os.cpu_count()
        message("distributing to {n_processes} processes".format(n_processes=n_processes))


        with mp.Pool(processes=n_processes, initializer=initialise(mesh_prefix,self.cell_dataframe)) as p:
            result_list = p.map(compute, scatter_list)


        # result_list  =[]
        # for s in scatter_list:
        #     initialise(mesh_prefix,self.cell_dataframe)
        #     result_list.append(compute(s))

        element_list = []
        for file in result_list:
            if not file is None:
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
                    if time_index == 0:
                        continue
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
                            "{field}_surf_c".format(field=field_name)].std(),
                        "surf_c_cv": self.cell_dataframe.loc[filter(self.cell_dataframe)][
                            "{field}_surf_c".format(field=field_name)].std()/self.cell_dataframe.loc[filter(self.cell_dataframe)][
                            "{field}_surf_c".format(field=field_name)].mean()

                    }
                    parameter_set = ParameterSet.deserialize_from_xml(file.find("Parameters/ParameterSet[@name='dynamic']"))

                    d.update(parameter_set.get_as_dictionary())

                    for v in g_values:
                        d.update({v.tag: float(v.text)})
                    result.append(d)

        return pd.DataFrame(result)

    def get_timing_dataframe(self) -> pd.DataFrame:

        records_path = os.path.join(self.path, "records")
        with open(os.path.join(records_path, "dump.json"), "r") as f:
            records = json.load(f)

        # df = pd.DataFrame(columns=["task","start","end","info"])

        result = []
        for group_name, task_group in records.items():

            for task in task_group:
                entry = {}
                entry["task"] = group_name
                entry["start"] = task["start"]
                entry["end"] = task["end"]

                for info_key, info_value in task["info"].items():
                    entry[info_key] = info_value

                result.append(entry)

        df = pd.DataFrame(result)

        offset = df["start"].min()
        df["start"] = df["start"].sub(offset)
        df["end"] = df["end"].sub(offset)
        df["duration"] = df["end"] - df["start"]
        df["name"] = df["task"].map(lambda x: x.split(":")[-1])


        return df

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

    def get_cell_dataframe(self, kde=False, time_indices = None, n_processes = 1):

        self.stateManager = st.StateManager(self.path)
        self.stateManager.load_xml()

        result: pd.DataFrame = self.stateManager.get_cell_ts_data_frame(time_indices = time_indices, n_processes = n_processes)

        try:
            result = result.groupby("id").apply(lambda x: x.ffill().bfill()).drop_duplicates()
        except:
            print("")

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


                    for type_name, kernel in kernels[time_index-1].items():
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


        # return result.drop_duplicates()

        return result

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

    def save_dataframes(self, extra_cell_constants = False):

        df = self.cell_dataframe
        ids = df["id_id"]
        if extra_cell_constants:
            cell_df_constant = df.loc[:, (df == df.iloc[0]).all()].iloc[0,:]
            cell_df = df.loc[:, (df != df.iloc[0]).any()]
            cell_df.to_hdf(os.path.join(self.path, "cell_df.h5"), key="df", mode="w")
            cell_df_constant.to_hdf(os.path.join(self.path, "cell_constants_df.h5"), key="df", mode="w")
        else:
            df.to_hdf(os.path.join(self.path, "cell_df.h5"), key="df", mode="w")

        self.global_dataframe.to_hdf(os.path.join(self.path, 'global_df.h5'), key="data", mode="w", data_columns=self.global_dataframe.columns)
        self.timing_dataframe.to_hdf(os.path.join(self.path,"timing_df.h5"), key="df", mode="w")

    def run_post_process(self, n_processes, render_paraview = False, extra_cell_constants = False, kde=False, make_images=False, time_indices = None):
        self.cell_dataframe = self.get_cell_dataframe(kde=kde, time_indices = time_indices, n_processes= n_processes)
        self.write_post_process_xml(n_processes, render_paraview = render_paraview, make_images=make_images, time_indices = time_indices)
        self.global_dataframe = self.get_global_dataframe()
        self.timing_dataframe = self.get_timing_dataframe()

        self.save_dataframes(extra_cell_constants = extra_cell_constants)


def get_concentration_conversion(unit_length_exponent: int):
    assert isinstance(unit_length_exponent, int)
    return 10 ** (-1 * (unit_length_exponent * 3 + 3))


def get_gradient_conversion(unit_length_exponent: int, target_exponent = -6):
    #converts gradient to nM/um for target_exponent = -6.
    assert isinstance(unit_length_exponent, int)
    # exp = (-1 * (unit_length_exponent * 3 + 3)) - 1
    # exp -= (6 + unit_length_exponent)  # to nM/um

    return 10 ** (-3-4 * unit_length_exponent + target_exponent)


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

class CV(PostProcessComputation):
    """Defines a custom computation"""

    def __init__(self):
        """this name will appear in output dataframe"""
        self.name = "CV"

    def __call__(self, u, grad, c_conv, grad_conv, mesh_volume, **kwargs):
        nodal_values = np.array(u.vector())
        sd: float = np.std(nodal_values) * c_conv
        mean_c = np.sum(nodal_values) * (1 / nodal_values.shape[0])*c_conv

        return sd/mean_c

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

    mesh = fcs.Mesh()
    with fcs.XDMFFile(os.path.join(mesh_prefix, compute_settings.mesh_path)) as f:
        f.read(mesh)

    mesh_volume = get_mesh_volume(mesh)

    function_space = fcs.FunctionSpace(mesh, "P", 1)

    u: fcs.Function = fcs.Function(function_space)
    with fcs.HDF5File(comm, os.path.join(compute_settings.path, compute_settings.file_path), "r") as f:
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

    p = ParameterSet.deserialize_from_xml(et.fromstring(compute_settings.parameters))

    V_vec: fcs.VectorFunctionSpace = fcs.VectorFunctionSpace(mesh, "P", 1)

    grad: fcs.Function = fcs.project(fcs.grad(u), V_vec, solver_type="gmres")

    c_conv = get_concentration_conversion(compute_settings.unit_length_exponent)
    g_conv = get_gradient_conversion(compute_settings.unit_length_exponent)

    if compute_settings.make_images:
        try:
            make_images(u, c_conv, compute_settings, visual_conv = compute_settings.additional_conversion)
        except RuntimeError as e:
            warning("could not make images for scanindex {si} at t = {ti}".format(si = scan_index, ti = time_index))


    if compute_settings.render_paraview:

        try:
            make_paraview_images(c_conv, compute_settings,  visual_conv = compute_settings.additional_conversion)
        except FileNotFoundError as e:
            warning("could not render paraview images for scanindex {si} at t = {ti}. pvbatch not found".format(si=scan_index, ti=time_index))
        except NameError as e:
            warning("could not render paraview images for scanindex {si} at t = {ti}. pvbatch not found".format(
                si=scan_index, ti=time_index))
        except AttributeError as e:
            warning("could not render paraview images for scanindex {si} at t = {ti}. Running data from old simulation?".format(
                si=scan_index, ti=time_index))


    for comp in compute_settings.computations:

        assert isinstance(comp, PostProcessComputation)
        try:
            if comp.add_xml_result:

                comp_result = comp(
                    u,
                    grad,
                    c_conv,
                    g_conv,
                    mesh_volume,
                    V=function_space,
                    V_vec=V_vec,
                    path = compute_settings.path,
                    solution_path = compute_settings.solution_path,
                    scan_index = compute_settings.scan_index,
                    time_index = compute_settings.time_index,
                    cell_df = compute_settings.cell_df,
                    p = p
                )

            else:
                comp(
                    u,
                    grad,
                    c_conv,
                    g_conv,
                    mesh_volume,
                    V=function_space,
                    V_vec=V_vec,
                    path=compute_settings.path,
                    scan_index=compute_settings.scan_index,
                    time_index=compute_settings.time_index,
                    cell_df=compute_settings.cell_df,
                    p=p)
        except Exception as e:
            message("could not perform post process computation: {name}".format(name = comp.name))
            raise e
            break

        result_element: et.Element = et.SubElement(global_results, comp.name)
        result_element.text = str(comp_result)

    tmp_path = compute_settings.tmp_path

    filename = tmp_path + "post_{r}.txt".format(r=str(random.randint(0, 2 ** 32)))
    while filename in os.listdir(tmp_path):
        filename = tmp_path + "post_{r}.txt".format(r=str(random.randint(0, 2 ** 32)))

    with open(filename, 'wb') as f:
        f.write(et.tostring(result))

    return filename


def initialise(prefix, _cell_dataframe):

    global mesh_prefix
    mesh_prefix = prefix

    global cell_dataframe
    cell_dataframe = _cell_dataframe


def get_color_dictionary(cell_df, cell_color_key, cell_colors, round_legend_labels=3):
    cell_color_key = cell_color_key


    if isinstance(cell_colors, Dict):
        categorical = True
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
    labels = []
    if len(color_dict) > 8:
        step = int(len(color_dict) / 8)
        items = list(color_dict.items())[0::step]
        legend_color_dict = dict(items)
    else:
        legend_color_dict = color_dict

    ms = plt.rcParams['legend.fontsize']
    ms = ms if not isinstance(ms,str) else 5

    for type, color in legend_color_dict.items():
        legend_items.append(
            plt.Line2D([0], [0], marker='o', color='w',markersize=ms,
                       markerfacecolor=color, markeredgewidth=0, linewidth=0)
        )
        labels.append(type)

    return color_dict, legend_items,labels, categorical


def get_rectangle_plane_mesh(u, res = (200,200)):

    coords = u.function_space().mesh().coordinates()
    coords_t = np.transpose(coords)

    x_limits = (
        np.min(coords_t[0]), np.max(coords_t[0])
    )
    y_limits = (
        np.min(coords_t[1]), np.max(coords_t[1])
    )
    rec_mesh = fcs.RectangleMesh(fcs.Point(x_limits[0], y_limits[0]), fcs.Point(x_limits[1], y_limits[1]), res[0], res[1])

    return rec_mesh, [x_limits,y_limits]

def make_paraview_images(conv_factor, compute_settings,  visual_conv = 1):

    conv_factor = conv_factor * visual_conv
    from thesis.main import __path__ as main_path

    try:
        for path in os.environ["PARAVIEW_PATH"].split(":"):
            if os.path.exists(path):
                paraview = path
                break
    except KeyError as e:
        warning("PARAVIEW_PATH env variable not set. cannot render images")
        return 0

    paraview_script = main_path._path[0] + "/paraview_render_script.py"
    pvbatch = os.path.join(paraview, "bin/pvbatch")

    import subprocess as sp
    import json
    if "type_name" in compute_settings.markers.keys():


        legend_title = compute_settings.legend_title
        cell_colors = compute_settings.cell_colors
        round_legend_labels = compute_settings.round_legend_labels

        color_dict, legend_items, labels, categorical = get_color_dictionary(cell_dataframe, "type_name",
                                                                             cell_colors,
                                                                             round_legend_labels=round_legend_labels)

        for type_name, entry in color_dict.items():
            if type_name in compute_settings.marker_lookup.keys():
                mesh_function_label = compute_settings.marker_lookup[type_name]
                if not str(mesh_function_label) in compute_settings.paraview_settings["lookup"].keys():
                    compute_settings.paraview_settings["lookup"][str(mesh_function_label)] = [type_name,color_dict[type_name],1]


        xdmf = os.path.join(compute_settings.path, compute_settings.solution_path)
        marker_path = os.path.join(compute_settings.path, compute_settings.markers["type_name"])
        settings_path = os.path.join(compute_settings.path, "paraview_settings_{f}.json".format(f = compute_settings.field))


        paraview_settings = compute_settings.paraview_settings

        # for k,v in compute_settings.paraview_settings.items():
        #         if isinstance(v,Dict):
        #             pass
        #         else:
        #             paraview_settings[k] = v

        with open(settings_path,"w") as f:
            paraview_settings["conversion_factor"] = conv_factor
            paraview_settings["field_name"] = "{f}({unit_name})".format(
                f = compute_settings.field, unit_name = compute_settings.unit_name
            )
            json.dump(paraview_settings, f, indent=1)

        img_path = "images/scan_{scan_index}/{field}/{time_index}/".format(
            field=compute_settings.field,
            scan_index=compute_settings.scan_index,
            time_index = compute_settings.time_index
        )
        img_path = os.path.join(compute_settings.path, img_path)

        os.makedirs(img_path,exist_ok=True)
        my_env = os.environ.copy()
        my_env["LD_LIBRARY_PATH"] = os.path.join(paraview,"lib")

        p = sp.Popen([pvbatch, paraview_script, xdmf, marker_path, img_path, settings_path], env = my_env)
        p.communicate()

def make_images(u, conv_factor, compute_settings, visual_conv = 1):

    conv_factor = conv_factor * visual_conv

    scan_index = compute_settings.scan_index
    time_index = compute_settings.time_index

    cell_color_key = compute_settings.cell_color_key
    legend_title = compute_settings.legend_title
    cell_colors = compute_settings.cell_colors
    round_legend_labels = compute_settings.round_legend_labels

    color_dict, legend_items, labels, categorical = get_color_dictionary(cell_dataframe, cell_color_key, cell_colors,
                                                                 round_legend_labels=round_legend_labels)

    cell_df = cell_dataframe.loc[
        (cell_dataframe["scan_index"] == scan_index) &
        (cell_dataframe["time_index"] == time_index)
        ]
    u.set_allow_extrapolation(True)
    rec_mesh, bbox = get_rectangle_plane_mesh(u)

    aspect = abs(bbox[0][0]-bbox[0][1]) / abs(bbox[1][0]-bbox[1][1])

    w = compute_settings.figure_width
    h = w / aspect

    from mpl_toolkits.axes_grid1 import Divider, Size

    pad_l = w * compute_settings.pad_l
    pad_r = w * compute_settings.pad_r

    fig = plt.figure(figsize=(w+pad_l+pad_r,h+pad_l + pad_r))

    vertical = [Size.Fixed(pad_l), Size.Fixed(h),Size.Fixed(pad_r)]
    height = [Size.Fixed(pad_l), Size.Fixed(w),Size.Fixed(pad_r)]

    divider = Divider(fig, (0, 0, 1, 1), height, vertical, aspect=False)
    pseudo_color_axes = fig.add_axes(divider.get_position(),axes_locator=divider.new_locator(
        nx=1,
        ny=1
    ))

    legend_fig = plt.figure(figsize=(8,8))

    legend_gs = GridSpec(1,2, legend_fig)
    legend_axes = [
        legend_fig.add_subplot(legend_gs[0]),
        # legend_fig.add_subplot(legend_gs[1]),
        # legend_fig.add_subplot(legend_gs[2])
    ]

    # for index, cell in cell_df.iterrows():
    #     p = [
    #         cell["x"],
    #         cell["y"]
    #     ]
    #
    #     key = cell[cell_color_key] if categorical else round(cell[cell_color_key], round_legend_labels)
    #     if not key in color_dict:
    #         color = "black"
    #         message("no color provided for key {k}".format(k=key))
    #     else:
    #         color = color_dict[key]
    #
    #     fig.gca().scatter(p[0],p[1],marker=5, color=color)
    #     # c = plt.Circle(p, 5, color=color)
    #     # fig.gca().add_artist(c)

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

    pseudo_color_axes.tricontourf(triang, slice_v, levels=100, norm = norm_cytokine)
    pseudo_color_axes.set_xlabel(r"x ($\mu m $)")
    pseudo_color_axes.set_ylabel(r"y ($\mu m $)")

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

        c = plt.Circle(p, 5, color=color, transform = fig.gca().transData)
        pseudo_color_axes.add_artist(c)



    if hasattr(compute_settings, "colorbar_range"):
        cb_cytokine = legend_fig.colorbar(mappable, ax = legend_axes[0],fraction = 0.5, aspect = 10)
        cb_cytokine.set_label("cytokine (nM)")
    else:
        cb_cytokine = legend_fig.colorbar(mappable)
        cb_cytokine.set_label("cytokine (nM)")

    if not categorical:
        vmin = min(color_dict.keys())
        vmax = max(color_dict.keys())

        norm = Normalize(vmin=vmin, vmax=vmax)
        cb = legend_fig.colorbar(ScalarMappable(norm=norm, cmap=cell_colors),ax = legend_axes[0], fraction = 0.5, aspect = 10)
        cb.set_label(cell_color_key)
    else:
        legend_axes[0].legend(handles=legend_items,labels=labels, title=legend_title )

    img_path = "images/scan_{scan_index}/{field}/".format(
        field=compute_settings.field,
        scan_index=scan_index
    )
    file_name = "{time_index}".format(
        time_index=time_index
    )
    os.makedirs(compute_settings.path + img_path, exist_ok=True)
    os.makedirs(os.path.join(compute_settings.path + img_path +"/legends/"),exist_ok=True)
    # fig.tight_layout()
    # legend_fig.tight_layout()
    legend_fig.savefig(compute_settings.path + img_path +"/legends/"+ file_name + ".pdf")
    fig.savefig(compute_settings.path + img_path + file_name + ".png", dpi=compute_settings.dpi)
