#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 12:56:59 2019

@author: Lukas Kiwitz
"""
import json
import multiprocessing as mp
import os
import random
from math import ceil
from typing import List, Dict

import KDEpy
import fenics as fcs
import lxml.etree as et
import mpi4py.MPI as MPI
import numpy as np
import pandas as pd

import MyError
import StateManager as st
from ParameterSet import ParameterSet
from myDictSorting import groupByKey
from my_debug import message


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
    def __init__(self)->None:


        self.file_path: str = ""
        """file path to volume file
        """
        self.field: str = ""
        self.dynamic: Dict = {}
        self.scan_index: int = 0
        self.time_index: float = 0
        self.tmp_path = ""


    def set_mesh(self,mesh: fcs.Mesh):
        self._mesh = mesh


    def set_u(self,u: fcs.Function):
        self._u: fcs.Function = u

    # noinspection PyPep8Naming,PyPep8Naming
    def set_V(self,V: fcs.FunctionSpace):
        self._V: fcs.FunctionSpace = V


    def set_boundary_markers(self,bm: fcs.MeshFunction):
        self._boundary_markers: fcs.MeshFunction = bm


    def get_mesh(self) ->fcs.Mesh:
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
        # self.pDicts = []
        self.cellDump = []
        self.out_tree_path = path + "postProcess.xml"
        self.path = path
        self.cell_dataframe: pd.DataFrame = pd.DataFrame()
        self.global_dataframe: pd.DataFrame = pd.DataFrame()
        self.cell_stats: pd.DataFrame = pd.DataFrame()
        self.unit_length_exponent: int = 1
        self.computations = [c(),grad(),SD()]

    def get_mesh_volume(self, mesh):

        sum = 0
        for cell in fcs.cells(mesh):
            sum += cell.volume()
        return sum
    # noinspection PyPep8Naming

    def compute(self, compute_settings: ComputeSettings) -> str:

        """
        performs computations according to compute_settings and return result as element tree string

        :param compute_settings:
        :return:
        """

         #converts concentration from nmol/(unit_length**3) to nmol/l

        result: et.Element = et.Element("File")
        global_results: et.Element = et.SubElement(result, "GlobalResults")


        result.set("field_name", str(compute_settings.field))
        result.set("dist_plot_path", str(compute_settings.file_path))


        result.append(compute_settings.parameters)
        result.set("scan_index", str(compute_settings.scan_index))
        result.set("time_index", str(compute_settings.time_index))
        result.set("time", str(compute_settings.time))

        mesh: fcs.Mesh = compute_settings.get_mesh()
        u: fcs.Function = compute_settings.get_u()
        mesh_volume = self.get_mesh_volume(mesh)
        mesh_volume_element = et.SubElement(global_results,"MeshVolume")
        mesh_volume_element.text = str(mesh_volume)

        p = ParameterSet("dummy",[])
        p.deserialize_from_xml(compute_settings.parameters)


        V_vec: fcs.VectorFunctionSpace = fcs.VectorFunctionSpace(mesh, "P", 1)
        V: fcs.FunctionSpace = fcs.FunctionSpace(mesh, "P", 1)
        grad: fcs.Function  = fcs.project(fcs.grad(u), V_vec, solver_type="gmres")

        for comp in compute_settings.computations:

            result_element: et.Element = et.SubElement(global_results, comp.name)
            result_element.text = str(comp(
                u,
                grad,
                get_concentration_conversion(self.unit_length_exponent),
                get_gradient_conversion(self.unit_length_exponent),
                mesh_volume,
                V=V,
                V_vec = V)
            )


        return et.tostring(result)

    # noinspection PyPep8Naming
    def job(self, compute_settings_list: List[ComputeSettings], output, thread_index: int, tmp_path: str):
        try:
            comm = MPI.COMM_WORLD
            local = comm.Dup()



            # boundary_markers = fcs.MeshFunction(
            #     "size_t", mesh, mesh.topology().dim() - 1)
            # with fcs.HDF5File(local, sub_domain_cache, "r") as f:
            #     f.read(boundary_markers, "/boundaries")
            # boundary_markers = boundary_markers

            result_list = []
            for n, compute_settings in enumerate(compute_settings_list):

                compute_settings.computations = self.computations

                mesh = fcs.Mesh()
                with fcs.XDMFFile(compute_settings.path +"/"+ compute_settings.mesh_path) as f:
                    f.read(mesh)
                mesh = mesh

                V = fcs.FunctionSpace(mesh, "P", 1)

                u: fcs.Function  = fcs.Function(V)
                with fcs.HDF5File(local, compute_settings.path +"/"+compute_settings.file_path, "r") as f:
                    f.read(u, "/" + compute_settings.field)
                if not compute_settings:
                    continue
                compute_settings.set_mesh(mesh)
                compute_settings.set_V(V)
                compute_settings.set_u(u)

                message(
                    "Process {thread}: Reading file {file} ({n}/{tot})".format(thread=thread_index, file=compute_settings.file_path, n=n, tot=len(compute_settings_list))
                )

                data_out: str = self.compute(compute_settings)
                result_list.append(str(data_out))
            message("Process {index} has finished computation".format(index=thread_index))
            filename = tmp_path + "post_{r}.txt".format(r=str(random.randint(0, 2 ** 32)))
            while filename in os.listdir(tmp_path):
                filename  = tmp_path+"post_{r}.txt".format(r=str(random.randint(0,2**32)))
            message("Process {index} writing results to file {f}".format(index=thread_index,f=filename))
            f = open(filename,'x')
            f.write(json.dumps(result_list))
            f.close()
            output.put(filename)

        except Exception as e:
            message(str(e))
            output.put(e)

    def write_post_process_xml(self, threads,debug=False):
        """
        runs compute-function for all scans an writes result to xml file

        """

        assert type(threads) == int

        # initializes state manager from scan log
        self.stateManager = st.StateManager(self.path)
        self.stateManager.load_xml()

        tmp_path = self.path+"tmp/"
        os.makedirs(tmp_path,exist_ok=True)

        scatter_list: List[ComputeSettings] = []


        # loads timestep logs
        for scan_index,scan_sample in enumerate(self.stateManager.element_tree.findall("ScanContainer/ScanSample")):
            parameters = scan_sample.find("Parameters")


            for field_step in scan_sample.findall("TimeSeries/Step/Fields/Field"):

                compute_settings: ComputeSettings = ComputeSettings()
                compute_settings.path = self.path
                compute_settings.file_path = field_step.get("dist_plot_path")
                compute_settings.field = field_step.get("field_name")
                compute_settings.mesh_path = field_step.get("mesh_path")

                compute_settings.parameters = parameters
                compute_settings.scan_index = scan_index
                compute_settings.time_index = field_step.getparent().getparent().get("time_index")
                compute_settings.time = field_step.getparent().getparent().get("time")

                scatter_list.append(compute_settings)

        message("distributing to {threads} threads".format(threads=threads))
        if debug:
            scatter_list = scatter_list[0:debug]
        size = ceil(len(scatter_list) / threads)
        partitioned_list = [scatter_list[x:x + size]
                           for x in range(0, len(scatter_list), size)]
        output = mp.Queue(threads)
        jobs = [mp.Process(target=self.job, args=(i, output, index, tmp_path)) for index, i in enumerate(partitioned_list)]
        for j in jobs:
            j.start()
        from time import time
        start = time()
        timeout = 24*60*60

        while True:
            if (time() - start) < timeout:
                running = False
                for j in jobs:
                    if j.is_alive():
                        running = True
                        break
                if not running:
                    message("collecting distributed tasks")
                    break
            else:
                raise MyError.SubProcessTimeout(timeout)
                break

        # file_list: List[str] = [output.get(True, 10) for j in jobs]
        file_list: List[str] = []
        from queue import Empty
        for job in jobs:
            try:
                time_out = 10 * 60
                result = output.get(True, time_out)
                file_list.append(result)
            except Empty as e:
                message("Could not retrieve result from queue. Machine busy?".format(t=time_out))

        for file in file_list:
            if not type(file) == str:
                message("A Worker fininshed with Error: {e}".format(e=file))
                file_list.remove(file)

        message("Collected results from {l} Processes".format(l=1 + len(file_list)))
        result_list = []
        for file in file_list:
            f = open(file, "r")
            result_list.append(json.load(f))
        message("successfully collected distributed task")
        flattend_list: List[str] = []

        for i in result_list:
            for o in i:
                flattend_list.append(et.XML(o[2:-1]))
        indexed_list = [
            {"scan_index": i.get("scan_index"),
             "time_index": i.get("time_index"),
             "time": i.get("time"),
             "entry": i}
            for i in flattend_list]
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
                    time.set("time",i["time"])
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
                        "surf_c": self.cell_dataframe.loc[filter(self.cell_dataframe)]["{field}_surf_c".format(field=field_name)].mean()
                    }
                    parameter_set = ParameterSet("dummy_set",[])
                    parameter_set.deserialize_from_xml(file.find("Parameters/ParameterSet[@name='dynamic']"))
                    d.update(parameter_set.get_as_dictionary())


                    for v in g_values:
                        d.update({v.tag:float(v.text)})
                    result.append(d)

        return pd.DataFrame(result)

    def get_kde_estimators(self, n ,ts, time_index, type_names):

        kernels = {}
        message("computing kde for time series: {n} and timestep {t}".format(n=n, t=time_index))
        for type_name in type_names:
            inital_cells = ts.loc[(ts["time_index"] == time_index) & (ts["type_name"] == type_name)]
            if inital_cells.shape[0] == 0:
                break
            elif inital_cells.shape[0] == 1:
                data = np.array([inital_cells["x"].iloc[0], inital_cells["y"].iloc[0], inital_cells["z"].iloc[0]])
                kernel = KDEpy.NaiveKDE("tri", bw=10e-2).fit(data)
            else:
                data = np.array([inital_cells["x"], inital_cells["y"], inital_cells["z"]]).T
                kernel = KDEpy.TreeKDE("tri", bw=10e-2).fit(data)

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
                for time_index,step in ts.groupby(["time_index"]):

                    kde_time_index = time_index-1 if time_index > 0 else 0
                    kernels.append(self.get_kde_estimators(scan_index,ts, kde_time_index, result["type_name"].unique()))

                for time_index, step in ts.groupby(["time_index"]):

                    for type_name, kernel in kernels[time_index].items():
                        positions = np.array([step["x"], step["y"], step["z"]]).T

                        scores = kernel.evaluate(positions).T
                        scores =  pd.Series(scores)
                        scores.index = step.index

                        step.insert(step.shape[1],"{type_name}_score".format(type_name=type_name),scores)

                    for type_name, kernel in kernels[0].items():
                        positions = np.array([step["x"], step["y"], step["z"]]).T

                        scores = kernel.evaluate(positions).T
                        scores =  pd.Series(scores)
                        scores.index = step.index

                        step.insert(step.shape[1],"{type_name}_score_init".format(type_name=type_name),scores)

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

            no_init = group#group.loc[~group["id"].isin(ids)]

            for old in pd.Series(group.columns).str.extract("(.*_score)").dropna()[0].unique():
                new = "{old}_norm".format(old=old)
                group.insert(group.shape[1],new, group[old] / float(no_init.mean()[old]))

            for old in pd.Series(group.columns).str.extract("(.*_score_init)").dropna()[0].unique():
                new = "{old}_norm".format(old=old)
                group.insert(group.shape[1],new, group[old] / float(no_init.mean()[old]))

            result = result.append(group)

        return result

    def get_stats(self):

        # temp_cell_df = self.cell_dataframe.drop(columns = ["x","y","z","id"])
        grouped = self.cell_dataframe.groupby(["type_name", "time_index", "scan_index"], as_index=False)

        des = grouped.describe(include = 'all')

        # des = des.reset_index()

        return des

    def make_dataframes(self,kde=False):
        self.cell_dataframe = self.get_cell_dataframe(kde=kde)
        self.global_dataframe = self.get_global_dataframe()
        # self.cell_stats = self.get_stats()

def get_concentration_conversion(unit_length_exponent: int):

    assert isinstance(unit_length_exponent,int)
    return 10**(-1 * (unit_length_exponent * 3 + 3))

def get_gradient_conversion(unit_length_exponent: int):

    assert isinstance(unit_length_exponent, int)
    exp = (-1 * (unit_length_exponent * 3 + 3)) -  1
    exp -= (6 + unit_length_exponent)# to nM/um

    return 10**(exp)

class SD:

    def __init__(self):
        self.name = "SD"

    def __call__(self, u, grad, c_conv, grad_conv, mesh_volume, **kwargs):

        sd: float = np.std(np.array(u.vector())) * c_conv

        return sd

class c:

    def __init__(self):
        self.name = "Concentration"

    def __call__(self, u, grad, c_conv, grad_conv, mesh_volume, **kwargs):


        concentration: float = (fcs.assemble(u * fcs.dX)/ mesh_volume)*c_conv
        return  concentration

class grad:

    def __init__(self):
        self.name = "Gradient"

    def __call__(self, u, grad, c_conv, grad_conv, mesh_volume, **kwargs):

        gradient: float = fcs.assemble(fcs.sqrt(fcs.dot(grad, grad)) * fcs.dX) * grad_conv / mesh_volume
        return gradient