#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 12:56:59 2019

@author: Lukas Kiwitz
"""
import json
import multiprocessing as mp
from copy import deepcopy
from math import ceil
from typing import List, Dict

import KDEpy
import fenics as fcs
import lxml.etree as et
import mpi4py.MPI as MPI
import numpy as np
import pandas as pd
from scipy.constants import N_A

import StateManager as st
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
        self.cell_data: List[Dict] = []
        self.dynamic: Dict = {}
        self.scan_index: int = 0
        self.time_index: float = 0

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
        self.pDicts = []
        self.cellDump = []
        self.out_tree_path = path + "postProcess.xml"
        self.path = path
        self.cell_dataframe: pd.DataFrame = pd.DataFrame()
        self.global_dataframe: pd.DataFrame = pd.DataFrame()
        self.cell_stats: pd.DataFrame = pd.DataFrame()

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
        result: et.Element = et.Element("file")
        global_results: et.Element = et.SubElement(result, "global")
        cell_results: et.Element = et.SubElement(result, "cell_results")

        result.set("field", str(compute_settings.field))
        result.set("path", str(compute_settings.file_path))
        result.set("dynamic", json.dumps(compute_settings.dynamic))
        result.set("scanIndex", str(compute_settings.scan_index))
        result.set("timeIndex", str(compute_settings.time_index))

        mesh: fcs.Mesh = compute_settings.get_mesh()
        boundary_markers: fcs.MeshFunction = compute_settings.get_boundary_markers()
        u: fcs.Function = compute_settings.get_u()
        mesh_volume = self.get_mesh_volume(mesh)
        mesh_volume_element = et.SubElement(global_results,"mesh_volume")
        mesh_volume_element.text = str(mesh_volume)


        V_vec: fcs.VectorFunctionSpace = fcs.VectorFunctionSpace(mesh, "P", 1)
        grad: fcs.Function  = fcs.project(fcs.grad(u), V_vec, solver_type="gmres")
        gradient: float = fcs.assemble(fcs.sqrt(fcs.dot(grad, grad)) * fcs.dX) * 1e8
        gradient = gradient/mesh_volume
        gradient_result: et.Element = et.SubElement(global_results, "gradient")
        gradient_result.text = str(gradient)

        concentration: float= fcs.assemble(u * fcs.dX) * 1e9
        concentration = concentration/mesh_volume
        concentration_result: et.Element = et.SubElement(global_results, "concentration")
        concentration_result.text = str(concentration)

        sd: float = np.std(np.array(u.vector()))
        sd_result: et.Element = et.SubElement(global_results,"sd")
        sd_result.text = str(sd)

        return et.tostring(result)

    # noinspection PyPep8Naming
    def job(self, compute_settings_list: List[ComputeSettings], ext_cache: str, sub_domain_cache: str, output,thread_index: int):
        try:
            comm = MPI.COMM_WORLD
            local = comm.Dup()

            mesh = fcs.Mesh()
            with fcs.XDMFFile(ext_cache + ".xdmf") as f:
                f.read(mesh)
            mesh = mesh

            # V = fcs.FunctionSpace(mesh, "P", 1)
            V = fcs.FunctionSpace(mesh, "P", 1)

            boundary_markers = fcs.MeshFunction(
                "size_t", mesh, mesh.topology().dim() - 1)
            with fcs.HDF5File(local, sub_domain_cache, "r") as f:
                f.read(boundary_markers, "/boundaries")
            boundary_markers = boundary_markers

            result_list = []
            for n, compute_settings in enumerate(compute_settings_list):

                u: fcs.Function  = fcs.Function(V)
                with fcs.HDF5File(local, compute_settings.file_path, "r") as f:
                    f.read(u, "/" + compute_settings.field)
                if not compute_settings:
                    continue
                compute_settings.set_mesh(mesh)
                compute_settings.set_V(V)
                compute_settings.set_boundary_markers(boundary_markers)
                compute_settings.set_u(u)

                message(
                    "reading file {file} ({n}/{tot})".format(file=compute_settings.file_path, n=n, tot=len(compute_settings_list))
                )

                data_out: str = self.compute(compute_settings)
                result_list.append(deepcopy(data_out))
            message("Thread {index} has finished computation".format(index=thread_index))
            output.put(result_list)
        except Exception as e:
            message(e)
            output.put([])

    def write_post_process_xml(self, threads,debug=False):
        """
        runs compute-function for all scans an writes result to xml file

        """

        # initializes state manager from scan log
        self.stateManager = st.StateManager(self.path)
        self.stateManager.loadXML()
        self.ext_cache = self.stateManager.elementTree.find(
            "/cellDump/field/ext_cache").text
        self.subdomaincache = self.stateManager.elementTree.find(
            "/cellDump/field/subdomains").text

        for s in self.stateManager.elementTree.findall("/scans/scan"):
            self.pDicts.append(self.stateManager.getParametersFromElement(s))

        for s in self.stateManager.elementTree.findall("/cellDump/field/cell"):
            patch = int(s.find("patch").text)
            center = json.loads(s.find("center").text)
            self.cellDump.append({"patch": patch, "center": center})

        scatter_list: List[ComputeSettings] = []

        #        for field in fields:
        cell_data = self.stateManager.elementTree.findall(
            "/cellDump/field[@name='il2']/cell")
        cell_data = [{"patch": int(i.find("patch").text), "center": json.loads(
            i.find("center").text)} for i in cell_data]
        # loads timestep logs
        for scan in self.stateManager.elementTree.findall("scans/scan"):
            dynamic = scan.findall("parameters/dynamic/parameter")
            dynamic = [{"name": i.get("name"), "value": i.text}
                       for i in dynamic]

            for step in scan.findall("timeSeries/field/step"):

                compute_settings: ComputeSettings = ComputeSettings()
                compute_settings.file_path = step.find("distPlotPath").text
                compute_settings.field = step.getparent().get("name")
                compute_settings.cell_data = cell_data
                compute_settings.dynamic = dynamic
                compute_settings.scan_index = scan.get("i")
                compute_settings.time_index = step.get("t")

                scatter_list.append(compute_settings)

        message("distributing to {threads} threads".format(threads=threads))
        if debug:
            scatter_list = scatter_list[0:debug]
        size = ceil(len(scatter_list) / threads)
        partitioned_list = [scatter_list[x:x + size]
                           for x in range(0, len(scatter_list), size)]
        output = mp.Queue(threads)
        jobs = [mp.Process(target=self.job, args=(i, self.ext_cache, self.subdomaincache, output, index)) for index, i in enumerate(partitioned_list)]
        for j in jobs:
            j.start()

        for j in jobs:
            try:
                j.join(600)
            except Exception:
                message("Join Timeout")
        message("collecting distributed tasks")
        result_list: List[str] = [output.get(True, 60) for j in jobs]
        message("successfully collected distributed task")
        flattend_list: List[str] = []

        for i in result_list:
            for o in i:
                flattend_list.append(et.XML(o))
        indexed_list = [
            {"scanIndex": i.get("scanIndex"),
             "timeIndex": i.get("timeIndex"),
             "entry": i}
            for i in flattend_list]
        indexed_list = groupByKey(indexed_list, ["scanIndex"])
        for i, e in enumerate(indexed_list):
            indexed_list[i] = groupByKey(e, ["timeIndex"])

        post_process_result = et.Element("postProcess")
        tree = et.ElementTree(element=post_process_result)
        for s in indexed_list:
            scan = et.SubElement(post_process_result, "scan")
            #            message(s[0])
            scan.set("i", str(s[0][0]["scanIndex"]))
            for t in s:
                for i in t:
                    time = et.SubElement(scan, "timeStep")
                    time.set("i", i["timeIndex"])
                    time.append(i["entry"])
        message("writing post process output to {p}".format(p=self.out_tree_path))
        tree.write(self.out_tree_path, pretty_print=True)

    def get_global_dataframe(self) -> pd.DataFrame:
        message("collecting global dataframe from post_process.xml")
        in_tree = et.parse(self.out_tree_path)
        for scan in in_tree.findall("/scan"):
            time_steps = np.unique([int(i.get("i"))
                                    for i in scan.findall("./timeStep")
                                    ])
        result = []
        assert not self.cell_dataframe.empty

        for step in in_tree.findall("/scan"):
            files = step.findall("./timeStep/file")
            for file in files:
                g_values = file.find("./global")
                scan_index = int(file.get("scanIndex"))
                timeIndex = int(file.get("timeIndex"))
                field_name = file.get("field")
                filter = lambda x: (x["t"] == timeIndex) & \
                                   (x["scan_index"] == scan_index)

                d = {
                    "scanIndex": scan_index,
                    "timeIndex": timeIndex,
                    "field_name": field_name,
                    "surf_c": self.cell_dataframe.loc[filter(self.cell_dataframe)]["surf_c_{field}".format(field=field_name)].mean()
                }
                for p in json.loads(file.get("dynamic")):
                    d[p["name"]] = float(p["value"])
                for v in g_values.getchildren():
                    d[v.tag] = float(v.text)
                result.append(d)

        return pd.DataFrame(result)

    def get_cell_dataframe(self, kde=False):

        self.stateManager = st.StateManager(self.path)
        self.stateManager.loadXML()

        result: pd.DataFrame = self.stateManager.get_cell_ts_data_frame()

        result["x"] = result["center"].apply(lambda x: x[0])
        result["y"] = result["center"].apply(lambda x: x[1])
        result["z"] = result["center"].apply(lambda x: x[2])
        result = result.drop(columns=["center"])

        result = result.groupby("id").apply(lambda x: x.ffill().bfill()).drop_duplicates()

        for name in self.stateManager.get_field_names():
            result["R_{f}".format(f=name)] = result["R_il2"].div(N_A ** -1 * 1e9)  # to mol/cell
            result["surf_c_{f}".format(f=name)] = result["surf_c_{f}".format(f=name)].mul(1e9)  # to nM
            result["surf_g_{f}".format(f=name)] = result["surf_g_{f}".format(f=name)].mul(1e8)  # to nM


        result["t"] = result["t"].apply(lambda x: float(x))
        result["scan_index"] = result["scan_index"].apply(lambda x: float(x))

        """---------------------------"""

        if kde:
            message("running kernel density estimation")
            r_grouped = result.groupby(["scan_index"], as_index=False)
            kde_result = pd.DataFrame()
            for n, ts in r_grouped:
                kernels = {}
                message("computing kde for time series: {n}".format(n=n))
                for type_name in result["type_name"].unique():
                    inital_cells = ts.loc[(ts["t"] == 1) & (ts["type_name"] == type_name)]
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

                for type_name, kernel in kernels.items():
                    positions = np.array([ts["x"], ts["y"], ts["z"]]).T
                    scores = kernel.evaluate(positions).T
                    ts["{type_name}_score".format(type_name=type_name)] = pd.Series(scores)
                kde_result = kde_result.append(ts)
                result = self._normalize_cell_score(kde_result)

        return result

    def _normalize_cell_score(self, x):
        ids = x.loc[
            (x["type_name"] != "Tn") &
            (x["t"] == 1)
            ]["id"].unique()
        x = x.groupby(["t", "scan_index"], as_index=False)

        result = pd.DataFrame()
        for group in x:
            i = group[0]
            group = group[1]

            no_init = group.loc[~group["id"].isin(ids)]
            # count = group.groupby(["type_name"],as_index=False).count()
            # count = count.drop(count.columns.drop(["type_name","n"]),axis=1)

            group["Tn_score_norm"] = group["Tn_score"] / float(no_init.mean()["Tn_score"])
            group["Tfh_score_norm"] = group["Tfh_score"] / float(no_init.mean()["Tfh_score"])
            group["Th1_score_norm"] = group["Th1_score"] / float(no_init.mean()["Th1_score"])
            result = result.append(group)

        return result

    def get_stats(self):

        temp_cell_df = self.cell_dataframe.drop(columns = ["x","y","z","id"])
        grouped = temp_cell_df.groupby(["type_name", "t", "scan_index"], as_index=True)


        des = grouped.describe(include = 'all')
        # des = des.reset_index()

        return des

    def make_dataframes(self,kde=False):
        self.cell_dataframe = self.get_cell_dataframe(kde=kde)
        self.global_dataframe = self.get_global_dataframe()
        self.cell_stats = self.get_stats()



