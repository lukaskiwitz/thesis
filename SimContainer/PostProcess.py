#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 12:56:59 2019

@author: kiwitz
"""
import StateManager as st

import fenics as fcs
import numpy as np
import json
import lxml.etree as et
import mpi4py.MPI as MPI
import multiprocessing as mp
from copy import deepcopy, copy
from math import ceil
import pandas as pd
from typing import List, Tuple, Dict

from myDictSorting import groupByKey

class ComputeSettings:
    """
    Class to store computation input and output

    Attributes
    :param file_path: file path to volume file
    :type file_path: str

    :param field: field name in h5 file
    :type field: str

    :param cell_data: cell data
    :type cell_data List[Dict]

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
        self.field: str = ""
        self.cell_data: List[Dict] = []
        self.dynamic: Dict = {}
        self.scan_index: int = 0
        self.time_index: float = 0
        self.compute_gradient: bool = False
        self.compute_concentration: bool = False
        self.compute_surface_concentration: bool = False

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

        if compute_settings.compute_gradient:
            V_vec: fcs.VectorFunctionSpace = fcs.VectorFunctionSpace(mesh, "P", 1)
            grad: fcs.Function  = fcs.project(fcs.grad(u), V_vec, solver_type="gmres")
            gradient: float = fcs.assemble(fcs.sqrt(fcs.dot(grad, grad)) * fcs.dX) * 10 ** 8
            gradient_result: et.Element = et.SubElement(global_results, "gradient")
            gradient_result.text = str(gradient)
        if compute_settings.compute_concentration:
            concentration: float= fcs.assemble(u * fcs.dX) * 10 ** 9
            concentration_result: et.Element = et.SubElement(global_results, "concentration")
            concentration_result.text = str(concentration)

        if compute_settings.compute_surface_concentration:
            for i, cell in enumerate(compute_settings.cell_data):
                cell_element: et.Element = et.SubElement(cell_results, "cell")
                patch: et.Element = et.SubElement(cell_element, "patch")
                center: et.Element = et.SubElement(cell_element, "center")
                patch.text = str(cell["patch"])
                center.text = json.dumps(list(cell["center"]))

                ds: fcs.Measure = fcs.Measure("ds", domain=mesh,
                                 subdomain_data=boundary_markers)
                v: float = (fcs.assemble(u * ds(cell["patch"])) / (4 * np.pi * 0.05 ** 2) * 10 ** 9)

                surface_concentration: et.Element = et.SubElement(cell_element, "surface_concentration")
                surface_concentration.text = str(v)

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

                print(
                    "reading file {file} ({n}/{tot})".format(file=compute_settings.file_path, n=n, tot=len(compute_settings_list))
                )

                data_out: str = self.compute(compute_settings)
                result_list.append(deepcopy(data_out))
            print("thread no {index}".format(index=thread_index))
            output.put(result_list)
        except Exception as e:
            print(e)

    def dump(self, path, threads):
        # initializes state manager from scan log
        self.stateManager = st.StateManager(path)
        self.stateManager.loadXML()
        self.ext_cache = self.stateManager.elementTree.find(
            "/cellDump/field/extCache").text
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
                compute_settings.compute_surface_concentration = True
                compute_settings.compute_concentration = True
                compute_settings.compute_gradient = True

                scatter_list.append(compute_settings)

        print("scatter")
        # scatter_list = scatter_list[0:4]
        size = ceil(len(scatter_list) / threads)
        partitioned_list = [scatter_list[x:x + size]
                           for x in range(0, len(scatter_list), size)]
        output = mp.Queue(threads)
        jobs = [mp.Process(target=self.job, args=(i, self.ext_cache, self.subdomaincache, output, index)) for index, i in enumerate(partitioned_list)]
        for j in jobs:
            j.start()

        for j in jobs:
            try:
                j.join(60)
            except Exception:
                print("Join Timeout")
        print("joined jobs")
        result_list: List[str] = [output.get(True, 60) for j in jobs]
        print("collected output from threads")
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
            #            print(s[0])
            scan.set("i", str(s[0][0]["scanIndex"]))
            for t in s:
                for i in t:
                    time = et.SubElement(scan, "timeStep")
                    time.set("i", i["timeIndex"])
                    time.append(i["entry"])

        tree.write(self.out_tree_path, pretty_print=True)

    def prep_data(self) -> pd.DataFrame:
        in_tree = et.parse(self.out_tree_path)
        frames = []
        for scan in in_tree.findall("/scan"):
            scan_index = float(scan.get("i"))
            time_steps = np.unique([int(i.get("i"))
                                   for i in scan.findall("./timeStep")])

            for t in time_steps:
                files = scan.findall("./timeStep[@i='{t}']/file".format(t=t))
                # Todo get number of fields and cells dynamicly
                offset = 3  # one file per field
                cell_results = np.empty((1500, len(files) + offset))
                # name_list = []
                for cellIndex, cell in enumerate(
                        files[0].findall("./cellResults/cell")
                ):

                    x = json.loads(cell.find("./center").text)[0]
                    cell_results[cellIndex, 0] = x
                    cell_results[cellIndex, 1] = t
                    cell_results[cellIndex, 2] = scan_index
                    name_list = []
                    for o, file in enumerate(files):
                        cell = file.findall("./cellResults/cell")[cellIndex]
                        field_name = file.get("field")
                        name_list.append(field_name)
                        cell_results[cellIndex, o + offset] = float(
                            cell.find("./surfaceConcentration").text
                        )

                cell_frame = pd.DataFrame(cell_results, columns=[
                                                                  "x", "time", "scanIndex"] + name_list)
                frames.append(cell_frame)
        # join dataframes
        result = frames[0]
        for i in range(len(frames) - 1):
            result = result.append(frames[i + 1])
        return result

    def prep_global_data(self) -> pd.DataFrame:
        in_tree = et.parse(self.out_tree_path)
        for scan in in_tree.findall("/scan"):
            time_steps = np.unique([int(i.get("i"))
                                    for i in scan.findall("./timeStep")
                                    ])
        result = []
        for step in in_tree.findall("/scan"):
            files = step.findall("./timeStep/file")
            for file in files:
                g_values = file.find("./global")
                d = {
                    "scanIndex": int(file.get("scanIndex")),
                    "timeIndex": int(file.get("timeIndex")),
                    "fieldName": file.get("field")
                }
                for p in json.loads(file.get("dynamic")):
                    d[p["name"]] = p["value"]
                for v in g_values.getchildren():
                    d[v.tag] = float(v.text)
                result.append(d)
        return pd.DataFrame(result)


PATH_LIST = [
    # {"path":"/extra/kiwitz/results_parameter_scan_Diffusion/","key":"D"},
    # {"path":"/extra/kiwitz/results_parameter_scan_fraction/","key":"fraction"},
    # {"path":"/extra/kiwitz/results_parameter_scan_kd/","key":"decay"},
    # {"path":"/extra/kiwitz/results_parameter_scan_kON/","key":"k_{on}"},
    # {"path":"/extra/kiwitz/results_parameter_scan_q_s/","key":"q Secretors"},
    # {"path":"/extra/kiwitz/results_parameter_scan_R_il2_f/","key":"R il2 f"},
    {"path":"/extra/kiwitz/results_parameter_scan_R_il2_s/","key":"R il2 s"}
             ]
for i in PATH_LIST:
    path = i["path"]
    pp = PostProcessor(path)
    pp.prep_global_data().to_hdf(path + 'global_dataframe.h5', key="data", mode="w")
    pp.prep_data().to_hdf(path + 'dataframe.h5', key="data", mode="w")
# pp.dump(PATH, 64)
