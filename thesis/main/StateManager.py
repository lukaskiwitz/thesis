#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 14:34:19 2019

@author: Lukas Kiwitz
"""
import json
import os
import time
import traceback
from copy import deepcopy
from functools import reduce
from typing import *

import lxml.etree
import lxml.etree as ET
import mpi4py.MPI as MPI
import numpy as np
import pandas as pd
from lxml import etree
from scipy.stats.mstats import gmean
from tqdm import tqdm

from thesis.main.Entity import Cell
from thesis.main.MyScenario import MyScenario
from thesis.main.ParameterSet import ParameterSet, GlobalCollections, GlobalParameters
from thesis.main.ScanContainer import ScanContainer, ScanSample
from thesis.main.SimContainer import SimContainer
from thesis.main.TaskRecord import ClassRecord
from thesis.main.my_debug import message, warning, critical

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


def outputParse(v):
    if type(v) == np.ndarray:
        return json.dumps(v.tolist())
    else:
        return json.dumps(v)


ElementTree = lxml.etree.ElementTree
Element = lxml.etree.Element


class StateManager:
    """
    This class manages the state of the simulation. Its main use is to run parameter scans in an organised manner and
    produces the xml files necessary for the use of PostProcessor.


    :ivar path:
    :ivar ruse:
    :ivar scan_folder_pattern:
    :ivar element_tree:
    :ivar scan_container:
    :ivar sim_container:
    :ivar T:
    :ivar N:
    :ivar compress_log_file:
    :ivar globarl_collections:
    :ivar global_parameters:
    :ivar record:
    :ivar progress_bar:
    :ivar eta_estimates:

    :vartype path: str
    :vartype ruse: pd.DataFrame
    :vartype scan_folder_pattern: str
    :vartype element_tree: ET.ElementTree
    :vartype scan_container: ScanContainer
    :vartype sim_container: SimConatiner
    :vartype T: List[float]
    :vartype dt: float
    :vartype N: int
    :vartype compress_log_file: bool
    :vartype global_collections: GlobalCollections
    :vartype global_parameters: GlobalParameters
    :vartype record: ClassRecord
    :vartype progress_bar:
    :vartype eta_estimates: List[float]
    """

    def __init__(self, path: str):

        """
        :param path: Path to simulation directoy
        :return None
        """
        self.path: str = path
        self.ruse = []

        self.scan_tree: MyScanTree = MyScanTree(path)
        self.scan_tree.compress_log_file: bool = True
        self.scan_folder_pattern: str = "scan_{n}/"

        self.scan_container: ScanContainer = None
        self.sim_container: SimContainer = None
        self.scenario: MyScenario = None
        self.T: List[float] = None
        self.dt: float = 1
        self.N: int = 10
        self.record: ClassRecord = ClassRecord("StateManager")
        self.progress_bar = None
        self.eta_estimates: List[float] = []
        self.marker_lookup: Mapping[str, int] = {}
        self.debug: bool = False

    def get_scan_folder(self, n: int) -> str:
        return self.path + self.scan_folder_pattern.format(n=n)

    def get_cell_ts_data_frame(self, time_indices: List[int] = None,scan_indicies = None, n_processes: int = 1 ) -> pd.DataFrame:

        def init():
            global path
            path = self.path

        result = []

        scans = self.scan_tree.get_scan_elements(scan_indicies)

        import multiprocessing as mp

        scatter_list = []
        for n, scan in enumerate(scans):
            scan = self.scan_tree.rebuild_timesteps(scan)

            scatter_list.append((n, ET.tostring(scan), time_indices, self.path))

        n_processes = n_processes if n_processes <= os.cpu_count() else os.cpu_count()

        message("State Manager: Distributing {total} scans to {n_processes} processes".format(n_processes=n_processes,
                                                                                           total=len(scans)))


        with mp.Pool(processes=n_processes, initializer=init) as p:
            result = p.map(target, scatter_list)
            while [] in result:
                result.remove([])

        # result = []
        # for sc in scatter_list:
        #     result.append(target(sc))
        result = reduce(lambda x, v: x + v, result, [])

        return pd.DataFrame(result)

    def get_scan_sample(self, i: int) -> ScanSample:

        scan_container = self.scan_tree.deserialize_from_element_tree()
        assert i <= len(scan_container.scan_samples) - 1
        sample = scan_container.scan_samples[i]
        return sample

    def apply_sample_flags(self, sc: SimContainer, scan_index: int):

        sample = self.get_scan_sample(scan_index)
        for f in sc.fields:
            f.remesh_scan_sample = sample.remesh_scan_sample
            f.remesh_timestep = sample.remesh_timestep

    def update_sim_container(self, sc: SimContainer, scan_index: int, model_index: int) -> ParameterSet:

        sc.path = os.path.join(self.get_scan_folder(scan_index), "m_{}/".format(model_index))
        sc.relative_path = os.path.relpath(sc.path, os.path.commonpath([self.get_scan_folder(scan_index), self.path]))
        sample = self.get_scan_sample(scan_index)

        assert hasattr(sc, "default_sample")

        sc.p.update(sc.default_sample.p)
        sc.p.update(sample.p)

        for f in sc.fields:
            f.apply_sample(sc.default_sample)
            f.apply_sample(sample)

        for e in sc.entity_list:
            e.p.update(sc.default_sample.p, override=True)
            e.p.update(sample.p, override=True)

        for entity_type in sc.default_sample.entity_types:
            sc.add_entity_type(entity_type)

        for entity_type in sample.entity_types:
            sc.add_entity_type(entity_type)

        for e in sc.entity_list:
            for cell_type in sc.default_sample.entity_types:
                if e.type_name == cell_type.name:
                    e.change_type = cell_type.name

        for e in sc.entity_list:
            for cell_type in sample.entity_types:
                if e.type_name == cell_type.name:
                    e.change_type = cell_type.name

        return deepcopy(sample.p)

    def estimate_time_remaining(self, sc: SimContainer, scan_index: int, time_index: int, S: List, T: List) -> None:

        sample = self.record.child_tasks["run"].child_tasks["scan_sample"]
        time_step = sample.child_tasks["SimContainer"].child_tasks["run"].child_tasks["step"]
        step_durations = [i["end"] - i["start"] for i in time_step.history]
        mean_step_duration = gmean(step_durations)

        n_steps = (len(T) - time_index)
        n_scans = len(S) - scan_index

        time_series_eta = time.time() + n_steps * mean_step_duration
        scan_eta = time.time() + (n_scans * len(T) * mean_step_duration + n_steps * mean_step_duration)
        self.eta_estimates.append([time.time(), scan_eta])

        np.save(os.path.join(self.get_records_path(), "eta"), np.array(self.eta_estimates))

        scan_eta = time.strftime("%H:%M:%S %m/%d/%y", time.localtime(scan_eta))
        time_series_eta = time.strftime("%H:%M:%S", time.localtime(time_series_eta))

        self.time_series_bar.update(1)
        self.scan_bar.postfix = "ETA: {eta}".format(eta=scan_eta)
        self.time_series_bar.postfix = "ETA: {eta}".format(eta=time_series_eta)
        # message("ETA: {eta}".format(eta = eta))

    def run(self, ext_cache: str = "", model_index=0) -> None:

        run_task = self.record.start_child("run")
        sample_task = run_task.start_child("scan_sample")

        self.scan_tree.serialize_to_element_tree(self.scan_container)
        self.scan_container = self.scan_tree.deserialize_from_element_tree()

        n_samples = len(self.scan_container.scan_samples)
        self.scan_bar = tqdm(desc="Parameter Scan", total=n_samples, dynamic_ncols=True)
        self.time_series_bar = tqdm(desc="Time Series   ", initial=0, dynamic_ncols=True)

        self.scan_tree.load_xml()
        self.scan_tree.write_element_tree()

        for scan_index in range(n_samples):
            if not sample_task.running:
                sample_task.start()
            self.time_series_bar.reset()
            scan_name = self.scan_container.scan_samples[scan_index].p.get_misc_parameter("scan_name",
                                                                                          "scan_name").get_in_sim_unit()
            sample_task.info.update({"scan_index": scan_index, "scan_name": scan_name})
            sample_task.update_child_info()
            try:
                def post_step(sc, time_index, t, T):
                    self.estimate_time_remaining(sc, scan_index, time_index, list(range(n_samples)), T)
                    self.scan_tree.add_time_step_to_element_tree(sc, model_index, scan_index, time_index, t)
                    self.scan_tree.write_element_tree()

                    import resource

                    rusage = [time_index, scan_index, resource.getrusage(resource.RUSAGE_SELF)]
                    if self.ruse is None:
                        self.ruse = [rusage]
                    else:
                        self.ruse.append(rusage)

                self.scan_container.t = 0

                get_sim_container_task = sample_task.start_child("build_sim_container")
                self.sim_container = self.scenario.get_sim_container(self.get_scan_sample(scan_index).p,
                                                                     model_index=model_index)
                self.sim_container.path = self.path

                self.sim_container.marker_lookup = self.marker_lookup
                get_sim_container_task.stop()
                sample_task.add_child(self.sim_container.record)

                self.apply_sample_flags(self.sim_container, scan_index)

                initialize_task = sample_task.start_child("initialize")
                if not ext_cache == "":
                    self.sim_container.set_ext_cache(ext_cache)
                self.sim_container.initialize()
                initialize_task.stop()

                sample_task.start_child("update_sim_container")
                self.update_sim_container(self.sim_container, scan_index, model_index)
                sample_task.stop_child("update_sim_container")

                sample_task.start_child("init_xdmf_files")
                self.sim_container.init_xdmf_files()
                sample_task.stop_child("init_xdmf_files")

                sample_task.start_child("pre_scan")
                self.pre_scan(self, scan_index)
                sample_task.stop_child("pre_scan")

                self.sim_container._post_step = post_step

                # todo not a permantent solution
                self.sim_container.pre_step = self.pre_step if hasattr(self,
                                                                       "pre_step") else self.sim_container.pre_step
                self.sim_container.post_step = self.post_step if hasattr(self,
                                                                         "post_step") else self.sim_container.post_step

                if self.T is None:
                    T = np.arange(0, self.N * self.dt, self.dt)
                    self.time_series_bar.total = len(T) - 1
                    self.sim_container.run(T)
                else:
                    self.time_series_bar.total = len(self.T) - 1
                    self.sim_container.run(self.T)

                sample_task.start_child("post_scan")
                self.post_scan(self, scan_index)
                sample_task.stop_child("post_scan")
                self.scan_bar.update(1)

            except Exception as e:

                warning("Scan {i} failed.".format(i=scan_index))
                critical(traceback.format_exc())
                warning("Continuing to next scan sample")
                self.scan_bar.update(1)
                sample_task.reset()
                if self.debug:
                    raise e
                else:
                    continue

            sample_task.stop()
            run_task.start_child("write_element_tree")
            self.scan_tree.write_element_tree()
            run_task.stop_child("write_element_tree")

        run_task.stop()
        self.scan_bar.close()
        self.time_series_bar.close()

        self.save_records()

    def get_records_path(self) -> str:

        records_path = os.path.join(self.path, "records")
        os.makedirs(records_path, exist_ok=True)

        return records_path

    def save_records(self) -> None:

        records = self.record.gather_records()
        records_path = self.get_records_path()

        for f in os.listdir(records_path):
            os.remove(os.path.join(records_path, f))

        ti = [r[0] for r in self.ruse]
        si = [r[1] for r in self.ruse]

        ruse = pd.DataFrame([r[2] for r in self.ruse])
        ruse.columns = [
            "ru_utime",
            "ru_stime",
            "ru_maxrss",
            "ru_ixrss",
            "ru_idrss",
            "ru_isrss",
            "ru_minflt",
            "ru_majflt",
            "ru_nswap",
            "ru_inblock",
            "ru_oublock",
            "ru_msgsnd",
            "ru_msgrcv",
            "ru_nsignals",
            "ru_nvcsw",
            "ru_nivcsw"]
        ruse["time_index"] = ti
        ruse["scan_index"] = si

        ruse.to_hdf(os.path.join(records_path, "ruse.h5"), key="df", mode="w")

        with open(os.path.join(records_path, "dump.json"), "w") as f:
            json.dump(records, f)

        np.save(os.path.join(records_path, "eta"), np.array(self.eta_estimates))

    def pre_scan(self, state_manager: 'StateManager', scan_index: int):
        pass

    def post_scan(self, state_manager: 'StateManager', scan_index: int):
        pass


# class ScanManager:
#
#     def __init__(self, path: str):
#
#         self.path: str = path
#         self.scan_folder_pattern: str = "scan_{n}/"
#         self.element_tree: ElementTree = None
#         self.scan_container: ScanContainer = None
#         self.sim_container = None
#         self.T = None
#         self.dt = 1
#         self.N = 10
#         self.compress_log_file = True
#         self.global_collections = GlobalCollections()
#         self.global_parameters = GlobalParameters()
#         self.record = ClassRecord("ScanManager")
#         self.progress_bar = None
#
#     def update_sim_container(self, sc, i) -> Dict:
#
#         scan_container = self.deserialize_from_element_tree()
#         sc.path = self.get_scan_folder(i)
#         assert i <= len(scan_container.scan_samples) - 1
#         sample = scan_container.scan_samples[i]
#
#         assert hasattr(sc, "default_sample")
#
#         sc.p.update(sc.default_sample.p)
#         sc.p.update(sample.p)
#
#         for f in sc.fields:
#             f.apply_sample(sc.default_sample)
#             f.apply_sample(sample)
#
#         for e in sc.entity_list:
#             e.p.update(sc.default_sample.p, override=True)
#             e.p.update(sample.p, override=True)
#
#         for entity_type in sc.default_sample.entity_types:
#             sc.add_entity_type(entity_type)
#
#         for entity_type in sample.entity_types:
#             sc.add_entity_type(entity_type)
#
#         return deepcopy(sample.p)


class MyScanTree:

    def __init__(self, path: str):

        self.path: str = path
        self.element_tree: ET.ElementTree = None
        self.compress_log_file: bool = True

        self.global_collections: GlobalCollections = GlobalCollections()
        self.global_parameters: GlobalParameters = GlobalParameters()

    def load_xml(self) -> None:

        """loads xml representation from file"""
        file_path = os.path.join(self.path, "log.xml")
        if os.path.exists(file_path):
            parser = etree.XMLParser(remove_blank_text=True)
            self.element_tree = etree.parse(file_path, parser)

    def get_scan_elements(self, scan_indicies = None):

        root = self.getroot()
        scans = []
        if scan_indicies is not None:
            for si in scan_indicies:
                scans = scans + (root.findall("ScanContainer/ScanSample[@scan_index = '{si}']".format(si=si)))
        else:
            scans = root.findall("ScanContainer/ScanSample")

        return scans

    def get_scan_element(self, scan_index = None):

        scan = self.getroot().find("./ScanContainer/ScanSample[@scan_index='{i}']".format(i=scan_index))

        return scan

    def get_model_elements(self, scan_indicies = None, model_indicies = None):

        scans = []
        models = []

        root = self.getroot()
        if scan_indicies is not None:
            for si in scan_indicies:
                scans = scans + (root.findall("ScanContainer/ScanSample[@scan_index = '{si}']".format(si = si)))
        else:
            scans = root.findall("ScanContainer/ScanSample")

            for scan in scans:
                if model_indicies is not None:
                    for mi in model_indicies:
                        models = models + scan.findall("Model[@model_index = '{mi}']".format(mi = mi))
                else:
                    models = models + scan.findall("Model")


        return models

    def get_model_element(self, scan_index:int, model_index:int) -> Element:

        return  self.get_scan_element(scan_index).find("Model[@model_index='{m}']".format(m=model_index))

    def create_model_element(self, scan_index:int, model_index:int) -> Element:

        model = ET.SubElement(self.get_scan_element(scan_index), "Model")
        model.set("model_index", str(model_index))

        return model

    def write_element_tree(self) -> None:

        try:
            f = os.path.join(self.path, "log.xml")
            message("writing element tree to {file}".format(file=f))
            self.element_tree.write(f, pretty_print=True)
        except Exception as e:
            message("Could not write element tree to file: {e}".format(e=e))

    def serialize_to_element_tree(self, scan_container) -> None:

        scan_folder_pattern = "scan_{n}"

        if scan_container is None:
            message("Scan Container not set. Running one sample with default parmeters")
            scan_container = ScanContainer()

        if len(scan_container.scan_samples) == 0:
            message("No Scan Sample found. Running one sample with default parmeters")
            scan_container.add_sample(ScanSample([], [], {}, scan_name="default"))

        root = ET.Element("Run")
        root.set("path", json.dumps(self.path))
        root.set("scan_folder_pattern", json.dumps(scan_folder_pattern))

        root.append(scan_container.serialize_to_xml(scan_folder_pattern))

        self.element_tree = ET.ElementTree(element=root)

        return scan_container

    def getroot(self):

        return self.element_tree.getroot()

    def deserialize_from_element_tree(self) -> ScanContainer:

        root = self.element_tree.getroot()
        # self.scan_folder_pattern = json.loads(root.get("scan_folder_pattern"))

        scan_container = ScanContainer()
        scan_container.deserialize_from_xml(root.find("ScanContainer"))

        return scan_container

    def get_field_names(self) -> List[str]:

        scans = self.element_tree.getroot().find("ScanContainer")
        names = np.unique(
            np.array([i.get("field_name") for i in scans.findall("./ScanSample/TimeSeries/Step/Fields/Field")]))
        return names

    def add_time_step_to_element_tree(self, sim_container, model_index: int, scan_index: int, time_index: int,
                                      time: float) -> None:

        """
        :param marker_paths: Dictionary of shape {marker_key :marker_name}
        """



        scan = self.get_model_element(scan_index,model_index)
        if scan is None:
            scan = self.create_model_element(scan_index,model_index)

        if scan.find("TimeSeries") is not None:
            time_series = scan.find("TimeSeries")
        else:
            time_series = ET.SubElement(scan, "TimeSeries")

        if self.compress_log_file:
            path = "./scan_{i}/timestep_logs/".format(i=scan_index)
            file_name = "step_{mi}_{si}_{ti}.xml".format(si=scan_index, ti=time_index, mi= model_index)
            os.makedirs(os.path.join(self.path, path), exist_ok=True)
            path = os.path.join(path, file_name)
            substitute_element = ET.SubElement(time_series, "Step")
            substitute_element.set("path", path)
            step = ET.Element("Step")


        else:

            step = ET.SubElement(time_series, "Step")
        step.set("time_index", str(time_index))
        step.set("time", str(time))

        fields = ET.SubElement(step, "Fields")

        for field in sim_container.fields:
            field_element = field.get_result_element(time_index, scan_index, sim_container.get_current_path(),
                                                     sim_container.markers, sim_container.marker_lookup)
            fields.insert(0, field_element)

        lookup_table_element = ET.SubElement(step, "MarkerLookupTable")
        for k, v in sim_container.marker_lookup.items():
            lookup_element = ET.SubElement(lookup_table_element, "MarkerLookup")
            lookup_element.set("key", str(k))
            lookup_element.set("value", str(v))

        cells = ET.SubElement(step, "Cells")
        for c in sim_container.entity_list:
            if isinstance(c, Cell):
                cell = ET.SubElement(cells, "Cell")

                cell.set("x", str(c.center[0] + c.offset[0]))
                cell.set("y", str(c.center[1] + c.offset[1]))
                cell.set("z", str(c.center[2] + c.offset[2]))
                cell.set("name", str(c.name))
                cell.set("entity_id", str(c.id))
                cell.set("type_name", str(c.type_name))
                if self.compress_log_file:
                    cell.append(c.p.serialize_to_xml(
                        global_collections=self.global_collections,
                        global_parameters=self.global_parameters
                    ))
                else:
                    cell.append(c.p.serialize_to_xml())

        old_e = time_series.find("GlobalCollections")
        if not old_e is None:
            time_series.remove(old_e)

        old_e = time_series.find("GlobalParameters")
        if not old_e is None:
            time_series.remove(old_e)

        time_series.insert(0, self.global_collections.serialize_to_xml())
        time_series.insert(0, self.global_parameters.serialize_to_xml())

        if self.compress_log_file:
            tree = ET.ElementTree(step)
            tree.write(os.path.join(self.path, path), pretty_print=True)

    def rebuild_timesteps(self, element: Element) -> None:
        """


        :param element: lxml element with TimeSeries/Step sub elements
        :return:
        """
        for step in element.findall(".//TimeSeries/Step"):
            path = step.get("path")
            if not path is None:
                et = ET.parse(os.path.join(self.path, path))
                step.getparent().replace(step, et.getroot())

        return element


def target(mp_input: Tuple[int, str, List[int], str]) -> List[pd.DataFrame]:

    try:
        n_scans, scan, time_indices,path = mp_input
        scan = ET.fromstring(scan)

        result = []
        for model in scan.findall("Model"):
            for step in model.findall("TimeSeries/Step"):

                # if "path" in step.attrib.keys():
                #     step = ET.parse(os.path.join(path, step.get("path")))
                #     step = step.getroot()

                time_index = step.get("time_index")
                if (not time_indices is None) and not (int(time_index) in time_indices):
                    continue
                time = step.get("time")

                message("State Manager: Computing timestep {n} for scan {scan_n}".format(n=time_index, scan_n=n_scans))

                for cell in step.findall("Cells/Cell"):

                    parameter_set = ParameterSet.deserialize_from_xml(cell.find("ParameterSet"))

                    p_temp = parameter_set.get_as_dictionary()
                    p = {}
                    import numbers
                    from typing import List

                    for k, v in p_temp.items():
                        if (isinstance(v, str) or isinstance(v, numbers.Number)):
                            p[k] = v
                        elif (isinstance(v, List)):
                            p[k] = pd.Series(v)
                            # pass
                    p["time_index"] = time_index
                    p["time"] = time
                    p["scan_index"] = scan.get("scan_index")
                    p["model_index"] = model.get("model_index")
                    p["x"] = float(cell.get("x"))
                    p["y"] = float(cell.get("y"))
                    p["z"] = float(cell.get("z"))
                    p["id"] = cell.get("entity_id")
                    p["type_name"] = cell.get("type_name")

                    result.append(p)
        return result
    except:
        import traceback
        print(traceback.format_exc())
        return []
