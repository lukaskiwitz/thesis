#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 14:34:19 2019

@author: Lukas Kiwitz
"""
import json
import os
import time
from copy import deepcopy
from typing import Dict

import lxml.etree as ET
import mpi4py.MPI as MPI
import numpy as np
import pandas as pd
import traceback
import scipy

from thesis.main.Entity import Cell
from thesis.main.ParameterSet import ParameterSet, GlobalCollections, GlobalParameters
from thesis.main.ScanContainer import ScanContainer, ScanSample
from thesis.main.my_debug import message, warning, critical
from tqdm import tqdm
from thesis.main.TaskRecord import TaskRecord, ClassRecord
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


def outputParse(v):
    if type(v) == np.ndarray:
        return json.dumps(v.tolist())
    else:
        return json.dumps(v)

class StateManager:
    """
    Used to store Simulation state as xml file. Either to save the state of a simulation to file or for usage in post processing
    """

    def __init__(self, path):
        """Docstring for constructor

        :param path: This is a test
        :return None
        """
        self.path = path
        self.ruse = []

        self.scan_folder_pattern = "scan_{n}/"
        self.element_tree = None
        self.scan_container = None
        self.sim_container = None
        self.scenario = None
        self.T = None
        self.dt = 1
        self.N = 10
        self.compress_log_file = True
        self.global_collections = GlobalCollections()
        self.global_parameters = GlobalParameters()
        self.record = ClassRecord("StateManager")
        self.progress_bar = None
        self.eta_estimates = []
        self.marker_lookup = {}

    def load_xml(self):
        """loads xml representation from file"""
        self.element_tree = ET.parse("{p}log.scan".format(p=self.path))

    def get_scan_folder(self, n):
        return self.path + self.scan_folder_pattern.format(n=n)

    def write_element_tree(self):

        if rank == 0:
            try:
                f = self.path + "log.scan"
                message("writing element tree to {file}".format(file=f))
                self.element_tree.write(f, pretty_print=True)
            except Exception as e:
                message("Could not write element tree to file: {e}".format(e=e))


    def serialize_to_element_tree(self):


        if not self.scan_container:
            message("Scan Container not set. Running one sample with default parmeters")
            self.scan_container = ScanContainer()

        if len(self.scan_container.scan_samples) == 0:
            message("No Scan Sample found. Running one sample with default parmeters")
            self.scan_container.add_sample(ScanSample([],[],{}, scan_name="default"))


        root = ET.Element("Run")
        root.set("path", json.dumps(self.path))
        root.set("scan_folder_pattern", json.dumps(self.scan_folder_pattern))

        root.append(self.scan_container.serialize_to_xml(self.scan_folder_pattern))

        self.element_tree = ET.ElementTree(element=root)


    def deserialize_from_element_tree(self):

        root = self.element_tree.getroot()
        self.scan_folder_pattern = json.loads(root.get("scan_folder_pattern"))

        scan = ScanContainer()
        scan.deserialize_from_xml(root.find("ScanContainer"))

        return scan

    def _addCellDump(self, sc, i):
        """
        i Scan index in xml file
        """
        run = self.element_tree.getroot()
        cellDump = ET.SubElement(run, "cellDump")

        for f in [sc.fields[0]]:
            field = ET.SubElement(cellDump, "field")
            field.set("name", str(f.field_name))
            extCache = ET.SubElement(field, "ext_cache")
            extCache.text = f.ext_cache
            subdomains = ET.SubElement(field, "subdomains")
            subdomains.text = sc.ext_boundary_markers
            for n in f.registered_entities:
                cell = ET.SubElement(field, "cell")
                patch = ET.SubElement(cell, "patch")
                center = ET.SubElement(cell, "center")

                patch.text = str(n["patch"])
                center.text = json.dumps(list(n["entity"].center))

    def _loadCellDump(self):
        self.cellDump = self.element_tree.getroot().find("/cellDump")

    def add_time_step_to_element_tree(self, sc, scan_index: int, time_step: int, time: float, result_path: Dict, marker_paths: Dict):
        scan = self.element_tree.getroot().find("ScanContainer/ScanSample[@scan_index='{s}']".format(s=scan_index))
        if scan.find("TimeSeries") is not None:
            time_series = scan.find("TimeSeries")
        else:
            time_series = ET.SubElement(scan, "TimeSeries")

        if self.compress_log_file:
            path = "./scan_{i}/timestep_logs/".format(i = scan_index)
            file_name = "step_{si}_{ti}.xml".format(si = scan_index, ti = time_step)
            os.makedirs(os.path.join(self.path, path),exist_ok=True)
            path = os.path.join(path, file_name)
            substitute_element = ET.SubElement(time_series, "Step")
            substitute_element.set("path", path)
            step = ET.Element("Step")


        else:
            step = ET.SubElement(time_series, "Step")
        step.set("time_index", str(time_step))
        step.set("time", str(time))

        lookup_table_element = ET.SubElement(step, "MarkerLookupTable")
        for k,v in self.sim_container.marker_lookup.items():
            lookup_element = ET.SubElement(lookup_table_element, "MarkerLookup")
            lookup_element.set("key",str(k))
            lookup_element.set("value",str(v))


        for marker_key, marker_path in marker_paths.items():
            marker_element = ET.SubElement(step,"Marker")
            marker_element.set("path",marker_path)
            marker_element.set("marker_key",marker_key)


        cells = ET.SubElement(step, "Cells")
        for c in sc.entity_list:
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
                        global_collections = self.global_collections,
                        global_parameters = self.global_parameters
                    ))
                else:
                    cell.append(c.p.serialize_to_xml())



        old_e = time_series.find("GlobalCollections")
        if not old_e is None:
            time_series.remove(old_e)

        old_e = time_series.find("GlobalParameters")
        if not old_e is None:
            time_series.remove(old_e)

        time_series.insert(0,self.global_collections.serialize_to_xml())
        time_series.insert(0, self.global_parameters.serialize_to_xml())

        fields = ET.SubElement(step,"Fields")

        for field_name, (distplot, sol, field_index) in result_path.items():

            d = os.path.split(os.path.split(sc.path)[0])[1]
            field = time_series.find("Field[@field_name='{n}']".format(n=field_name))

            if not field:
                field = ET.SubElement(fields, "Field")

            dynamic_mesh =  sc.fields[field_index].moving_mesh or sc.fields[field_index].remesh or sc.fields[field_index].ale
            if dynamic_mesh:
                field.set("mesh_path",os.path.join(d,sc.fields[field_index].get_mesh_path(time_step-1, local=True)))
            else:
                field.set("mesh_path", sc.fields[field_index].get_mesh_path(time_step-1, local=True))

            field.set("dynamic_mesh",str(dynamic_mesh))
            field.set("field_name", field_name)

            if sol is None:
                field.set("success",str(False))
                continue
            else:
                field.set("success", str(True))
                field.set("dist_plot_path", "scan_{i}/".format(i = scan_index) + distplot)
                field.set("solution_path", "scan_{i}/".format(i = scan_index) + sol)

        if self.compress_log_file:
            tree = ET.ElementTree(step)
            tree.write(os.path.join(self.path, path),pretty_print = True)

    def get_field_names(self):

        scans = self.element_tree.getroot().find("ScanContainer")
        names = np.unique(np.array([i.get("field_name") for i in scans.findall("./ScanSample/TimeSeries/Step/Fields/Field")]))
        return names

    def rebuild_tree(self,scan):
        for step in scan.findall("TimeSeries/Step"):
            path = step.get("path")
            if not path is None:
                et = ET.parse(os.path.join(self.path, path))
                step.getparent().replace(step, et.getroot())

    def get_cell_ts_data_frame(self, time_indices = None,n_processes = 1, **kwargs):


        def init():
            global path
            path = self.path


        root = self.element_tree.getroot()
        result = []
        if "scan_index" in kwargs:

            scans = root.findall("./ScanContainer/ScanSample[@scan_index='{i}']".format(i=kwargs["scan_index"]))

        else:
            scans = root.findall("ScanContainer/ScanSample")


        import multiprocessing as mp

        scatter_list = []
        for n_scans, scan in enumerate(scans):

            # for field in scan.findall("TimeSeries/Field"):
            self.rebuild_tree(scan)
            scatter_list.append((n_scans,ET.tostring(scan),time_indices))

        n_processes = n_processes if n_processes <= os.cpu_count() else os.cpu_count()

        message("State Manager: Distributing {total} scans to {n_processes} processes".format(n_processes=n_processes,total = len(scans)))
        with mp.Pool(processes=n_processes,initializer = init) as p:
            result = p.map(target, scatter_list)
            while [] in result:
                result.remove([])

            result = list(np.array(result).ravel())

        return pd.DataFrame(result)

    def get_scan_sample(self, i):

        scan_container = self.deserialize_from_element_tree()
        assert i <= len(scan_container.scan_samples) - 1
        sample = scan_container.scan_samples[i]
        return sample

    def update_sim_container(self, sc, i) -> Dict:

        sc.path = self.get_scan_folder(i)
        sample = self.get_scan_sample(i)

        assert hasattr(sc,"default_sample")

        sc.p.update(sc.default_sample.p)
        sc.p.update(sample.p)

        for f in sc.fields:
            f.apply_sample(sc.default_sample)
            f.apply_sample(sample)

            if sample.dynamic_mesh:
                f.moving_mesh = False
                f.remesh = True

        for e in sc.entity_list:
            e.p.update(sc.default_sample.p, override = True)
            e.p.update(sample.p, override = True)

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

    def estimate_time_remaining(self, sc, scan_index, time_index,S,T):


        sample = self.record.child_tasks["run"].child_tasks["scan_sample"]
        time_step = sample.child_tasks["SimContainer"].child_tasks["run"].child_tasks["step"]
        step_durations = [i["end"] - i["start"] for i in time_step.history]
        mean_step_duration = scipy.stats.gmean(step_durations)

        n_steps = (len(T) - time_index)
        n_scans = len(S) - scan_index

        time_series_eta = time.time() + n_steps * mean_step_duration
        scan_eta = time.time() + (n_scans * len(T) * mean_step_duration + n_steps * mean_step_duration)
        self.eta_estimates.append([time.time(), scan_eta])


        np.save(os.path.join(self.get_records_path(), "eta"), np.array(self.eta_estimates))

        scan_eta = time.strftime("%H:%M:%S %m/%d/%y",time.localtime(scan_eta))
        time_series_eta = time.strftime("%H:%M:%S",time.localtime(time_series_eta))



        self.time_series_bar.update(1)
        self.scan_bar.postfix = "ETA: {eta}".format(eta = scan_eta)
        self.time_series_bar.postfix = "ETA: {eta}".format(eta=time_series_eta)
        # message("ETA: {eta}".format(eta = eta))

    def run(self, ext_cache = ""):

        run_task = self.record.start_child("run")
        sample_task = run_task.start_child("scan_sample")

        self.serialize_to_element_tree()

        n_samples = len(self.scan_container.scan_samples)
        self.scan_bar = tqdm(desc="Parameter Scan", total=n_samples, dynamic_ncols=True)
        self.time_series_bar = tqdm(desc="Time Series   ", initial=0, dynamic_ncols=True)
        self.write_element_tree()


        for scan_index in range(n_samples):
            if not sample_task.running:
                sample_task.start()
            self.time_series_bar.reset()
            scan_name = self.scan_container.scan_samples[scan_index].p.get_misc_parameter("scan_name","scan_name").get_in_sim_unit()
            sample_task.info.update({"scan_index": scan_index, "scan_name":scan_name})
            sample_task.update_child_info()
            try:
                def post_step(sc, time_index, t, T):
                    self.estimate_time_remaining(sc, scan_index, time_index, range(n_samples), T)
                    result_paths = sc.save_fields(int(time_index - 1))
                    marker_paths = sc.save_markers(time_index)
                    self.add_time_step_to_element_tree(sc, scan_index, time_index, t, result_paths, marker_paths)
                    self.write_element_tree()

                    import resource

                    rusage = [time_index, scan_index,resource.getrusage(resource.RUSAGE_SELF)]
                    if self.ruse is None:
                        self.ruse = [rusage]
                    else:
                        self.ruse.append(rusage)

                self.scan_container.t = 0

                get_sim_container_task = sample_task.start_child("build_sim_container")
                self.sim_container = self.scenario.get_sim_container(self.get_scan_sample(scan_index).p)
                self.sim_container.path = self.path
                self.sim_container.marker_lookup = self.marker_lookup
                get_sim_container_task.stop()
                sample_task.add_child(self.sim_container.record)

                initialize_task = sample_task.start_child("initialize")
                if not ext_cache == "":
                    self.sim_container.set_ext_cache(ext_cache)
                    self.sim_container.initialize(load_subdomain=True, file_name="mesh")
                else:
                    self.sim_container.initialize()
                initialize_task.stop()

                sample_task.start_child("update_sim_container")
                self.update_sim_container(self.sim_container, scan_index)
                sample_task.stop_child("update_sim_container")

                sample_task.start_child("init_xdmf_files")
                self.sim_container.init_xdmf_files()
                sample_task.stop_child("init_xdmf_files")

                sample_task.start_child("pre_scan")
                self.pre_scan(self, scan_index)
                sample_task.stop_child("pre_scan")

                self.sim_container._post_step = post_step

                #todo not a permantent solution
                self.sim_container.pre_step = self.pre_step if hasattr(self,"pre_step") else self.sim_container.pre_step
                self.sim_container.post_step = self.post_step if hasattr(self,"post_step") else self.sim_container.post_step



                if self.T is None:
                    T = np.arange(0, self.N * self.dt, self.dt)
                    self.time_series_bar.total = len(T)-1
                    self.sim_container.run(T)
                else:
                    self.time_series_bar.total = len(self.T)-1
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


            sample_task.stop()
            run_task.start_child("write_element_tree")
            self.write_element_tree()
            run_task.stop_child("write_element_tree")


        run_task.stop()
        self.scan_bar.close()
        self.time_series_bar.close()

        self.save_records()

    def get_records_path(self):

        records_path = os.path.join(self.path, "records")
        os.makedirs(records_path, exist_ok=True)

        return records_path

    def save_records(self):

        records = self.record.gather_records()
        records_path = self.get_records_path()


        for f in os.listdir(records_path):
            os.remove(os.path.join(records_path,f))

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

        ruse.to_hdf(os.path.join(records_path, "ruse.h5"),key = "df", mode = "w")

        with open(os.path.join(records_path, "dump.json"), "w") as f:
            json.dump(records, f)

        np.save(os.path.join(records_path, "eta"), np.array(self.eta_estimates))

    def pre_scan(self, state_manager, scan_index):
        pass

    def post_scan(self, state_manager, scan_index):
        pass


def target(mp_input):

    try:
        n_scans, scan, time_indices = mp_input
        scan = ET.fromstring(scan)

        result = []

        for step in scan.findall("TimeSeries/Step"):

            if "path" in step.attrib.keys():
                step = ET.parse(os.path.join(path, step.get("path")))
                step = step.getroot()

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

class ScanManager:


    def __init__(self, path):

        self.path = path
        self.scan_folder_pattern = "scan_{n}/"
        self.element_tree = None
        self.scan_container = None
        self.sim_container = None
        self.T = None
        self.dt = 1
        self.N = 10
        self.compress_log_file = True
        self.global_collections = GlobalCollections()
        self.global_parameters = GlobalParameters()
        self.record = ClassRecord("ScanManager")
        self.progress_bar = None

    def update_sim_container(self, sc, i) -> Dict:


        scan_container = self.deserialize_from_element_tree()
        sc.path = self.get_scan_folder(i)
        assert i <= len(scan_container.scan_samples) - 1
        sample = scan_container.scan_samples[i]

        assert hasattr(sc,"default_sample")

        sc.p.update(sc.default_sample.p)
        sc.p.update(sample.p)

        for f in sc.fields:
            f.apply_sample(sc.default_sample)
            f.apply_sample(sample)

        for e in sc.entity_list:
            e.p.update(sc.default_sample.p, override = True)
            e.p.update(sample.p, override = True)

        for entity_type in sc.default_sample.entity_types:
            sc.add_entity_type(entity_type)

        for entity_type in sample.entity_types:
            sc.add_entity_type(entity_type)


        return deepcopy(sample.p)