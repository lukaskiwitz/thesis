#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 14:34:19 2019

@author: Lukas Kiwitz
"""
import json
import time
from copy import deepcopy
from typing import Dict

import lxml.etree as ET
import mpi4py.MPI as MPI
import numpy as np
import pandas as pd

from Entity import Cell
from ParameterSet import ParameterSet
from ScanContainer import ScanContainer, ScanSample
from my_debug import message, total_time

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
        self.scan_folder_pattern = "scan_{n}/"
        self.element_tree = None
        self.scan_container = None
        self.sim_container = None
        self.T = [0]
        self.dt = 1


    def load_xml(self):
        """loads xml representation from file"""
        self.element_tree = ET.parse("{p}log.scan".format(p=self.path))

    def get_scan_folder(self, n):
        return self.path + self.scan_folder_pattern.format(n=n)

    def write_element_tree(self):
        if not rank == 0:
            # message("not rank 0")
            return None
        f = self.path + "log.scan"
        message("writing element tree to {file}".format(file=f))
        self.element_tree.write(f, pretty_print=True)

    def serialize_to_element_tree(self):


        if not self.scan_container:
            message("Scan Container not set. Running one sample with default parmeters")
            self.scan_container = ScanContainer()

        if len(self.scan_container.scan_samples) == 0:
            message("No Scan Sample found. Running one sample with default parmeters")
            self.scan_container.add_sample(ScanSample([],[],{}))


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

    def add_time_step_to_element_tree(self, sc, scan_index: int, time_step: int, time: float, result_path: Dict):
        scan = self.element_tree.getroot().find("ScanContainer/ScanSample[@scan_index='{s}']".format(s=scan_index))
        if scan.find("TimeSeries"):
            time_series = scan.find("TimeSeries")
        else:
            time_series = ET.SubElement(scan, "TimeSeries")

        step = ET.SubElement(time_series, "Step")
        step.set("time_index", str(time_step))
        step.set("time", str(time))


        cells = ET.SubElement(step, "Cells")
        for c in sc.entity_list:
            if isinstance(c, Cell):
                cell = ET.SubElement(cells, "Cell")

                cell.set("x", str(c.center[0]))
                cell.set("y", str(c.center[1]))
                cell.set("z", str(c.center[2]))
                cell.set("name", str(c.name))
                cell.set("entity_id", str(c.id))
                cell.set("type_name", str(c.type_name))
                cell.append(c.p.serialize_to_xml())

        fields = ET.SubElement(step,"Fields")
        for field_name, (distplot, sol, field_index) in result_path.items():
            field = time_series.find("Field[@field_name='{n}']".format(n=field_name))
            if not field:
                field = ET.SubElement(fields, "Field")
                field.set("mesh_path",sc.fields[field_index].get_mesh_path())
                field.set("field_name", field_name)
                field.set("dist_plot_path", "scan_{i}/".format(i = scan_index) + distplot)
                field.set("solution_path", "scan_{i}/".format(i = scan_index) + sol)

    def get_field_names(self):

        scans = self.element_tree.getroot().find("ScanContainer")
        names = np.unique(np.array([i.get("field_name") for i in scans.findall("./ScanSample/TimeSeries/Step/Fields/Field")]))
        return names

    def get_cell_ts_data_frame(self, **kwargs):

        root = self.element_tree.getroot()
        result = []
        if "scan_index" in kwargs:

            scans = root.findall("./ScanContainer/ScanSample[@scan_index='{i}']".format(i=kwargs["scan_index"]))

        else:
            scans = root.findall("ScanContainer/ScanSample")

        for n_scans, scan in enumerate(scans):

            message("State Manager: Collecting cell dataframe for {n} of {total} scans".format(n=n_scans, total=len(scans)))

            #for field in scan.findall("TimeSeries/Field"):
            for step in scan.findall("TimeSeries/Step"):
                time_index = step.get("time_index")
                time = step.get("time")

                message("State Manager: Computing timestep {n} for scan {scan_n}".format(n=time_index, scan_n=n_scans))

                for cell in step.findall("Cells/Cell"):

                    parameter_set = ParameterSet("dummy",[])
                    parameter_set.deserialize_from_xml(cell.find("ParameterSet"))

                    p_temp = parameter_set.get_as_dictionary()
                    p = {}
                    import numbers
                    for k,v in p_temp.items():
                        if  (isinstance(v,str) or isinstance(v,numbers.Number)):
                           p[k] = v

                    p["time_index"] = time_index
                    p["time"] = time
                    p["scan_index"] = scan.get("scan_index")

                    p["x"] = float(cell.get("x"))
                    p["y"] = float(cell.get("y"))
                    p["z"] = float(cell.get("z"))
                    p["id"] = cell.get("entity_id")
                    p["type_name"] = cell.get("type_name")

                    result.append(p)
        return pd.DataFrame(result)

    def update_sim_container(self, sc, i) -> Dict:

        scan_container = self.deserialize_from_element_tree()
        sc.path = self.get_scan_folder(i)
        assert i <= len(scan_container.scan_samples) - 1
        sample = scan_container.scan_samples[i]

        sc.p.update(sample.p)

        for f in sc.fields:
            f.apply_sample(sample)

        for e in sc.entity_list:
            e.p.update(sample.p)

        for entity_type in sample.entity_types:
            sc.add_entity_type(entity_type)

        return deepcopy(sample.p)

    def run(self):


        self.serialize_to_element_tree()

        for scan_index in range(len(self.scan_container.scan_samples)):

            def post_step(sc, time_index, t, T):
                result_paths = sc.save_fields(time_index)
                self.add_time_step_to_element_tree(sc, scan_index, time_index, t, result_paths)

            start = time.process_time()

            self.update_sim_container(self.sim_container, scan_index)
            self.sim_container.init_xdmf_files()
            self.pre_scan(self, scan_index)

            self.sim_container._post_step = post_step

            self.sim_container.run(self.T, self.dt)

            self.post_scan(self, scan_index)

            end = time.process_time()
            total_time(end - start, pre="Total time ", post=" for scan sample {n}".format(n=scan_index))
        self.write_element_tree()

    def pre_scan(self, state_manager, scan_index):
        pass

    def post_scan(self, state_manager, scan_index):
        pass