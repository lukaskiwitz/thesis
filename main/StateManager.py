#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 14:34:19 2019

@author: Lukas Kiwitz
"""
import lxml.etree as ET
import numpy as np
import json
import itertools
from copy import deepcopy
import pandas as pd

import mpi4py.MPI as MPI

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
        self.scanFolderPattern = "{path}scan_{n}/"
        self.elementTree = None

    def loadXML(self):
        """loads xml representation from file"""
        self.elementTree = ET.parse("{p}log.scan".format(p=self.path))

    def getScanFolder(self, n):
        return self.scanFolderPattern.format(n=n, path=self.path)

    def writeElementTree(self):
        if not rank == 0:
            # print("not rank 0")
            return None
        f = self.path+"log.scan"
        print("writing element tree to {file}".format(file=f))
        self.elementTree.write(f, pretty_print=True)

    def output_parameter_dict(self,p,root_name):
        root = ET.Element(root_name)
        for k, v in p.items():
            par = ET.SubElement(root, "parameter")
            par.set("type", str(type(v)))
            par.set("name", k)
            par.text = outputParse(v)
        return root

    def scan_log(self, scan, p):
        root = ET.Element("run")
        scans = ET.SubElement(root, "scans")

        for n, s in enumerate(scan):
            scan = ET.SubElement(scans, "scan")
            scan.set("i", str(n))
            path = ET.SubElement(scan, "path")
            path.text = self.getScanFolder(n)

            parameters: ET.Element = ET.SubElement(scan, "parameters")

            parameters.append(self.output_parameter_dict(p,"constant"))
            parameters.append(self.output_parameter_dict(s, "dynamic"))

            number = ET.SubElement(scan, "number")
            number.text = str(n)
            timeSeries = ET.SubElement(scan, "timeSeries")

        self.elementTree = ET.ElementTree(element=root)
#        self.writeElementTree()

    def addCellDump(self, sc, i):
        """
        i Scan index in xml file
        """
        run = self.elementTree.getroot()
        cellDump = ET.SubElement(run, "cellDump")

        for f in [sc.fields[0]]:
            field = ET.SubElement(cellDump, "field")
            field.set("name", str(f.fieldName))
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

    def loadCellDump(self):
        self.cellDump = self.elementTree.getroot().find("/cellDump")

    def addTimeStep(self, i, n, t, fieldName="field", displot="", sol="",cell_list=[]):
        scans = self.elementTree.getroot().find("scans")
        scan = scans.find("./scan[@i='{i}']".format(i=i))

        timeSeries = scan.find("timeSeries")
        field = timeSeries.find("./field[@name='{f}']".format(f=fieldName))

        if field == None:
            field = ET.SubElement(timeSeries, "field")
            field.set("name", str(fieldName))
        step = ET.SubElement(field, "step")
        step.set("n", str(n))
        step.set("t", str(t))
        distPlotPath = ET.SubElement(step, "distPlotPath")
        distPlotPath.text = displot

        solutionPath = ET.SubElement(step, "solutionPath")
        solutionPath.text = sol

        cells = ET.SubElement(step,"cells")
        for c in cell_list:
            cell = ET.SubElement(cells,"cell")
            cell.append(self.output_parameter_dict(c, "properties"))


        self.writeElementTree()

    def getParametersFromElement(self, element):
        p = {}
        for i in element.iter(tag="parameter"):
            t = i.get("type")
            if t == "<class 'numpy.ndarray'>":
                l = json.loads(i.text)
                p[i.get("name")] = np.array(l)
            else:
                p[i.get("name")] = json.loads(i.text)
        return p
    def get_cell_ts_data_frame(self):
        root = self.elementTree.getroot()
        # result = pd.DataFrame()
        result = []
        for scan in root.findall("scans/scan"):
            for field in scan.findall("timeSeries/field"):
                for step in field.findall("step"):
                    t = step.get("t")
                    n = step.get("n")
                    print(t)
                    for cell in step.findall("cells/cell"):
                        p = self.getParametersFromElement(cell.find("properties"))
                        p["t"] = t
                        p["n"] = n


                        result.append(p)
        return pd.DataFrame(result)

    def updateSimContainer(self, sc, i):
        scans = self.elementTree.getroot().find("scans")
        scan = next(itertools.islice(scans.iter(tag="scan"), i, i+1, 1))

        p = self.getParametersFromElement(
            scan.find("parameters").find("constant"))
        s = self.getParametersFromElement(
            scan.find("parameters").find("dynamic"))
        p.update(s)
        for f in sc.fields:
            f.p = deepcopy(p)
            f.solver.p = deepcopy(p)
            f.outer_domain.p.update(s)
        return deepcopy(p)
