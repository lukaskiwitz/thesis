#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 14:34:19 2019

@author: kiwitz
"""
import lxml.etree as ET
import numpy as np
import json
import itertools
from copy import deepcopy

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
    def __init__(self, path):
        self.path = path
        self.scanFolderPattern = "{path}scan_{n}/"
        self.elementTree = None

    def loadXML(self):
        self.elementTree = ET.parse("{p}log.scan".format(p=self.path))

    def getScanFolder(self, n):
        return self.scanFolderPattern.format(n=n, path=self.path)

    def writeElementTree(self):
        if not rank == 0:
            print("not rank 0")
            return None
        f = self.path+"log.scan"
        print("writing element tree to {file}".format(file=f))
        self.elementTree.write(f, pretty_print=True)

    def scan_log(self, scan, p):
        root = ET.Element("run")
        scans = ET.SubElement(root, "scans")

        for n, s in enumerate(scan):
            scan = ET.SubElement(scans, "scan")
            scan.set("i", str(n))
            path = ET.SubElement(scan, "path")
            path.text = self.getScanFolder(n)

            parameters = ET.SubElement(scan, "parameters")

            constant = ET.SubElement(parameters, "constant")
            for k, v in p.items():
                par = ET.SubElement(constant, "parameter")
                par.set("type", str(type(v)))
                par.set("name", k)
                par.text = outputParse(v)
            dynamic = ET.SubElement(parameters, "dynamic")
            for k, v in s.items():
                par = ET.SubElement(dynamic, "parameter")
                par.set("type", str(type(v)))
                par.set("name", k)
                par.text = outputParse(v)
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
            extCache = ET.SubElement(field, "extCache")
            extCache.text = f.extCache
            subdomains = ET.SubElement(field, "subdomains")
            subdomains.text = sc.extBoundaryMarkers
            for n in f.registeredEntities:
                cell = ET.SubElement(field, "cell")
                patch = ET.SubElement(cell, "patch")
                center = ET.SubElement(cell, "center")

                patch.text = str(n["patch"])
                center.text = json.dumps(list(n["entity"].center))

#            dumpList.append({"patch":i["patch"],"center":i["entity"].center})
    def loadCellDump(self):
        cellDump = self.elementTree.getroot().find("/cellDump")

    def addTimeStep(self, i, n, t, fieldName="field", displot="", sol=""):
        scans = self.elementTree.getroot().find("scans")
        scan = scans.find("./scan[@i='{i}']".format(i=i))
#        scan  = next(itertools.islice(scans.iter(tag="scan"),i,i+1,1))

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
            f.outerDomain.p.update(s)
        return deepcopy(p)
