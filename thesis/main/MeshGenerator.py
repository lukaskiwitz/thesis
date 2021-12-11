#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 12:06:41 2019

@author: Lukas Kiwitz
"""

import logging
import os

import dolfin as dlf
# from dolfin import *
import fenics as fcs
import meshio
import numpy as np
# import h5py
import pygmsh

import thesis.main.Entity as Entity
from thesis.main.my_debug import message

module_logger = logging.getLogger('main.MeshGenerator')

class DomainTypeError(Exception):
    def __init__(self, text):
        self.text = text

    def __str__(self):
        return "Unknown Domain Type {t}".format(t=self.text)


class MeshGenerator:
    outerDomain = None
    entityList = []

    dim = 2

    def __init__(self, **kwargs):
        if "outer_domain" in kwargs:
            self.outerDomain = kwargs["outer_domain"]

        self.logger = logging.getLogger('main.MeshGenerator.MeshGenerator')

    def meshGen(self, p, mesh_path, subdomain_path, load_mesh=False, load_subdomain=False, ):

        os.makedirs(os.path.dirname(mesh_path), exist_ok=True)
        os.makedirs(os.path.dirname(subdomain_path), exist_ok=True)

        geom = pygmsh.opencascade.Geometry(
            characteristic_length_min=p.get_misc_parameter("min_char_length", "numeric").get_in_sim_unit(),
            characteristic_length_max=p.get_misc_parameter("max_char_length", "numeric").get_in_sim_unit()
        )

        if isinstance(self.outerDomain, Entity.DomainCube):
            p1 = self.outerDomain.p1
            p2 = self.outerDomain.p2
            domain = geom.add_box(p1, np.array(p2) - np.array(p1))
        elif isinstance(self.outerDomain, Entity.DomainSphere):
            c = self.outerDomain.center
            r = self.outerDomain.radius
            domain = geom.add_ball(c, r)
        else:
            raise DomainTypeError(type(self.outerDomain))
        entities = []

        for i in self.entityList:
            r = i["entity"].radius
            c = i["entity"].center
            ball = geom.add_ball(c, r)
            geom.add_physical(ball)
            entities.append(ball)
        if not load_mesh or not os.path.isfile(mesh_path):
            if len(entities) > 0:
                geom.boolean_difference([domain], entities)
            mesh = pygmsh.generate_mesh(geom)
            meshio.write(mesh_path, meshio.Mesh(points=mesh.points, cells={"tetra": mesh.cells["tetra"]}))

        else:
            message("loading Mesh from: " + mesh_path, self.logger)
        mesh = dlf.Mesh()

        with dlf.XDMFFile(dlf.MPI.comm_world, mesh_path) as f:
            f.read(mesh)
        message("mesh loaded", self.logger)

        if load_subdomain and os.path.isfile(subdomain_path):
            message("loading subdomain from " + subdomain_path, self.logger)

            boundary_markers = fcs.MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
            with fcs.HDF5File(fcs.MPI.comm_world, subdomain_path, "r") as f:
                f.read(boundary_markers, "/boundaries")
            for i, o in enumerate(self.entityList):
                domain_patches = self.outerDomain.get_subdomains()
                a = domain_patches[-1]["patch"] + 1 if len(domain_patches) > 0 else 1
                o["entity"].get_compiled_subdomain().mark(boundary_markers, i + a)
                o["patch"] = i + a

        else:
            boundary_markers = fcs.MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
            boundary_markers.set_all(0)
            message("boundaries marked", self.logger)

            for i in self.outerDomain.get_subdomains():
                i["entity"].get_subdomain().mark(boundary_markers, i["patch"])
            message("outer domain set", self.logger)
            for i, o in enumerate(self.entityList):
                domain_patches = self.outerDomain.get_subdomains()
                a = domain_patches[-1]["patch"] + 1 if len(domain_patches) > 0 else 1
                o["entity"].get_compiled_subdomain().mark(boundary_markers, i + a)
                o["patch"] = i + a
            message("loop complete", self.logger)
            if load_subdomain:
                with fcs.HDF5File(fcs.MPI.comm_world, subdomain_path, "w") as f:
                    f.write(boundary_markers, "/boundaries")
        return mesh, boundary_markers
