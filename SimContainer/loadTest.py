#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 10:27:44 2019

@author: lukas
"""
import dolfin as dlf


mesh = dlf.Mesh()
with dlf.XDMFFile(dlf.MPI.comm_world,"./cache/meshCache_il2.xdmf") as f:
    f.read(mesh)