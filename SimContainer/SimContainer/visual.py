#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 10:09:06 2019

@author: kiwitz
#"""
from distPlot_dump_test import distPlot_dump
from distPlot import distPlot
import mpi4py.MPI as MPI
import seaborn as sns

path = "/extra/kiwitz/results/"
extCache = "/extra/kiwitz/extCache/"
fields = ["il2","il6"]

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

l = ["secHigh_dirichlet_fractionLow",
"secLow_dirichlet_fractionHigh",
"secHigh_noFLux_fractionLow",
"secLow_noFlux_fractionHigh",
"secHigh_fLux_fractionLow",
"secLow_flux_fractionHigh"]
#l = l[2:-1]
#l = l[0:1]
threshholds = [
        [0.2,0.5]
        ]
for i in l:
    pass
    distPlot_dump(path+i+"/",extCache,fields,10)
for i in l:  
    pass
#    distPlot(path+i+"/",threshholds[0],cutoff=[0.2,10])
#distPlot("/extra/kiwitz/test/data/",threshholds[0])