#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 10:09:06 2019

@author: kiwitz
#"""
from distPlot_dump import distPlot_dump
#from distPlot import distPlot
import mpi4py.MPI as MPI

path = "/extra/kiwitz/results/"
extCache = "/extra/kiwitz/extCache/"
fields = ["il6"]

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

l = ["secHigh_dirichlet_fractionLow",
"secLow_dirichlet_fractionHigh",
"secHigh_noFLux_fractionLow",
"secLow_noFlux_fractionHigh",
"secHigh_fLux_fractionLow",
"secLow_flux_fractionHigh"]
l = l[1:2]
threshholds = [
        [0.05,6]
        ]
for i in l:
    distPlot_dump(path+i+"/",extCache,fields)
#if rank == 0:    
#    distPlot("/extra/results/"+l[0]+"/",threshholds[0])