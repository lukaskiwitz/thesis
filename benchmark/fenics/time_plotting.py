#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 16:59:44 2019

@author: lukas
"""
import numpy as np
import matplotlib.pyplot as plt

#cg = np.load("cg.npy")
precond = ["none","amg","hypre_amg","icc","ilu","sor"]
solver = ["lu","bicgstab","gmres","tfqmr"]
data = []
for c in precond:
    data.append(np.load("gmres_"+c+".npy"))

for d in data:
    plt.plot(d[:,0],d[:,1])
    
plt.legend(precond)
plt.xticks([16,32,64,128])
plt.xlabel("resolution")
plt.ylabel("process time in ns")
plt.savefig("preconditioner_performance.png",dpi=200)
