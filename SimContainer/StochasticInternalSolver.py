#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 18:00:36 2019

@author: lukas
"""
import random
import numpy as np
import matplotlib.pyplot as plt
class BooleanInternalSolver:
    def __init__(self,dist):
        self.dist = dist
        self.random = random.Random()
        self.random.seed = 1000
    def log(self,p):
        pass
    def step(self,dT,p):
        draw = self.random.random()
        uptake = -p["flux_cytokine"]+p["q"]
        i = np.argmin(np.abs(self.dist[0]-uptake))
        plt.plot([uptake],[draw],"+")
        if draw > self.dist[1][i]:
            print("switch at: "+str(uptake))
            print(str(draw)+" and "+str(self.dist[1][i]))
        else:
            print("none at: "+str(uptake))
x = np.linspace(0,1,100)
s = 5
y = 10*x*np.exp(-s*x)
plt.plot(x,y)
sol = BooleanInternalSolver([x,y])

p = {"q":1,
     "flux_cytokine":0.5
     }
for i in range(10):
    sol.step(1,p)
