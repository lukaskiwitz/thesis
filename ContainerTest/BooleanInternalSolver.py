#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 12:22:58 2019

@author: kiwitz
"""
import pandas as pd
class BooleanInternalSolver:
    def __init__(self):
        self.internal = {"level":0}
    def log(self,p):
        df = {"flux":p["flux_cytokine"],"level":self.internal["level"],"R":p["R"]}
        return df
        
    def step(self,dT,p):
        if "flux_cytokine" not in p:
            p["flux_cytokine"] = 0
#        print(p["flux_cytokine"])
#        print("level={l}".format(l=self.internal["level"]))
#        print("-------------------------")
        uptake = -p["flux_cytokine"]+p["q"]
        self.internal["level"] += dT*(uptake - p["decay"]*self.internal["level"])
        
        
#        print("level={l}".format(l=self.internal["level"]))
        if self.internal["level"] > p["threshold"] and p["R"] > p["low"]:
            p["R"] = p["R"]*(1-dT*0.1)
        elif self.internal["level"] < p["threshold"] and p["R"] < p["high"]:
           p["R"] = p["R"]*(1+dT*0.1)
    def timeToStateChange(self,p):
        uptake = p["flux_cytokine"]
        
        diff = p["threshold"] - self.internal["level"]
        if diff > 0 and not (p["flux_cytokine"] == 0):
            nextTimeStep= diff/uptake
            nextTimeStep *= 0.25
        else:
            nextTimeStep = False
        return nextTimeStep
        
        