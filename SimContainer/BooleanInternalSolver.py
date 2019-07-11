#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 12:22:58 2019

@author: kiwitz
"""

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
        uptake = round(p["flux_cytokine"],14)
        self.internal["level"] += dT*(uptake - p["decay"]*self.internal["level"])
        
        
#        print("level={l}".format(l=self.internal["level"]))
        if self.internal["level"] > p["threshold"] and p["R"] > p["low"]:
            p["R"] = p["R"]*(1-dT*p["receptor_decay"])
        elif self.internal["level"] < p["threshold"] and p["R"] < p["high"]:
           pass#p["R"] = p["R"]*(1+dT*p["receptor_decay"])
    def timeToStateChange(self,p):
        uptake = p["flux_cytokine"]
        
        diff = p["threshold"] - self.internal["level"]
        if diff > 0 and not (p["flux_cytokine"] == 0):
            nextTimeStep= diff/uptake
            nextTimeStep *= 0.25
        else:
            nextTimeStep = False
        return nextTimeStep
        
        