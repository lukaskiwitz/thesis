#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 12:22:58 2019

@author: kiwitz
"""
import pandas as pd
class BooleanInternalSolver:
    def __init__(self):
        self.internal = {"level":0,"threshold":2e-12,"decay":0.8}
    def log(self,p):
        df = {"flux":p["flux_cytokine"],"level":self.internal["level"],"R":p["R"]}
        return df
        
    def step(self,dT,p):
        if "flux_cytokine" not in p:
            p["flux_cytokine"] = 0
#        print(p["flux_cytokine"])
#        print("level={l}".format(l=self.internal["level"]))
#        print("-------------------------")
        uptake = p["flux_cytokine"]
        self.internal["level"] += dT*(uptake - self.internal["decay"]*self.internal["level"])
#        self.internal["level"] += dT*uptake - 
        
        
#        print("level={l}".format(l=self.internal["level"]))
        if self.internal["level"] > self.internal["threshold"] and p["R"] > 10**2:
#            print("-------------------switch---------------------")
            p["R"] = p["R"]*(1-0.5)
        elif p["R"] < 10**4:
           p["R"] = p["R"]*(1+0.5)
    def timeToStateChange(self,p):
        uptake = p["flux_cytokine"]
        
        diff = self.internal["threshold"] - self.internal["level"]
        if diff > 0 and not (p["flux_cytokine"] == 0):
            nextTimeStep= diff/uptake
            nextTimeStep *= 0.25
        else:
            nextTimeStep = False
        return nextTimeStep
        
        