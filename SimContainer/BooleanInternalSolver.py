#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 12:22:58 2019

@author: kiwitz
"""

class BooleanInternalSolver:
    internal = {"level":0,"threshold":0.1}
    
    def step(self,dT,p):
        if "flux_cytokine" not in p:
            p["flux_cytokine"] = 0
#        print("-----------------------_")
#        print(p["flux_cytokine"])
#        print("-------------------------")
        
        self.internal["level"] += dT*p["flux_cytokine"]
        print("level={l}".format(l=self.internal["level"]))
        if self.internal["level"] > self.internal["threshold"]:
            print("-------------------switch---------------------")
            p["R"] = 10**2
        else:
            p["R"] = 10**4
    def timeToStateChange(self,p):
        diff = self.internal["threshold"] - self.internal["level"]
        if diff > 0 and not (p["flux_cytokine"] == 0):
            nextTimeStep= diff/p["flux_cytokine"]
            nextTimeStep *= 0.25
        else:
            nextTimeStep = False
        return nextTimeStep
        
        