#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 12:22:13 2019

@author: kiwitz
"""

import MySubDomain as SD

class Entity:
    pass
class Cell(Entity):
    def __init__(self,center,radius,bcDict):
        self.center = center
        self.radius = radius
        self.bcDict = bcDict
        super().__init__()
    def getSubDomain(self):
            return SD.CellSubDomain(self.center,self.radius)
    def getBcDict(self):
        return self.bcDict
    
class DomainEntity(Entity):
    def __init__(self):
        super().__init__()
        
class DomainSphere(DomainEntity):
    def __init__(self,center,radius,bcDict):
        self.center = center
        self.radius = radius
        self.bcDict = bcDict
        super().__init__()
    def getSubDomain(self):
            return SD.OuterSphere(self.center,self.radius)
    def getBcDict(self):
        return self.bcDict