
from PySteppables import *
import CompuCell
import sys

from PlayerPython import *
from math import *
import numpy as np


    
class secretion(RunBeforeMCSSteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        RunBeforeMCSSteppableBasePy.__init__(self,_simulator,_frequency)
        
        self.valueList = []
        self.k_on = 111.7
        self.q = 10
        self.R = pow(10,2)
    def addValueToList(self,tupel):
        
        present = False
        
        for i in range(len(self.valueList)):
            if self.valueList[i]["point"] == tupel[0]:
                present = i
                break
            
        if present:
            self.valueList[present]["values"].append(tupel[1])
        else:
            self.valueList.append({"point":tupel[0],"values":[tupel[1]]})
        
    def getValues(self,p):
        for i in range(len(self.valueList)):
            if self.valueList[i]["point"] == p:
                return self.valueList[i]["values"]
        return False
    def calcNext(self,v):
        
        
        print(v)
        u = v[-1]
        return (self.q-self.k_on*u*self.R)
    def step(self,mcs):
        self.substance = CompuCell.getConcentrationField(self.simulator, "substance")
        attrSecretor=self.getFieldSecretor("substance")
        list = []
        for cell in self.cellList:
            if cell.type==1:
                pixelList = self.getCellBoundaryPixelList(cell)
                for boundaryPixelTrackerData in pixelList:
                    p = boundaryPixelTrackerData.pixel
                    self.addValueToList(((p.x,p.y,p.z),self.substance[p.x,p.y,p.z]))
        for x, y, z in self.everyPixel(1, 1, 1):
            p = (x,y,z)
            v = self.getValues(p)
            if v:
                self.substance[x,y,z] =  -1*self.calcNext(v)
            else:
                self.substance[x,y,z] = 0
    def finish(self):
        # this function may be called at the end of simulation - used very infrequently though
        return