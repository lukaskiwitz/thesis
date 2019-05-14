from PySteppables import *
import CompuCell
import sys
from PlayerPython import *
from math import *

from util import *

class secretion(SecretionBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.runBeforeMCS = 1
    def start(self):
        pass
    def step(self,mcs):
        self.infg = CompuCell.getConcentrationField(self.simulator, "INFg")
        self.il2 = CompuCell.getConcentrationField(self.simulator, "il-2")
        self.il4 = CompuCell.getConcentrationField(self.simulator, "il-4")
        
        for x, y, z in self.everyPixel(1, 1, 1):
            cell = self.cellField[x, y, z]
            fields = {"infg":self.infg[x, y, z],"il-2":self.il2[x, y, z],"il-4":self.il4[x, y, z]}
            if cell:
                if cell.type == self.CD4:
                    sec = cd4Secretion(fields)
                    addUptake(cell,sec)
                    self.infg[x,y,z] = sec["infg"]
                    self.il2[x,y,z] = sec["il-2"]
                    self.il4[x,y,z] = sec["il-4"]
                elif cell.type == self.TH1:
                    sec = th1Secretion(fields)
                    addUptake(cell,sec)
                    self.infg[x,y,z] = sec["infg"]
                    self.il2[x,y,z] = sec["il-2"]
                    self.il4[x,y,z] = sec["il-4"]
                elif cell.type == self.TH2:
                    sec = th2Secretion(fields)
                    addUptake(cell,sec)
                    self.infg[x,y,z] = sec["infg"]
                    self.il2[x,y,z] = sec["il-2"]
                    self.il4[x,y,z] = sec["il-4"]
                elif cell.type == self.SOURCE:
                    sec = sourceSecretion(fields)
                    self.infg[x,y,z] = sec["infg"]
                    self.il2[x,y,z] = sec["il-2"]
                    self.il4[x,y,z] = sec["il-4"]
            else:
                self.infg[x, y, z] = 0.0
                self.il2[x, y, z] = 0.0
                self.il4[x, y, z] = 0.0
class cellSwitch(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        
    def start(self):
        print "cellSwitch: This function is called once before simulation"
        
    def step(self,mcs):
        for cell in self.cellList:
            if cell.type == self.CD4:
                if "infg" in cell.dict and cell.dict["infg"] > 10:
                    cell.type = self.TH1
            elif cell.type == self.TH1:
                pass
            elif cell.type == self.TH2: 
                pass
