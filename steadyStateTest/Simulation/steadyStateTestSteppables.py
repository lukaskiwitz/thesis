Km = 10
vmax = 50
threshold = 1000
secrConst = 50
from PySteppables import *
import CompuCell
import sys
class steadyStateTestSteppable(SteppableBasePy):

    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):
        # any code in the start function runs before MCS=0
        pass
    def finish(self):
        # Finish Function gets called after the last MCS
        pass
        
class secretion(SecretionBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.runBeforeMCS = 1
        
    def start(self):
        self.field = CompuCell.getConcentrationField(self.simulator, "cytokine")
        uptake = 1
        for x, y, z in self.everyPixel(1, 1, 1):
            cell = self.cellField[x, y, z]
            if cell and cell.type == self.SOURCE:
                self.field[x, y, z] = -secrConst
            elif cell and cell.type == self.SINK:
                if "level" in cell.dict and cell.dict["level"] > threshold:
                    cell.type = self.NEUTRAL
                else:
                    s0 = self.field[x,y,z]
                    uptake = vmax * s0 /(Km + s0)
                    self.field[x, y, z] = uptake
                    if "level" in cell.dict:
                        cell.dict["level"] += uptake
                    else:
                        cell.dict["level"] = 0
            else:
                self.field[x, y, z] = 0.0
        
    def step(self, mcs):
        uptake = 1
        for x, y, z in self.everyPixel(1, 1, 1):
            cell = self.cellField[x, y, z]
            if cell and cell.type == self.SOURCE:
                self.field[x, y, z] = -secrConst
            elif cell and cell.type == self.SINK:
                if "level" in cell.dict and cell.dict["level"] > threshold:
                    cell.type = self.NEUTRAL
                else:
                    s0 = self.field[x,y,z]
                    uptake = vmax * s0 /(Km + s0)
                    self.field[x, y, z] = uptake
                    if "level" in cell.dict:
                        cell.dict["level"] += uptake
                    else:
                        cell.dict["level"] = 0
            else:
                self.field[x, y, z] = 0.0