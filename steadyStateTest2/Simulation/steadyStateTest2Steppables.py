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
            if self.il2[x,y,z] < 0:
                print "negative"
            #fields = {"infg":self.infg[x, y, z],"il-2":self.il2[x, y, z],"il-4":self.il4[x, y, z]}
            if cell:
                cell.dict["fields"] = {"infg":self.infg[x, y, z],"il-2":self.il2[x, y, z],"il-4":self.il4[x, y, z]}
                
                if cell.type == self.CD4:
                    nextSecretion = cd4Secretion(self.getSteppableByClassName("SBML"),cell)
                elif cell.type == self.SOURCE:
                    nextSecretion= {"infg":0,"il-2":-10000,"il-4":0}
                else:
                    nextSecretion = {"infg":0,"il-2":0,"il-4":0}
                    
                self.infg[x, y, z] = nextSecretion["infg"]
                self.il2[x, y, z] = nextSecretion["il-2"]
                self.il4[x, y, z] = nextSecretion["il-4"]
            else:
                self.infg[x, y, z] = 0.0
                self.il2[x, y, z] = 0.0
                self.il4[x, y, z] = 0.0
class cellSwitch(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.runBeforeMCS = 1
    def start(self):
        print "cellSwitch: This function is called once before simulation"
        
    def step(self,mcs):
        for cell in self.cellList:
            if cell.type == self.CD4:
                if "switch" in cell.dict and cell.dict["switch"] > 10:
                    print "switch-------------------------------------------------------"
                    cell.type = self.TH1
            elif cell.type == self.TH1:
                pass
            elif cell.type == self.TH2: 
                pass

class SBML(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.runBeforeMCS = 1
    def start(self):
        modelFile = './sbml.xml'  
        
        initialConditions = {}
        initialConditions['C'] = 0
        initialConditions['RK'] = 1
        
        stepSize = 0.5
        self.addSBMLToCellTypes(_modelFile=modelFile, _modelName='receptor', _types=[self.CD4],
                                _stepSize=stepSize, _initialConditions=initialConditions)
    def step(self,mcs):
        self.timestepSBML()
#         for cell in self.cellList:
#             if cell.type == self.CD4:
#                 pass
#             elif cell.type == self.SOURCE:
#                 cell.dict["nextSecretion"] = {"infg":0,"il-2":-10000,"il-4":0}
#             else:
#                 cell.dict["nextSecretion"] = {"infg":0,"il-2":0,"il-4":0}
            
            
            
