from PySteppables import *
import CompuCell
import sys
from PlayerPython import *
from math import *

def sourceSecretion(fields):
    return {"infg":-100,"il-2":0,"il-4":0}
def cd4Secretion(sbml,cell):
    state = sbml.getSBMLState(_modelName='receptor', _cell=cell)
    #print state.items()
    RK = state["RK"]
    result = ["nextSecretion"] = {"infg":0,"il-2":RK,"il-4":0}
    
    fields = cell.dict["fields"]
    C = -1*fields["il-2"]
    self.setSBMLValue(_modelName='receptor', _valueName='C', _value=C,_cell=cell)  
    cell.dict["switch"] = state["C_i"]
    return result
def th1Secretion(fields):
    return {"infg":0,"il-2":0,"il-4":0}
def th2Secretion(fields):
    return {"infg":0,"il-2":0,"il-4":0}
    
    
def addUptake(cell,sec):
    for k,v in sec.items():
        if v > 0:
            if k in cell.dict:
                cell.dict[k] += v
            else: 
                cell.dict[k] = v
        
   