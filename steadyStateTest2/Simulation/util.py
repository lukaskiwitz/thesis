from PySteppables import *
import CompuCell
import sys
from PlayerPython import *
from math import *

def sourceSecretion(fields):
    return {"infg":-100,"il-2":0,"il-4":0}
def cd4Secretion(fields):
    return {"infg":fields["infg"],"il-2":fields["il-2"],"il-4":fields["il-4"]}
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
        
   