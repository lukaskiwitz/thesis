
import sys
from os import environ
from os import getcwd
import string

sys.path.append(environ["PYTHON_MODULE_PATH"])


import CompuCellSetup


sim,simthread = CompuCellSetup.getCoreSimulationObjects()
        
# add extra attributes here
        
CompuCellSetup.initializeSimulationObjects(sim,simthread)
# Definitions of additional Python-managed fields go here
        
#Add Python steppables here
steppableRegistry=CompuCellSetup.getSteppableRegistry()
        
from steadyState_benchmarkSteppables import secretion
instanceOfsecretion=secretion(_simulator=sim,_frequency=1)
steppableRegistry.registerSteppable(instanceOfsecretion)


CompuCellSetup.mainLoop(sim,simthread,steppableRegistry)
        
        