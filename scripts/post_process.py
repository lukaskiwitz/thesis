#!/home/kiwitz/anaconda3/envs/fenics/bin/python

import sys
from my_debug import message
from PostProcess import PostProcessor
from my_debug import message

if len(sys.argv) > 1:
    path = sys.argv[1]
    options = {
            "surface_concentration":False,
            "concentration":True,
            "gradient":True
        }
    pp = PostProcessor(path)
    message(path)
    pp.dump(path, 16, options_dict=options)
    pp.prep_global_data().to_hdf(path + 'global_dataframe.h5', key="data", mode="w")
    # pp.prep_data().to_hdf(path + 'dataframe.h5', key="data", mode="w")
else:
    message("post process needs a cli parameter")