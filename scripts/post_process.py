#!/home/kiwitz/anaconda3/envs/fenicsproject/bin/python

import sys
from PostProcess import PostProcessor

if len(sys.argv) > 1:
    path = sys.argv[1]
    options = {
            "surface_concentration":False,
            "concentration":True,
            "gradient":True
        }
    pp = PostProcessor(path)
    print(path)
    pp.dump(path, 32, options_dict=options)
    pp.prep_global_data().to_hdf(path + 'global_dataframe.h5', key="data", mode="w")
    # pp.prep_data().to_hdf(path + 'dataframe.h5', key="data", mode="w")
else:
    print("post process needs a cli parameter")