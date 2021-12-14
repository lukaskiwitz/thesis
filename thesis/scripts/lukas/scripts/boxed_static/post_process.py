import sys

from parameters import path
from thesis.main.PostProcess import PostProcessor

"""number of threads can be passed as first cli argument"""
if len(sys.argv) > 1:
    threads = int(sys.argv[1])
else:
    threads = 4

"""
setting filepath to look for sim results. This is setup so that it works on the itb computers.
"""

pp = PostProcessor(path)
pp.unit_length_exponent = -6

"""appends a custom calculation.
default computations are defined in PostProcess.py"""
pp.image_settings = {
    "cell_colors": "Dark2",
    "cell_color_key": "type_name",
    "round_legend_labels": 4,
    "legend_title": "",
    "dpi": 350,
    # "colorbar_range":(0.6,1.3)
}


"""carries out the operations in pp.computations in parallel and stores the result in xml file"""
pp.run_post_process(threads, make_images=True, kde=False)

pp.global_dataframe.to_hdf(path + 'global_df.h5', key="data", mode="w")
pp.cell_dataframe.to_hdf(path + "cell_df.h5", key="df", mode="w")
