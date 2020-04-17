import sys

from parameters import path

from PostProcess import PostProcessor

if len(sys.argv) > 1:
    threads = int(sys.argv[1])
else:
    threads = 4


pp = PostProcessor(path)
pp.unit_length_exponent = -6

pp.write_post_process_xml(threads)

pp.make_dataframes(kde=True)
pp.global_dataframe.to_hdf(path + 'global_df.h5', key="data", mode="w")
pp.cell_dataframe.to_hdf(path+"cell_df.h5", key="df", mode="w")
pp.cell_stats.to_hdf(path+"cell_stats_df.h5", key="df", mode="w")