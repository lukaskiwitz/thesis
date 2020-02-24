import getpass
from PostProcess import PostProcessor
import sys

if len(sys.argv) > 1:
    threads = int(sys.argv[1])
else:
    threads = 4

user = getpass.getuser()
path = "/extra/{u}/scan_example_small/".format(u=user)

pp = PostProcessor(path)
# pp.write_post_process_xml(threads)

pp.make_dataframes(kde=True)
pp.global_dataframe.to_hdf(path + 'global_df.h5', key="data", mode="w")
pp.cell_dataframe.to_hdf(path+"cell_df.h5", key="df", mode="w")
pp.cell_stats.to_hdf(path+"cell_stats_df.h5", key="df", mode="w")