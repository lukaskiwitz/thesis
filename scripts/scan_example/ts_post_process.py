from PostProcess import PostProcessor
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

path = "/extra/kiwitz/scan_example_large/"
#
# pp = PostProcessor(path)
# pp.write_post_process_xml(32)
# #
# pp.get_global_dataframe().to_hdf(path + 'global_df.h5', key="data", mode="w")
# #
# pp.get_cell_dataframe(kde=True).to_hdf(path+"cell_df.h5", key="df", mode="w")
# #
#
#


f = "il6"
global_df = pd.read_hdf(path+"global_df.h5", mode="r")

cell_df = pd.read_hdf(path+"cell_df.h5",mode="r")

plt.subplot(2,2,1)
# plt.figure()
sns.lineplot(x="t",y="surf_c_{f}".format(f=f),data=cell_df)
# plt.show()
#

plt.subplot(2,2,2)
sns.lineplot(x="timeIndex", y="concentration", data=global_df,hue="field_name")

plt.subplot(2,2,3)
sns.lineplot(x="timeIndex", y="gradient", data=global_df,hue="field_name")
plt.show()




# plt.show()
# plt.subplot(2,2,3)
# plt.figure()
# sns.scatterplot(data=cell_df.loc[(cell_df["t"] == 1)],x="id",y="surf_c_{f}".format(f=f))

# plt.subplot(2,2,4)
# sns.distplot(cell_df.loc[cell_df["t"] == 1]["surf_c_{f}".format(f=f)],kde=False)
# plt.show()