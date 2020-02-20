import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from my_debug import message
large_df = pd.read_hdf("/home/kiwitz/master_global.h5",key="df",mode="r")
# large_df = large_df.loc[large_df["timeIndex"] == 10]

# plt.subplot(2,2,1)
# sns.lineplot(
#     x="scanIndex",y="concentration",hue="scan_name",data=large_df.loc[large_df["field_name"] == "il2"]
#     ,legend=False)
# plt.subplot(2,2,2)
# sns.lineplot(
#     x="scanIndex",y="concentration",hue="scan_name",data=large_df.loc[large_df["field_name"] == "il6"]
#     ,legend=False)
# plt.subplot(2,2,3)
# sns.lineplot(
#     x="scanIndex",y="concentration",hue="scan_name",data=large_df.loc[large_df["field_name"] == "infg"]
#     ,legend=False)
# plt.show()
#
# plt.figure()
# plt.subplot(2,2,1)
# sns.lineplot(
#     x="scanIndex",y="gradient",hue="scan_name",data=large_df.loc[large_df["field_name"] == "il2"]
#     ,legend=False)
# plt.subplot(2,2,2)
# sns.lineplot(
#     x="scanIndex",y="gradient",hue="scan_name",data=large_df.loc[large_df["field_name"] == "il6"]
#     ,legend=False)
# plt.subplot(2,2,3)
# sns.lineplot(
#     x="scanIndex",y="gradient",hue="scan_name",data=large_df.loc[large_df["field_name"] == "infg"]
#     ,legend=False)
# plt.show()

# plt.savefig(IMGPATH + subfolder+".pdf", dpi=600)


plt.figure()
sns.lineplot(
    x="timeIndex",y="gradient",hue="scan_name",data=large_df.loc[
        (large_df["field_name"] == "infg")&
        (large_df["scanIndex"] == 2)

        ]
    ,legend=False)
plt.show()