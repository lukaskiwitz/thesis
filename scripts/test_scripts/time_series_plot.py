import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.constants import N_A
from copy import deepcopy


import StateManager as ST

st = ST.StateManager("/extra/kiwitz/test_results/")
st.loadXML()
result: pd.DataFrame = st.get_cell_ts_data_frame()


result["center"] = result["center"].apply(lambda x: (x[0],x[1],x[2]))
result = result.groupby("id").apply(lambda x: x.ffill().bfill()).drop_duplicates()
r = deepcopy(result)
r["R_il2"] = result["R_il2"].div(N_A**-1*10e9)# to mol/cell
r["R_il6"] = result["R_il6"].div(N_A**-1*10e9)# to mol/cell
r["surf_c_il2"] = result["surf_c_il2"].mul(10e9) # to nM
r["surf_c_il6"] = result["surf_c_il6"].mul(10e9) # to nM
r["surf_c_infg"] = result["surf_c_infg"].mul(10e9) # to nM
r["t"] = r["t"].apply(lambda x: float(x))

r.insert(0,"x",r["center"].apply(lambda x: x[0]))


sns.set_context("paper",font_scale=1,rc={
            "lines.markersize":0,
            "lines.linewidth":1
            }
            )

g = r.groupby(["t","type_name"],as_index=False)
g = g.count()
g = g.drop(g.columns.drop(["t","type_name","n","x"]),axis=1)

g2 = r.groupby(["t","x","type_name"],as_index=False)
g2 = g2.mean()
g2 = g2.drop(g2.columns.drop(["x",
                              "t",
                              "type_name",
                              "surf_c_il2",
                              "surf_c_il6",
                              "surf_c_infg",
                              "surf_g_il2",
                              "surf_g_il6",
                              "surf_g_infg"]),axis=1)
plt.figure(1)
plt.subplot(2,2,1)
sns.lineplot(x="t",y="n",data=g, hue="type_name",ci="sd",sort=False,style="type_name",legend=False)

plt.ylim((0,g["n"].max()*(1.2)))

plt.subplot(2,2,2)
sns.lineplot(x="t",y="surf_c_il2",data=g2, hue="type_name",style="type_name",ci="sd",sort=False,legend=False)

plt.subplot(2,2,3)
sns.lineplot(x="t",y="surf_c_il6",data=g2, hue="type_name",style="type_name",ci="sd",sort=False,legend=False)

plt.subplot(2,2,4)
sns.lineplot(x="t",y="surf_c_infg",data=g2, hue="type_name",style="type_name",ci="sd",sort=False,legend="brief")

plt.show()
#
# plt.figure(2)
# plt.subplot(2,2,1)
# sns.lineplot(x="t",y="n",data=g, hue="type_name",ci="sd",sort=False,style="type_name",legend=False)
#
# plt.ylim((0,g["n"].max()*(1.2)))
#
# plt.subplot(2,2,2)
# sns.lineplot(x="t",y="surf_g_il2",data=g2, hue="type_name",style="type_name",ci="sd",sort=False,legend=False)
#
# plt.subplot(2,2,3)
# sns.lineplot(x="t",y="surf_g_il6",data=g2, hue="type_name",style="type_name",ci="sd",sort=False,legend=False)
#
# plt.subplot(2,2,4)
# sns.lineplot(x="t",y="surf_g_infg",data=g2, hue="type_name",style="type_name",ci="sd",sort=False,legend=False)
#
# plt.show()

# plt.savefig("/home/kiwitz/ts_result.pdf",dpi=600)

# def f( x):
#     # C = p["k_on"] * x[0] * p["surface_concentration_il2"]
#     C = x
#
#
#     K = 1000
#     v1 = 3000
#     v0 = 150
#     k_deg = 0.5
#     n = 12
#     # print(R / (N_A ** -1 * 10e9))
#     return (v1 * (C ** n) / (K ** n + C ** n))
#
# x = np.linspace(0,4000,1000)
# y = [f(i) for i in x]
# plt.plot(x,y)
# plt.show()
