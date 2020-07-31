import pandas as pd
import numpy as np
import scipy as scp
import scipy.stats as stats
import matplotlib.pyplot as plt
from copy import deepcopy
import seaborn as sns

types = ["default","Treg","effector"]
global il2

def get_A():
    A = []

    def f(dt):

        Kc = 0.01
        x = stats.gamma.cdf(dt,2,loc = 0,scale = 1)*il2/(Kc+ il2)*dt
        return x


    f0 = lambda x: 0
    A = [
        [f0, f0, f0],
        [f0, f0, f0],
        [f, f0, f0],
    ]

    A = pd.DataFrame(A,columns=types,index=types)
    return A


A = get_A()
n = 100
cells = pd.DataFrame({
    "id":range(n),
    "t":np.zeros(n),
    "type":np.random.randint(0,len(types),n,dtype=int)})

G = np.zeros((len(types),n))

def step(cells):


    result = pd.DataFrame(columns=cells.columns)
    t = cells["t"].max()
    cells_last = cells.loc[cells["t"] == t]

    dt = 0.1

    for i,c in cells_last.iterrows():


        for k,g in A[types[int(c["type"])]].items():
            c_new = deepcopy(c)

            globals()["il2"] = (0.1*len(cells.loc[cells_last["type"] == 2]))/(1 + 0.001*len(cells.loc[cells_last["type"] == 1]))


            G[types.index(k),i] += g(dt)

            if G[types.index(k),i] > np.random.uniform(0,1):

                c_new["type"] = types.index(k)
                G[types.index(k),i] = 0
                break

        c_new["t"] = t + dt
        result = result.append(c_new)

    return result

result = cells

for t in range(100):

    cells = cells.append(step(cells))


counts = cells.groupby(["t","type"]).count()
counts = counts.reset_index()

sns.lineplot(x = "t",y = "id",style = "type",data = counts)
plt.show()




