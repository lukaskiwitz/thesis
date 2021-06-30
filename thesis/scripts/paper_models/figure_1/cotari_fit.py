import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import  curve_fit




# R_M = 1e3
# ec50_min = 0.0015 *1e3
# ec50_max = 0.8 * 1e3
# N = 0.8
# n_il2 = 4

def f(R , ec50_min, ec50_max, N, R_M):

    return  (ec50_max - ec50_min) * (1 - R ** N / (R ** N + R_M ** N)) + ec50_min

x = np.logspace(1,8,100)

data = pd.read_csv("cotari_ec50.csv")

p, cov = curve_fit(f,data["x"],data["y"],
                   [0.0015 *1e3, 0.8 * 1e3, 0.8, 1e3],
                   bounds = (0,np.inf))


plt.plot(x,f(x, p[0],p[1],p[2],p[3]))
plt.scatter(data["x"],data["y"])

p = [np.round(i,2) for i in p]

plt.figtext(0.5,0.7,"EC50_min = {min} pM"
                    "\nEC50_max = {max} pM"
                    "\nn = {n}"
                    "\nR_M = {rm:e}".format(min = p[0],max = p[1], n = p[2], rm = p[3]))

plt.xlim([1e1,1e8])
plt.ylim([1,1e3])
plt.loglog()
plt.savefig("cotari_ec50_fit.pdf")
plt.show()