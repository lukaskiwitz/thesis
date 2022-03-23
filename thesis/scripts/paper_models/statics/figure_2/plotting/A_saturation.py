import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from plotting_rc import rc_ticks

sns.set_theme(context = "talk", style = "ticks", rc = rc_ticks)
fig,ax = plt.subplots()


def satu_uptake(k_endo, R, KD, u):
    return k_endo * R * u / (KD + u)

def lin_uptake(k_on, R, u):
    return k_on * R * u

N = 3
k_on = 3.1e7 # 1/(M*s)
k_endo = 0.0011
KD = 7.437e-12

c_range = np.linspace(0,100e-12,1000) #M

plt.plot(c_range * 1e12, lin_uptake(k_on,1.5e3, c_range), label = "linear")
plt.plot(c_range * 1e12, satu_uptake(k_endo,1.5e3, KD, c_range), label = "saturated")

plt.ylabel("uptake rate (mol./s)")
plt.xlabel("concentration (pM)")
plt.legend()
saving_string = r"/home/brunner/Documents/Current work/2022_03_04/uptake.pdf"

plt.savefig(saving_string, bbox_inches="tight")
plt.tight_layout()
plt.show()