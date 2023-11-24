import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks

rc_ticks['figure.figsize'] = (1.7, 0.7)
sns.set_theme(context = "talk", style = "ticks", rc = rc_ticks)
fig,ax = plt.subplots()


def satu_uptake(k_endo, R, KD, u):
    return k_endo * R * u / (KD + u)

def lin_uptake(k_on, R, u):
    return k_on * R * u

N = 3
k_on = 3.1e7 # 1/(M*s)
k_endo = 0.00046
KD = 7.437e-12

c_range = np.linspace(0,50e-12,1000) #M

plt.plot(c_range * 1e12, satu_uptake(k_endo, 1.5e3, KD, c_range), label = "saturated", color="#aa0000ff", alpha=1)
plt.plot(c_range * 1e12, lin_uptake(k_on, 1.5e3, c_range), label = "linear", color="blue", alpha=0.8)

plt.ylabel("uptake rate" "\n" r"(mol./s)")
plt.xlabel("concentration (pM)")
plt.xlim((0,50))
plt.ylim((0,1))
plt.xticks([0,25,50])
plt.yticks([0,0.5,1])
saving_string = r"/home/brunner/Documents/Current work/2023_08_04/uptake.pdf"

plt.savefig(saving_string, bbox_inches="tight", transparent=True)
plt.tight_layout()
plt.show()