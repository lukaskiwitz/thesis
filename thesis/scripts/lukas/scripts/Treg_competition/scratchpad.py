import matplotlib.pyplot as plt
import numpy as np

from thesis.cellBehaviourUtilities.cell_solver import kineticSolver

c_naive = (0.9058823529411765, 0.1607843137254902, 0.5411764705882353)
c_treg = (0.9019607843137255, 0.6705882352941176, 0.00784313725490196)

s = kineticSolver()
x = np.linspace(0, 8e3, 1000)
y = s.EC50_calculation(E_max=125e-12, E_min=0, k=860, N=1, R=x)
plt.plot(x, y, color=c_naive)
ylim = plt.gca().get_ylim()
plt.vlines(100, *ylim, color=c_naive)

y = s.EC50_calculation(E_max=125e-12, E_min=0, k=0.25 * 860, N=1, R=x)
plt.plot(x, y, color=c_treg)
plt.vlines(1000, *ylim, color=c_treg)

plt.legend(["naive", "treg"])
plt.show()
