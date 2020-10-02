import numpy as np
from scipy.constants import N_A
from sympy import symbols, solve
from scipy.integrate import solve_ivp

R = 1e4


# p.get_physical_parameter("R", "IL-2").set_in_post_unit(dt*2*R)
def func(t, y):
    R = y[0]
    dR_dt = 5-0.1*R
    return [dR_dt]


y0 = [R]

T = np.linspace(0,10000,500)
result = [R]
for t_idx, t in enumerate(T[1:]):
    t_space = np.linspace( T[t_idx], t, 2)
    # print(t_space)
    sol = solve_ivp(func, t_span=[T[t_idx], t], t_eval=t_space, y0=y0, method="LSODA", max_step=1)
    result.append(sol.y[0][-1])
    y0 = [sol.y[0][-1]]
print("done")
# print(result)
import matplotlib.pyplot as plt
#
plt.plot(T, result)
plt.show()