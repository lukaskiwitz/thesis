import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import N_A

p = {
    "k_on": 1e9 * 111.6 / 60 ** 2,  # 111.6 per hour per Mol
    "rho": 0.05,  # mu
    "D": (10 ** 0.5 * 0.01) ** 2,  # mu² per s
    "R_h": 4000 * N_A ** -1 * 1e9,#1.6605390671738466e-24,
    "R_l": 100 * N_A ** -1 * 1e9,
    "kd": 0,#0.1/(60*2),
    "q_h": 10 * N_A ** -1 * 1e9,
    "q_l": 1 * N_A ** -1 * 1e9,
    "fraction": 1,
    "sigma" : 0.3,
    "R_v1": 5e10 *  N_A ** -1,
    "K_R": 0.28e-9 #0.02 #0.08
}

N = 3
gammas = [0.1,10]
# gamma = p["gamma"] #10

eta = 1 / 72000  # 1/s
c_0 = 20.55e-12  # M

# k_factor = p.get_physical_parameter("k_factor", "IL-2").get_in_sim_unit()
k_factor = 2
k = k_factor * c_0  # M
R_start = 20000

c_range = np.linspace(0,500e-12,100)

for gamma in gammas:
    nenner = (gamma * c_0 ** N + k ** N) / (c_0 ** N + k ** N)
    alpha = R_start * eta / nenner  # R/t
    result = []
    for il2 in c_range:
        result.append(alpha * (gamma * il2 ** N + k ** N) / (il2 ** N + k ** N))
    plt.plot(c_range, result)
plt.show()