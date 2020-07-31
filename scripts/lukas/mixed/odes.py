import fenics as fcs
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import N_A, pi
from scipy.integrate import odeint

t = 0

x_list = []
u1_l = []
u2_l = []
A = (5 * 4 * pi ** 2)
l = np.sqrt(10)
T = np.linspace(0, 10, 1000)
# assigner = fcs.FunctionAssigner(M,M)

k_on = (111.6) / (60 ** 2) * 1e9 * 1e15 * 1 / N_A
k_iR = 0
k_off = 1
k_rec = 1
k_iC = 10
k_deg = 0
v_trans = 15000 / (60 ** 2) / A * l

D = fcs.Constant(10)  # mu² per s
q = 10 / A * l

r1 = 1000 / A * l
r2 = 1000 / A * l
c1 = 0
c2 = 0

e1 = 0
e2 = 0


# u1 = 1

def get_ap(f):
    def ode(x, t, f):
        r1, c1, e1 = x
        u1 = np.sin(f * t) + 1
        x_n = [
            v_trans - (k_on * u1 + k_iR) * r1 + k_off * c1 + k_rec * e1,
            k_on * u1 * r1 - (k_off + k_iC) * c1,
            k_iC * c1 - (k_rec + k_deg) * e1
        ]
        return x_n

    result = odeint(ode, [r1, c1, e1], T, args=(f,))
    e = result[3]
    e = e[int(len(e) * 0.5):len(e)]
    # plt.plot(result)
    # plt.show()
    return np.max(e) - np.min(e)


ap_list = []
for f in np.linspace(1, 1000, 100):
    ap_list.append(get_ap(f))
plt.plot(ap_list)
plt.show()
