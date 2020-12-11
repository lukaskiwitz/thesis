try:
    import fenics as fcs
except RuntimeError:
    import os
    os.environ['PATH'] = '/home/brunner/anaconda3/envs/Lukas2/bin:/home/brunner/.local/bin:/home/brunner/anaconda3/condabin:/usr/local/bin:/usr/bin:/bin:/usr/local/games:/usr/games:/opt/puppetlabs/bin'
    import fenics as fcs

import getpass
import random
import sys
import os

sys.path.append("/home/brunner/thesis/thesis/main/")
sys.path.append("/home/brunner/thesis/thesis/scenarios/")

import numpy as np
from scipy.constants import N_A

import matplotlib.pyplot as plt
import mpi4py.MPI as MPI
from scipy.constants import N_A

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

N = 3

gammas = [1e-2, 1e2]
# gamma = p["gamma"] #10

eta = 1 / 72000  # 1/s
c_0 = 8.6e-12 # M
myRange = np.linspace(0.01e-18, 50e-12, 1000)

# k_factor = p.get_physical_parameter("k_factor", "IL-2").get_in_sim_unit()
# k_factor = 1
# k = k_factor * c_0  # M
R_start = 1e4

# kmax = 1e6*eta

#################################### hill_full #########################################

#
# results = [[] for x in range(len(gammas))]
# # results = [[] for x in range(2)]
# for g,gamma in enumerate(gammas):
#     for n,N in enumerate([3]):
#         # if gamma > 1:
#         #     kmin = 1e2 * eta
#         #     kmax = gamma * kmin
#         # else:
#         #     kmin = 1e6 * eta
#         #     kmax = gamma * kmin
#         kmin = 1e4/gamma * eta
#         kmax = 1e4*gamma * eta
#         from sympy import symbols, solve
#         k_x = symbols("k_x")
#         try:
#             k = solve((kmin * k_x ** N + kmax * 8.6e-12 ** N) / ((k_x ** N + 8.6e-12 ** N) * eta) - 1e4)[0]
#         except:
#             k = 1
#
#         for il2 in myRange:
#             # kmin = kmax/gamma
#             # nenner = (gamma * c_0 ** N + k ** N) / (c_0 ** N + k ** N)
#             # alpha = a * R_start * eta / nenner  # R/t
#             results[g].append((kmin * k ** N + kmax * il2 ** N) / (( k ** N + il2 ** N) * eta))
#             # results[g].append(3600 * (alpha * (gamma * il2 ** N + k ** N) / (il2 ** N + k ** N) - eta*R_start))
#
# import seaborn as sns
# sns.set(font_scale=2, rc={"lines.linewidth":2})
# sns.set_style("white")
# labels = ["negative feedback", "positive feedback"]
# for i in range(len(results)):
#     plt.semilogy(myRange*1e12, results[i], label=labels[i])
# plt.axvline(8.6, color="black")
# plt.xlabel("Concentration in pM")
# plt.ylabel("Receptors")
# # plt.title("g = " + str(gammas[0]))
# plt.legend()
# plt.savefig("/home/brunner/Documents/Current work/25092020/hill_full.png",bbox_inches='tight')
# plt.show()
#
# myRange = np.linspace(0, 1, 1000)
# results2 = []
# kon = 100
# R = 1e4
# K_c = 0.0036
# for il2 in myRange:
#     results2.append(kon*R*K_c*il2/(K_c + il2) * 1/3600)
# plt.plot(myRange*1e12, results2)
# plt.show()
# gamma = 1e-2
# from sympy import symbols, solve
# k_x = symbols("k_x")
# if gamma > 1:
#     kmin = 1e4 * eta
#     kmax = gammas[0] * kmin
# else:
#     kmin = 1e4 * eta
#     kmax = gammas[0] * kmin
# k = solve((kmin * k_x ** N + kmax * 8.6e-12 ** N) / (( k_x ** N + 8.6e-12 ** N) * eta) - 1e4)[0]
# k = 10.85281383342742e-12
#

 #################################### hill_halfed #########################################

# print(k)
# gammas = [1e-2,1e2]
# results2 = [[] for x in range(len(gammas))]
# for g,gamma in enumerate(gammas):
#     for n,N in enumerate([3]):
#         for il2 in myRange:
#             if gamma > 1:
#                 kmin = 1e4 * eta
#                 kmax = gamma * kmin
#                 k = 8*8e-12
#             else:
#                 kmin = 1e4*eta
#                 kmax = gamma * kmin
#                 k = 8e-12
#             results2[g].append((kmin * k ** N + kmax * il2 ** N) / (( k ** N + il2 ** N) * eta))
# import seaborn as sns
# sns.set(font_scale=2, rc={"lines.linewidth":2})
# sns.set_style("white")
# labels = ["negative feedback", "positive feedback"]
# for i in range(len(results2)):
#     plt.semilogy(myRange*1e12, results2[i], label=labels[i])
# plt.axvline(8, color="black")
# plt.xlabel("Concentration in pM")
# plt.ylabel("Receptors")
# # plt.title("g = " + str(gammas[0]))
# plt.legend()
# plt.savefig("/home/brunner/Documents/Current work/25092020/hill_halfed.png",bbox_inches='tight')
# plt.show()

#################################### hill_pSTAT5 #########################################

myRange = np.linspace(0,0.19,500)
from sympy import symbols, solve
k_x = symbols("k_x")
fig = plt.figure(figsize=(10,8))
results = [[] for x in range(len(gammas))]
# results = [[] for x in range(2)]
for g,gamma in enumerate(gammas):
    for n,N in enumerate([3]):
        R_k_max = 4e-11 # max pSTAT5 EC50 in M
        R_k_min = 2e-12 # min pSTAT55 EC50 in M
        R_k = 1e4*N_A ** -1 * 1e9 # Receptors at half pSTAT5 maximum
        R_N = 1.2 # hill coefficient

        R_il2 = 1e2*N_A ** -1 * 1e9
        EC50 = (R_k_max * R_k ** R_N + R_k_min * R_il2 ** R_N) / (R_k ** R_N + R_il2 ** R_N)

        R_mean = 1e4 * N_A ** -1 * 1e9
        kmin = R_mean/gamma * eta
        kmax = R_mean*gamma * eta

        typical_pSTAT5 = c_0**3/(c_0**3 + EC50**3)
        k = ((typical_pSTAT5 ** N * (eta * R_mean - kmax)) / (kmin - eta * R_mean)) ** (1 / N)
        # k = solve((kmin * k_x ** N + kmax * typical_pSTAT5 ** N) / (( k_x ** N + typical_pSTAT5 ** N) * eta) - R_mean)[0]
        print(k)

        for pSTAT5 in myRange:
            # kmin = kmax/gamma
            # nenner = (gamma * c_0 ** N + k ** N) / (c_0 ** N + k ** N)
            # alpha = a * R_start * eta / nenner  # R/t
            results[g].append((kmin * k ** N + kmax * pSTAT5 ** N) / (( k ** N + pSTAT5 ** N) * eta))
            # results[g].append(3600 * (alpha * (gamma * il2 ** N + k ** N) / (il2 ** N + k ** N) - eta*R_start))

import seaborn as sns
sns.set(font_scale=2.5, rc={"lines.linewidth":3})
sns.set_style("white")
labels = ["negative feedback", "positive feedback"]
colors = ["blue", "red"]
for i in range(len(results)):
    plt.semilogy(myRange, np.array(results[i])/(N_A ** -1 * 1e9), label=labels[i], color=colors[i])
# plt.axvline(8.6, color="black")
plt.xlabel("pSTAT5 signal")
plt.ylabel("Receptors")
# plt.title("g = " + str(gammas[0]))
plt.legend()
plt.savefig("/home/brunner/Documents/Current work/30102020/hill_pSTAT5.png",bbox_inches='tight')
plt.show()