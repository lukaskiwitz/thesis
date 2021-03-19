import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.lines as mlines
from scipy.constants import N_A
from scipy.integrate import solve_ivp
import seaborn as sns

sns.set_style("ticks")
sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 5})
fig, ax = plt.subplots(figsize=(8, 7))


labels = ["negative feedback", "positive feedback"]
colors = ["blue", "red"]

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
a = [1/x for x in reversed(np.arange(10,110,20))]
b = [x for x in np.arange(10,110,20)]
gammas_list =[a,b]

# gammas = [1/10,10]
# gamma = p["gamma"] #10

eta = 1 / 72000  # 1/s
c_0 = 8.55e-12  # M


R_start = 1e4

# il2 = c_0  # post is nM, neeeded in M
R_il2 = 1e2 * N_A ** -1 * 1e9

R_k_max = 4e-11  # max pSTAT5 EC50 in M
R_k_min = 2e-12  # min pSTAT55 EC50 in M
R_k = 1e4 * N_A ** -1 * 1e9  # Receptors at half pSTAT5 maximum
R_N = 1.2  # hill coefficient
EC50 = (R_k_max * R_k ** R_N + R_k_min * R_il2 ** R_N) / (R_k ** R_N + R_il2 ** R_N)
R_mean = 1e4 * N_A ** -1 * 1e9


new_handles = []
c_range = np.linspace(0,200e-12,1000)
pSTAT5_list = []
for gammas in gammas_list:
    for g,gamma in enumerate(gammas):
        if gammas[0] < 1:
            palette = [x for x in reversed(sns.color_palette("Blues", len(gammas)))]
            norm = plt.Normalize(np.min(gammas), np.max(gammas))
            cmap = mpl.cm.Blues
            sm = plt.cm.ScalarMappable(cmap=cmap.reversed(), norm=norm)
            sm.set_array([])
            # plt.axes().get_legend().remove()
            # ax.figure.colorbar(sm)
        elif gammas[0] > 1:
            palette = sns.color_palette("Reds", len(gammas))
            norm = plt.Normalize(np.min(gammas), np.max(gammas))
            sm = plt.cm.ScalarMappable(cmap="Reds", norm=norm)
            sm.set_array([])
            # plt.axes().get_legend().remove()
            # ax.figure.colorbar(sm)
        else:
            raise ValueError
        result = []
        kmin = R_mean / gamma * eta
        kmax = R_mean * gamma * eta
        typical_pSTAT5 = c_0 ** 3 / (c_0 ** 3 + EC50 ** 3)
        # # k = solve((kmin * k_x ** N + kmax * typical_pSTAT5 ** N) / ((k_x ** N + typical_pSTAT5 ** N) * eta) - R_mean)[0]
        A = (typical_pSTAT5 ** N * (eta * R_mean - kmax))
        B = (kmin - eta * R_mean)
        k = (A / B) ** (1 / N)

        # k = ((c_0 ** N * (eta * R_mean - kmax)) / (kmin - eta * R_mean)) ** (1 / N)
        #
        for il2 in c_range:
            pSTAT5 = il2** 3 / (il2 ** 3 + EC50 ** 3)
            result.append((kmin * k ** N + kmax * pSTAT5 ** N) / ((k ** N + pSTAT5 ** N) * eta) * 1/(N_A ** -1 * 1e9)) #*gamma

        # plt.semilogy(c_range*1e12, result, color=palette[g])
        plt.semilogy(c_range** 3 / (c_range ** 3 + EC50 ** 3), result, color=palette[g])

        if g == 0:
            if gamma < 1:
                new_handles.append(mlines.Line2D([], [], color=palette[0], label=str(gamma)))
            else:
                new_handles.append(mlines.Line2D([], [], color=palette[-1], label=str(gamma)))
# plt.legend()
new_labels = ["neg. feedback", "pos. feedback"]
plt.legend(handles=new_handles, labels=new_labels)
ax.set_xlabel("pSTAT5 signal")
plt.ylabel("Receptors")
plt.ylim((60,2e6))
# plt.xlim((-0.02,0.41))

# ax2 = ax.twiny()
# def xticks_function(c_range, tick):
#     length = len(c_range)-1
#     print(tick*length)
#     a = [int(x) for x in (tick*length)]
#     return [round(x,0) for x in c_range[a]]
#
# new_tick_locations = np.array([0.,.25, .5, 1.])
# ax2.set_xlim(ax.get_xlim())
# ax2.set_xticks(new_tick_locations)
# ax2.set_xticklabels(xticks_function(c_range*1e12,new_tick_locations))
# ax2.set_xlabel(r"Concentration (pM)")

# plt.xlabel("Concentration (pM)")

plt.savefig("/home/brunner/Documents/Current work/04122020/hill_pSTAT5.png",bbox_inches='tight')
plt.show()

# plt.plot(c_range*1e12, c_range** 3 / (c_range ** 3 + EC50 ** 3))
# plt.ylabel("pSTAT5")
# plt.xlabel("Concentration (pM)")
# plt.savefig("/home/brunner/Documents/Current work/04122020/hill_c_pSTAT5.png",bbox_inches='tight')
# plt.show()