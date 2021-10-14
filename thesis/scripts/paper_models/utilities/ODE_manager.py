import numpy as np
from scipy.constants import N_A
from joblib import Parallel, delayed
from copy import deepcopy

class ODE_manager:


    def __init__(self,p):
        self.p = deepcopy(p)
        self.lognorm_mean = None
        self.lognorm_variance = None
        self.draws = None


    def get_volume_in_m(self):
        # return 2.2641509999999996e-12
        return ((self.p["N_cells"] ** (1/3) * self.p["d"] + self.p["m"])**3 - (self.p["N_cells"] * 4/3 * np.pi * self.p["rho"]**3)) * 1e-18

    def molar_to_molecules(self, x):
        V = self.get_volume_in_m()
        return x * (N_A * V * 1e3) # 1e3 bcs M = mol/l = 1000 mol/m^3

    def molecules_to_molar(self, x):
        V = self.get_volume_in_m()
        return x / (N_A * V * 1e3)

    def ng_per_milliliter_to_molar(self, value, molar_weight):
        value = 1e-9/1e-3 * value # to g/L
        return value/molar_weight

    def uptake(self, R, IL2, type = "linear"):
        if type == "linear":
            return self.p["k_on"] / self.molar_to_molecules(1) * R * IL2
        elif type == "I_saturation":
            return self.p["k_endo"] * R * IL2 / (self.p["K_D"] + IL2)


    def exp_time_range(self, dt, max_T, length):
        max_T = max_T*dt
        myRange = np.arange(0, length)

        def exp_func(x, a, b, c):
            return np.round(a * np.exp(b * x) + c, 4)

        a = 2 * dt
        c = -a
        b = np.log((max_T - c) / a) / (length - 1)
        return(exp_func(myRange, a, b, c))

    def set_R_feedback(self):
        if self.p["gamma"] < 1:
            self.p["k_min"] = 1e5 * self.p["nu"]
            self.p["k_max"] = 1e5 * self.p["nu"] * self.p["gamma"]
        elif self.p["gamma"] > 1:
            self.p["k_min"] = 1e2 * self.p["nu"]
            self.p["k_max"] = 1e2 * self.p["nu"] * self.p["gamma"]
        elif self.p["gamma"] == 1:
            self.p["k_min"] = self.p["k_max"] = self.p["R_mean"] * self.p["nu"]
            self.p["K_m"] = None
        return self.p

    def set_R_feedback_old(self):
        if self.p["gamma"] < 1:
            self.p["k_min"] = 1e5 * self.p["eta"]
            self.p["k_max"] = 1e4 * self.p["gamma"] * self.p["eta"]
        elif self.p["gamma"] > 1:
            self.p["k_min"] = 1e2 * self.p["eta"]
            self.p["k_max"] = 5e3 * self.p["gamma"] * self.p["eta"]
        elif self.p["gamma"] == 1:
            self.p["k_min"] = self.p["k_max"] = self.p["R_mean"] * self.p["eta"]

        if self.p["gamma"] != 1:
            # k = solve((self.p["k_min"]* k_x ** N + self.p["k_max"]* c_0 ** N) / ((k_x ** N + c_0 ** N) * self.p["eta"]) - R_mean)[0]
            self.p["K_m"] = ((self.p["c_0"] ** self.p["N_R"] * (self.p["eta"] * self.p["R_mean"] - self.p["k_max"])) / (
                    self.p["k_min"] - self.p["eta"] * self.p["R_mean"])) ** (1 / self.p["N_R"])
        else:  # if self.p["gamma"]~ 1
            self.p["K_m"] = None
        return self.p

    def set_R_lognorm(self, mean, variance, length = 1):
        tmp_sigma = np.sqrt(np.log(variance ** 2 / mean ** 2 + 1))
        mean = np.log(mean) - 1 / 2 * tmp_sigma ** 2
        return np.random.lognormal(mean, tmp_sigma, length)

    def set_R_lognorm_states(self, chain_length, mean, variance, length = 1, write_draws = True):
        if self.lognorm_mean == mean and self.lognorm_variance == variance and len(self.draws) == length:
            Rs = self.draws
        else:
            tmp_sigma = np.sqrt(np.log(variance ** 2 / mean ** 2 + 1))
            tmp_mean = np.log(mean) - 1 / 2 * tmp_sigma ** 2
            Rs =  np.random.lognormal(tmp_mean, tmp_sigma, length)
            if write_draws == True:
                self.lognorm_variance = variance
                self.lognorm_mean = mean
                self.draws = Rs
        return np.array([[draw] * (chain_length) for draw in Rs])