import numpy as np
from scipy.integrate import solve_ivp

def halftime_estimation(eta, tau, R_start = 1.5e3, gamma=100):
    def func(t, y, dict, dummy=0):
        R0, R1, R2, R3, R4, R5 = y

        decay = dict["eta"]

        dR0_dt = ((dict["k_min"] * dict["k"] ** dict["N"] + dict["k_max"] * dict["feedback_k"] ** dict["N"]) / (
                dict["k"] ** dict["N"] + dict["feedback_k"] ** dict["N"]) - decay * R0)

        dR1_dt = dict["lambda"] * R0 - dict["lambda"] * R1
        dR2_dt = dict["lambda"] * R1 - dict["lambda"] * R2
        dR3_dt = dict["lambda"] * R2 - dict["lambda"] * R3
        dR4_dt = dict["lambda"] * R3 - dict["lambda"] * R4
        dR5_dt = dict["lambda"] * R4 - dict["lambda"] * R5
        return [dR0_dt, dR1_dt, dR2_dt, dR3_dt, dR4_dt, dR5_dt]

    dict = {"N": 3, "gamma": gamma, "eta": eta, "il2": 1e5, "feedback_k": 1e5,
            "k_min": R_start * eta, "k_max": R_start * eta * gamma, "k": 0.5, "lambda": np.log(2) / (tau * 3600)}

    t1 = 0
    dt = 150 * 3600
    y0 = [R_start] * 6
    t_space = np.linspace(t1, t1 + dt, 20000)
    result = solve_ivp(func, t_span=[t1, t1 + dt], t_eval=t_space, y0=y0, method="LSODA", args=(dict, 0))

    plot_y = result.y[-1]
    half_value = (result.y[-1].max() - result.y[-1].min()) / 2 + result.y[-1].min()

    # print("half time (h):", t_space[np.argmin(np.abs(plot_y - half_value))] / 3600)
    return t_space[np.argmin(np.abs(plot_y - half_value))] / 3600