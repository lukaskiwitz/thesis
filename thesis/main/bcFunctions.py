#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 19:42:20 2019

@author: lukas
"""
from thesis.main.my_debug import message


def cellBC(u, p, field_quantity, area=1):
    """
    Defines integral boundary conditions for the cell.
    """

    R = p["R"]
    q = p["q"]
    k_on = p["k_on"]
    D = p["D"]

    uptake = k_on * R * u

    if "bc_type" in p.keys():

        v = p["bc_type"]

        if v == "linear":
            pass
        elif v == "R_saturation":
            Kc = p["Kc"]
            uptake = k_on * R * Kc * u / (Kc + u)
        elif v == "patrick_saturation":
            k_off = p["k_off"]
            k_endo = p["k_endo"]  # 1/s
            KD = k_off / k_on

            try:
                uptake = k_endo * R * u / (KD + u)
            except TypeError:
                KD = float(KD)
                k_endo = float(k_endo)
                uptake = k_endo * R * u / (KD + u)

        elif v == "amax_saturation":
            Kc = p["Kc"]
            amax = p["amax"]
            uptake = amax * u / (Kc + u)
        else:
            raise Exception
    else:
        message("bc type not in parameters. Using linear boundary condition")

    secretion = q
    linear = secretion / (D * area)
    billinear = -uptake / (D * area)

    return linear + billinear
