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
    Kc = p["Kc"]
    amax = p["amax"]

    uptake = k_on * R * u

    if "bc_type" in p.keys():

        v = p["bc_type"]

        if v == "linear":
            pass
        elif v == "R_saturation":
            uptake = k_on * R * Kc * u / (Kc + u)
        elif v == "amax_saturation":
            uptake = amax * u / (Kc + u)
        else:
            raise Exception
    else:
        message("bc type not in parameters. Using linear boundary condition")

    secretion = q
    linear = secretion / (D * area)
    billinear = -uptake / (D * area)

    return linear + billinear
