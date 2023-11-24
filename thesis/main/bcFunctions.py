#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 19:42:20 2019

@author: lukas
"""

import logging

from thesis.main.my_debug import message

module_logger = logging.getLogger(__name__)

def cellBC(I, p, field_quantity, area=1):
    """
    Defines integral boundary conditions for the cell.
    """
    R = p["R"]
    q = p["q"]
    k_on = p["k_on"]
    k_off = p["k_off"]
    D = p["D"]
    # print(k_on.__float__())
    try:
        KD = p["KD"]
    except:
        KD = k_off / k_on

    uptake = k_on * R * I

    if "bc_type" in p.keys():

        v = p["bc_type"]

        if v == "linear":
            pass
        elif v == "R_saturation":
            Kc = p["Kc"]
            uptake = k_on * R * Kc * I / (Kc + I)
        elif v == "patrick_saturation":
            k_endo = p["k_endo"]  # 1/s
            uptake = k_endo * R * I / (KD + I)
        elif v == "k_off_saturation":
            KD = k_off / k_on
            k_endo = p["k_endo"]  # 1/s
            uptake = k_endo * R * I / (KD + I)

        elif v == "amax_saturation":
            Kc = p["Kc"]
            amax = p["amax"]
            uptake = amax * I / (Kc + I)
        else:
            raise Exception
    else:
        message("bc type not in parameters. Using linear boundary condition", module_logger)

    secretion = q
    linear = secretion / (D * area)
    billinear = -uptake / (D * area)

    return linear + billinear
