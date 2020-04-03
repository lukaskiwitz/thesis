#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 19:42:20 2019

@author: lukas
"""

import fenics as fcs
import numpy as np
from my_debug import message
def cellBC(u,p, field_quantity, area = 1):
    """
    Defines the flux boundary conditions for the cell.
    """

    R = p.get_physical_parameter_by_field_quantity("R",field_quantity).get_in_sim_unit()
    q = p.get_physical_parameter_by_field_quantity("q",field_quantity).get_in_sim_unit()
    k_on = p.get_physical_parameter_by_field_quantity("k_on",field_quantity).get_in_sim_unit()

    D = fcs.Constant(p.get_physical_parameter_by_field_quantity("D",field_quantity).get_in_sim_unit())

    R = fcs.Constant(R)
    q = fcs.Constant(q)
    k_on = fcs.Constant(k_on)
    a = fcs.Constant(area)
    
    
    return (q-u*k_on*R)/(D*a)

def outerBC_il2(u,p , area = 1):
    """
    Defines the flux boundary condition on the outer boundary.
    """

    R = p.get_physical_parameter("R", "IL-2").get_in_sim_unit()
    q = p.get_physical_parameter("q", "IL-2").get_in_sim_unit()
    k_on = p.get_physical_parameter("k_on", "IL-2").get_in_sim_unit()
    D = fcs.Constant(p.get_physical_parameter("D", "IL-2").get_in_sim_unit())
    R = fcs.Constant(R)
    q = fcs.Constant(q)
    k_on = fcs.Constant(k_on)
    a = fcs.Constant(area)
    
    return (q-u*k_on*R)/(D*a)

def outerBC_il2_unitTest(u,p):
   

    k_on = fcs.Constant(p["k_on"])
    D = fcs.Constant(p["D"])
    R = fcs.Constant(p["R_il2_N"])
    q = fcs.Constant(p["q_il2"])
    a  = fcs.Constant(4*np.pi*p["radius"]**2)
    
    return (q-u*k_on*R)/(D*a)
    