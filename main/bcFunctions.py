#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 19:42:20 2019

@author: lukas
"""

import fenics as fcs
import numpy as np
def cellBC_il2(u,p):
    """
    Defines the flux boundary conditions for the cell.
    Can be passed to the solver in the "Rec" field of the "bcDict" dictionary
    """
    R = p["R_il2"]
    q = p["q_il2"]
#    q = p["q"]
    k_on = p["k_on"]
    D = fcs.Constant(p["D"])
    R = fcs.Constant(R)
    q = fcs.Constant(q)
    k_on = fcs.Constant(k_on)
    a = fcs.Constant(4*np.pi*p["rho"]**2)
    
    
    return (q-u*k_on*R)/(D*a)

def cellBC_il6(u,p):
    """
    Defines the flux boundary conditions for the cell.
    Can be passed to the solver in the "Rec" field of the "bcDict" dictionary
    """
    R = p["R_il6"]
    q = p["q_il6"]
#    q = p["q"]
    k_on = p["k_on"]
    D = fcs.Constant(p["D"])
    R = fcs.Constant(R)
    q = fcs.Constant(q)
    k_on = fcs.Constant(k_on)
    a = fcs.Constant(4*np.pi*p["rho"]**2)
    
    return (q-u*k_on*R)/(D*a)


def cellBC_infg(u, p):
    """
    Defines the flux boundary conditions for the cell.
    Can be passed to the solver in the "Rec" field of the "bcDict" dictionary
    """
    R = p["R_infg"]
    q = p["q_infg"]
    k_on = p["k_on"]
    D = fcs.Constant(p["D"])
    R = fcs.Constant(R)
    q = fcs.Constant(q)
    k_on = fcs.Constant(k_on)
    a = fcs.Constant(4 * np.pi * p["rho"] ** 2)

    return (q - u * k_on * R) / (D * a)

def outerBC_il2(u,p):
    """
    Defines the flux boundary condition on the outer boundary.
    Can be passed to the solver in the "Rec" field of the "bcDict" dictionary
    """

    k_on = fcs.Constant(p["k_on"])
    D = fcs.Constant(p["D"])
    R = fcs.Constant(p["R_il2_b"])
    q = fcs.Constant(p["q_il2_b"])
    a  = fcs.Constant(4*np.pi*p["rho"]**2)
    
    return (q-u*k_on*R)/(D*a)

def outerBC_il2_unitTest(u,p):
   

    k_on = fcs.Constant(p["k_on"])
    D = fcs.Constant(p["D"])
    R = fcs.Constant(p["R_il2_N"])
    q = fcs.Constant(p["q_il2"])
    a  = fcs.Constant(4*np.pi*p["radius"]**2)
    
    return (q-u*k_on*R)/(D*a)

def outerBC_il6(u,p):
    """
    Defines the flux boundary condition on the outer boundary.
    Can be passed to the solver in the "Rec" field of the "bcDict" dictionary
    """
    k_on = fcs.Constant(p["k_on"])
    D = fcs.Constant(p["D"])

    R = fcs.Constant(p["R_il6_b"])
    q = fcs.Constant(p["q_il6_b"])
    a = fcs.Constant(4*np.pi*p["rho"]**2)
    
    
    return (q-u*k_on*R)/(D*a)

    