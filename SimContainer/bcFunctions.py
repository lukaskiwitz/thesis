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
    R = p["R_il2"]/(4*np.pi*p["rho"]**2)
    q = p["q_il2"]/(4*np.pi*p["rho"]**2)
#    q = p["q"]
    k_on = p["k_on"]
    D = fcs.Constant(p["D"])
    R = fcs.Constant(R)
    q = fcs.Constant(q)
    k_on = fcs.Constant(k_on)
    
    
    return (q-u*k_on*R)/D

def cellBC_il6(u,p):
    """
    Defines the flux boundary conditions for the cell.
    Can be passed to the solver in the "Rec" field of the "bcDict" dictionary
    """
    R = p["R_il6"]/(4*np.pi*p["rho"]**2)
    q = p["q_il6"]/(4*np.pi*p["rho"]**2)
#    q = p["q"]
    k_on = p["k_on"]
    
    R = fcs.Constant(R)
    q = fcs.Constant(q)
    k_on = fcs.Constant(k_on)
    
    
    return (q-u*k_on*R)

def outerBC_il2(u,p):
    """
    Defines the flux boundary condition on the outer boundary.
    Can be passed to the solver in the "Rec" field of the "bcDict" dictionary
    """

    k_on = fcs.Constant(p["k_on"])
    D = fcs.Constant(p["D"])
    R = fcs.Constant(p["R_il2"])/(4*np.pi*p["rho"]**2)
    q = fcs.Constant(p["q_il2"]/(4*np.pi*p["rho"]**2))
    rd = R*fcs.Constant(1e0) #receptor density
    return fcs.Expression("near(x[0],-dd) ? 1 : 0",degree=1,dd=p["dd"])*(q-u*k_on*rd)/D

def outerBC_il6(u,p):
    """
    Defines the flux boundary condition on the outer boundary.
    Can be passed to the solver in the "Rec" field of the "bcDict" dictionary
    """
    k_on = fcs.Constant(p["k_on"])
    
    R = fcs.Constant(p["R_il6"])/(4*np.pi*p["rho"]**2)
    q = fcs.Constant(p["q_il6"]/(4*np.pi*p["rho"]**2))
    
    rd = R*fcs.Constant(1e0) #receptor density
    return fcs.Expression("near(x[0],-dd) ? 1 : 0",degree=1,dd=p["dd"])*(q-u*k_on*rd)
    