#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 19:42:20 2019

@author: lukas
"""

def cellBC(u,p):
    """
    Defines the flux boundary conditions for the cell.
    Can be passed to the solver in the "Rec" field of the "bcDict" dictionary
    """
    R = p["R"]
    q = p["q"]
    k_on = p["k_on"]
    
    return (q-u*k_on*R)
def outerBC(u,p):
    """
    Defines the flux boundary condition on the outer boundary.
    Can be passed to the solver in the "Rec" field of the "bcDict" dictionary
    """
    k_on = p["k_on"]
    rho = p["rho"]
    L = p["L"]
    N = p["N"]
    R_resp = p["R_resp"]
    
    return k_on*u*(L+rho)*N*R_resp