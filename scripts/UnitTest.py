#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 09:29:30 2019

@author: kiwitz
"""
import FieldProblem as fp
import Entity
import MySolver
import numpy as np
import matplotlib.pyplot as plt
import BC as bc
import SimContainer as SC
from bcFunctions import cellBC_il2,outerBC_il2, outerBC_il2_unitTest
from copy import deepcopy
import os
import time
import math
from scipy.constants import N_A
def u_exact_1(x,p):
    R = p["R_il2"]
    q = p["q_il2"]
    k_on = p["k_on"]
    rho = p["rho"]
    D= p["D"]
    
    return q*rho/(np.linalg.norm(x)*(k_on*R+D*4*np.pi*rho))

def du_dr_exact_1_(x,p):
    R = p["R_il2"]
    q = p["q_il2"]
    k_on = p["k_on"]
    rho = p["rho"]
    D= p["D"]
    
    return (1/math.log(np.linalg.norm(x)))*q*rho/((k_on*R+D*4*np.pi*rho))

def u_exact_2(x,p):
    """
    analytic solution for plotting with matplotlib
    """
    R = p["R_il2"]
    q = p["q_il2"]
    k_on = p["k_on"]
    rho = p["rho"]
    D= p["D"]
    r = np.linalg.norm(x)
    L = p["L"]
    N = p["N"]
    R_resp = p["R_resp"]
    
    
    A = q*rho/(k_on*r)
    num = 4*np.pi*D*r*(L + rho) + k_on*(L-r+rho)*N*R_resp
    denom = k_on*L*R*N*R_resp + 4*np.pi*D*rho*(L+rho)*(R+N*R_resp)
    return A * (num/denom)

def test1(p_global):
    path = "/extra/kiwitz/results/unit_test/test1/"
    p = p_global
    
    p_domain = deepcopy(p_global)
    p_domain.update({
             "R_il2":0,
             "q_il2":0})
    p_cell = deepcopy(p_global)
    p_cell.update(
            {"R_il2":p_global["high"],
             "q_il2":p_global["q_high"]
            }
            )
    domain = Entity.DomainSphere([0,0,0],p_global["radius"],p_domain,[
#           bc.outerDirichletBC(u_exact_1(p_global["radius"],p),"on_boundary",fieldQuantity="il2"),
            bc.outerDirichletBC("0","on_boundary",fieldQuantity="il2")
#           bc.outerIntegral(lambda u,p:1*u,"on_boundary",fieldQuantity="il2")
            ])
    
    """IL-2"""
#    solver_il2 = MySolver.PoissonSolver()
#    
#    fieldProblem_il2 = fp.FieldProblem()
#    fieldProblem_il2.fieldName = "il2"
#    fieldProblem_il2.fieldQuantity = "il2"
##    fieldProblem_il2.extCache = "/extra/kiwitz/results/unit_test/test1/meshCache_il2"
#    
#    
#    fieldProblem_il2.setSolver(solver_il2)
#    fieldProblem_il2.p = deepcopy(p_global)
#    fieldProblem_il2.setOuterDomain(domain)
#    
#    sc = SC.SimContainer()
#    
#    
#    cell = Entity.Cell([0,0,0],0.05,
#       [
#           bc.Integral(cellBC_il2,fieldQuantity="il2")
#       ])
#    cell.name = "cell"
#    cell.p = p_cell
#    sc.addEntity(cell)
#    
#    sc.addField(fieldProblem_il2)
#    
#    sc.path = path
#    sc.initialize(max_char_length=0.05,min_char_length=0.0001)
#    
#    
#    for n,i in enumerate(range(1)):
#    
#        start = time.process_time()
#        sc.step(1)
#        end = time.process_time()
#        print("time: "+str(end-start)+"s for step number "+str(n))
#        sc.saveFields(n)
#    
#    u_1 = sc.fields[0].getField()
    x_1= np.linspace(p_global["rho"],p_global["radius"],100)
#    x_1= np.linspace(p_global["rho"],0.2,100)
#    y_1 = [10**9*u_1([0,0,i]) for i in x_1]
    y_1_e = [10**9*u_exact_1(i,p_cell) for i in x_1]
    
#    r_1 = np.abs(np.subtract(y_1,y_1_e))
    
    plt.figure(1)
#    plt.plot(x_1,y_1,color="tab:red")
    plt.plot(x_1,y_1_e,"--",color="tab:red")
#    plt.plot(x_1,r_1,color="tab:orange")
#    plt.xlim(0,0.2)
    plt.legend(["FEM_low","exact solution","error"])
    
    tk = np.arange(0,x_1[-1]-x_1[0]+0.25,0.25)
    
    plt.xticks(tk+0.05,np.round(tk*100))
    plt.xlabel("distance from cell boundary [$\mu m$]")
    plt.ylabel("[IL-2][nM]")    
#    plt.savefig(path+"plot.png")
    
def test3(p_global):
    path = "/extra/kiwitz/results/unit_test/test1/"
    p = p_global
    
    p_domain = deepcopy(p_global)
    p_domain.update({
             "R_il2":0,
             "q_il2":0})
    p_cell = deepcopy(p_global)
    p_cell.update(
            {"R_il2":p_global["high"],
             "q_il2":p_global["q_high"]
            }
            )
    domain = Entity.DomainSphere([0,0,0],p_global["radius"],p_domain,[
#           bc.outerDirichletBC("0","on_boundary",fieldQuantity="il2"),
           bc.outerIntegral(lambda u,p:-u,"on_boundary",fieldQuantity="il2")
            ])
    
    """IL-2"""
    solver_il2 = MySolver.PoissonSolver()
    
    fieldProblem_il2 = fp.FieldProblem()
    fieldProblem_il2.fieldName = "il2"
    fieldProblem_il2.fieldQuantity = "il2"
    fieldProblem_il2.extCache = "/extra/kiwitz/results/unit_test/test1/meshCache_il2"
    
    
    fieldProblem_il2.setSolver(solver_il2)
    fieldProblem_il2.p = deepcopy(p_global)
    fieldProblem_il2.setOuterDomain(domain)
    
    sc = SC.SimContainer()
    
    
    cell = Entity.Cell([0,0,0],0.05,
       [
           bc.Integral(cellBC_il2,fieldQuantity="il2")
       ])
    cell.name = "cell"
    cell.p = p_cell
    sc.addEntity(cell)
    
    sc.addField(fieldProblem_il2)
    
    sc.path = path
    sc.initialize(max_char_length=0.05,min_char_length=0.0001)
    
    
    for n,i in enumerate(range(1)):
    
        start = time.process_time()
        sc.step(1)
        end = time.process_time()
        print("time: "+str(end-start)+"s for step number "+str(n))
        sc.saveFields(n)
    
    u_1 = sc.fields[0].getField()
    x_1= np.linspace(p_global["rho"],p_global["radius"],100)
#    x_1= np.linspace(p_global["rho"],0.2,100)
    y_1 = [10**9*u_1([0,0,i]) for i in x_1]
    y_1_e = [10**9*u_exact_1(i,p_cell) for i in x_1]
    
    r_1 = np.abs(np.subtract(y_1,y_1_e))
    
    plt.figure(1)
    plt.plot(x_1,y_1,color="tab:blue")
#    plt.plot(x_1,y_1_e,"--",color="tab:red")
#    plt.plot(x_1,r_1,color="tab:orange")
#    plt.xlim(0,1)
    plt.ylim(0,0.15)
    plt.legend([r'$u=0$',"exact solution","robin"])
    
#    tk = np.arange(0,x_1[-1]-x_1[0]+0.25,0.25)
    tk = np.arange(0,x_1[-1]-x_1[0],0.25)
    
    plt.xticks(tk+0.05,np.round(tk*100))
    plt.xlabel("distance from cell boundary [$\mu m$]")
    plt.ylabel("[IL-2][nM]")    
    plt.savefig(path+"plot.png")
   
def test2(p_global):
    path = "/extra/kiwitz/results/unit_test/test2/"
 
    p_global["L"] = p_global["radius"]
    
    p_domain = deepcopy(p_global)
    p_domain.update({
             "R_il2":p_global["R_resp"],
             "q_il2":0})
    
    domain = Entity.DomainSphere([0,0,0],p_global["radius"],p_domain,[
#           bc.outerIntegral(lambda u,p: fcs.Constant(1),"true",fieldQuantity="il2"),
           bc.outerIntegral(outerBC_il2_unitTest,"1",fieldQuantity="il2")
            ])
    
    """IL-2"""
#    solver_il2 = MySolver.PoissonSolver()
#    
#    fieldProblem_il2 = fp.FieldProblem()
#    fieldProblem_il2.fieldName = "il2"
#    fieldProblem_il2.fieldQuantity = "il2"
#    
#    
#    fieldProblem_il2.setSolver(solver_il2)
#    fieldProblem_il2.p = deepcopy(p_global)
#    fieldProblem_il2.setOuterDomain(domain)
#    
#    sc = SC.SimContainer()
#    
    p_cell = deepcopy(p_global)
    p_cell.update(
            {"R_il2":p_global["high"],
             "q_il2":p_global["q_high"]
            }
            )
#    cell = Entity.Cell([0,0,0],0.05,
#       [
#           bc.Integral(cellBC_il2,fieldQuantity="il2")
#       ])
#    cell.name = "cell"
#    cell.p = p_cell
#    sc.addEntity(cell)
#    
#    sc.addField(fieldProblem_il2)
#    
#    sc.path = path
#    sc.initialize(max_char_length=0.001,min_char_length=0.000001)
#    
#    
#    for n,i in enumerate(range(1)):
#    
#        start = time.process_time()
#        sc.step(1)
#        end = time.process_time()
#        print("time: "+str(end-start)+"s for step number "+str(n))
#        sc.saveFields(n)
#    u_2 = sc.fields[0].getField()
    x_2 = np.linspace(p_global["rho"],p_global["radius"],100)
#    
#    
#    
#    y_2 = [10**9*u_2([0,0,i]) for i in x_2]
    y_2_e = [10**9*u_exact_2(i,p_cell) for i in x_2]
    
#    r_2 = np.abs(np.subtract(y_2,y_2_e))
    
    
#    plt.figure(2)
#    plt.plot(x_2,y_2,color="tab:blue")
    plt.plot(x_2,y_2_e,"--",color="tab:blue",)
#    plt.plot(x_2,r_2,color="tab:orange")
    plt.legend(["FEM_low","FEM_high","exact solution"])
    tk = np.arange(0,x_2[-1]-x_2[0]+0.05,0.05)
    plt.xticks(tk+0.05,np.round(tk*100))
    plt.xlabel("distance from cell boundary [$\mu m$]")
    plt.ylabel("[IL-2][nM]")
    os.makedirs(path,exist_ok=True)
    plt.savefig(path+"plot.png")
p_global = {
             "R_il2":0,
             "q_il2":0,
             "k_on": 10e9*111.6/60**2,#111.6 per hour
             "rho": 0.05,#mu
             "D":(10**0.5*0.01)**2,#mu² per s
             "high":100*N_A**-1*10e9,
             "low":1*N_A**-1*10e9,
             "kd":0,#0.1/(60*2),
             "q_high":10*N_A**-1*10e9,
             "radius":1,
             "N":8,
#             "x":x,
#               "y":y,
#               "z":z,
               "d_x":1,
               "d_y":1,
               "d_z":1
            }

"""u_1"""
test1(p_global)




"""u_2"""
p_global.update({"radius":0.2,
                 "L":0.2,
                 "R_il2_N":p_global["N"]*p_global["low"],
                 "R_resp":p_global["low"]
                 })


test2(p_global)

