#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 13:56:50 2019

@author: chrx
"""
from __future__ import print_function

import fenics as fcs
import BC

class MySolver:
    pass
class PoissonSolver(MySolver):
    """
    class to solve the stationary diffusion equation in 2D/3D with nonlinear boundary condition
    
    Attributes
    ----------
    solver: fenics.NonlinearVariationalSolver
        instance of fenics nonlinear solver 
    p: Dict[str,float]
        stores modell parameters: R,q,k_on,D,rho
    u: fenics.MeshFunction
        stores solution
    neumann: List
        stores neumann BCs; set by self.compileSolver() 
    dirichlet: List[fenics.DirichletBC]:
        stores dirichlet BCs; set by self.compileSolver()
    nonLinBC:
        stores nonlinear BCs; set by self.compileSolver()
    dim: {2,3}
        mesh dimension
    """
    
    
    
    def __init__(self):
        self.dirichlet = []
        self.integralBC = []
        self.subdomains = []
        self.dim = 2
        self.p = {}
        self.mesh = None
        self.boundary_markers = None
        self.fieldQuantity = ""
        """
        Initalizes required fields.
        
        Parameters
        ----------
        mesh:dolfin.Mesh 
            domain mesh; intended to be genrated by myMeshing.MeshGenerator
        subdomains: List[myEntites.MySubDomain]
            list of all subdomains (including outer boundary) to read out boundary conditions
        boundary_markers: fenics.MeshFunction
            stores numbering of boundary patches
        dim : {2,3}
            mesh dimension
            
        Returns
        -------
        None
        """
        
        #defines functionspace v over mesh with first order Lagrange elemnts
        super().__init__()
    def log(self):
        return {"fieldQuantity":self.fieldQuantity,
                "p":self.p
                }
    def compileSolver(self):
        """
        compiles Solver; must be run befor solve()
        Detailed:
            creates u,v, and V; 
            defines boundary condition from subdomains list
            defines boundary integral ds() from MeshFunction "boundary_markers"
            defines problem and instantiates fenics solver
            
        """
        self.E_u = fcs.FiniteElement("P",self.mesh.ufl_cell(),1)
        self.E_c = fcs.FiniteElement("R",self.mesh.ufl_cell(),0)
        
#        self.V = fcs.FunctionSpace(self.mesh,self.E_u * self.E_c)
        self.V = fcs.FunctionSpace(self.mesh,self.E_u)

        
        #define trialfunction u, testfunction v from functionspace V
        
#        m = fcs.Function(self.V)
        u = fcs.Function(self.V)
        m = u
#        (u,c) = fcs.split(m)
#        (v_u,v_c) = fcs.TestFunction(self.V)
        (v_u) = fcs.TestFunction(self.V)



#        u = fcs.Function(self.V)
#        v = fcs.TestFunction(self.V)
        
        #define constants
        f = fcs.Constant(0)
        
        #checks wether diffusioncoefficient is given as paramenter
        
        if "D" in self.p:
            D = fcs.Constant(self.p["D"])
        else:
            pass
        if "kd" in self.p:
            kd = fcs.Constant(self.p["kd"])
        else:
            pass

        #defines ds as object of fencis class Measure, to give piecwise boundary integral as ds(i) for piece i
        ds = fcs.Measure("ds", domain=self.mesh, subdomain_data=self.boundary_markers)
        
        # iterates over subdomain list to set boundary condition according "bcDict" field
        
        self.dirichlet = []
        self.integralBC = []
        for i in self.subdomains:
            e = i["entity"]
            patch = i["patch"]
            bc = e.getBC(self.fieldQuantity)
            if type(bc) == BC.DirichletBC:
                #Dirichlet BC
                self.dirichlet.append(bc.getBC(self.V,self.boundary_markers,patch))
            if type(bc) == BC.Integral:
                self.integralBC.append(bc.getBC(u)*v_u*ds(patch))
        
        #Defines variational form of poisson equation with boundary integrals as sums
        F_u = -D*(fcs.dot(fcs.grad(u), fcs.grad(v_u))*fcs.dx)  - u*kd*v_u*fcs.dx + f*v_u*fcs.dx + D*(sum(self.integralBC))
#        F_c = (c*v_u + u*v_c)*fcs.dx  - mean*v_c*fcs.dx
        F = F_u#+ F_c
        #Defines nonlinear variational problem as F == 0 with Dirichlet BC and Jacobian given as the Gateaux derivative of F
        # with respect to u
        problem = fcs.NonlinearVariationalProblem(F,m, self.dirichlet,J=fcs.derivative(F, m))
        
        #instantiates fenics solver class as a field of "PoissonSolver"
        self.solver = fcs.NonlinearVariationalSolver(problem)
        
        #sets field u to be trialfunction
        self.u = u
        
    def solve(self):
        """
        Wrapper around fenics solve()
        
        Returns
        -------
        u:fenics.MeshFunction
            solution
        """
        #calls fenics solver; renames u for proper vtk output and returns solution u
        self.solver.solve()
#        u,c  = self.m.split()
        self.u.rename(self.fieldQuantity,self.fieldQuantity)
#        self.u = u
        return self.u
    def __del__(self):
        pass