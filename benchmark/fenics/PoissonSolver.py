#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 13:56:50 2019

@author: chrx
"""
from __future__ import print_function

import fenics as fcs

class PoissonSolver:
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
    neumann = []
    dirichlet = []
    nonLinBC = []
    subdomains = []
    dim = 2
    p = {}
    
    def __init__(self,mesh,subdomains,boundary_markers,dim):
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
        self.mesh = mesh
        self.dim = dim
        self.subdomains = subdomains
        self.boundary_markers = boundary_markers
        
        #defines functionspace v over mesh with first order Lagrange elemnts
        self.V = fcs.FunctionSpace(self.mesh,"P",1)
        
    def compileSolver(self):
        """
        compiles Solver; must be run befor solve()
        Detailed:
            creates u,v, and V; 
            defines boundary condition from subdomains list
            defines boundary integral ds() from MeshFunction "boundary_markers"
            defines problem and instantiates fenics solver
            
        """
        
        #define trialfunction u, testfunction v from functionspace V
        u = fcs.Function(self.V)
        v = fcs.TestFunction(self.V)
        
        #define constants
        f = fcs.Constant(0)
        
        #checks wether diffusioncoefficient is given as paramenter
        if "D" in self.p:
            D = fcs.Constant(self.p["D"])
        else:
            D = fcs.Constant(1)
            
        #defines ds as object of fencis class Measure, to give piecwise boundary integral as ds(i) for piece i
        ds = fcs.Measure("ds", domain=self.mesh, subdomain_data=self.boundary_markers)
        
        # iterates over subdomain list to set boundary condition according "bcDict" field
        self.neumann = []
        self.dirichlet = []
        self.nonLinBC = []
        for i in self.subdomains:
            if "D" in i.bcDict:
                #Dirichlet BC
                value = i.bcDict["D"]
                value = fcs.Expression(str(value),degree=2)
                bc = fcs.DirichletBC(self.V, value ,self.boundary_markers,i.patch)
                self.dirichlet.append(bc)
            if "N" in i.bcDict:
                #Neumann BC
                value = i.bcDict["N"]
                value = fcs.Expression(str(value),degree=2)
                self.neumann.append(value*v*ds(i.patch))
            if "Rec" in i.bcDict:
                #nonlinear Secretion according to Function object q(u,p) passed in dict entry
                q = i.bcDict["Rec"]
                p = i.bcDict["p"] if "p" in i.bcDict else self.p
                self.nonLinBC.append(q(u,p)*v*ds(i.patch))
        #Defines variational form of poisson equation with boundary integrals as sums
        F = -D*(fcs.dot(fcs.grad(u), fcs.grad(v))*fcs.dx) + f*v*fcs.dx + D*(sum(self.neumann) - sum(self.nonLinBC))
        
        #Defines nonlinear variational problem as F == 0 with Dirichlet BC and Jacobian given as the Gateaux derivative of F
        # with respect to u
        problem = fcs.NonlinearVariationalProblem(F,u, self.dirichlet,J=fcs.derivative(F, u))
        
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
        
        self.u.rename('u','u')
        return self.u
    def __del__(self):
        pass