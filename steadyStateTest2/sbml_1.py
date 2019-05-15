#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import tellurium as te

r = te.loada("""
             model test
             RK + C -> RK + C_i; vmax * C/(Km + C)             
             RK = 1
             C = 1
             vmax = 1
             Km = 1
             end
             """)
result = r.simulate(0,10,100)
r.plot()

sbml_file = open('./sbml.xml',"w")
sbml_file.write(r.getSBML())
sbml_file.close()