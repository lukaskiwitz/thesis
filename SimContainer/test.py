# -*- coding: utf-8 -*-

import FieldProblem as fp
import Entity
import MySolver

cell = Entity.Cell([0,0,0],0.5,{"D":1})
domain = Entity.DomainSphere([0,0,0],1,{"N":-1})
solver = MySolver.PoissonSolver()
fieldProblem = fp.FieldProblem()

fieldProblem.setSolver(solver)
fieldProblem.addEntity(cell)
fieldProblem.generateMesh(domain,16)
fieldProblem.step()


