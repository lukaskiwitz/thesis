#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 12:06:52 2019

@author: Lukas Kiwitz
"""

import fenics as fcs
import numpy as np
from dolfin import near, Point
from mshr import Circle, Sphere, Rectangle, Box


def insideRectangle(p1, p2, x):
    x1 = p1[0]
    x2 = p2[0]
    y1 = p1[1]
    y2 = p2[1]
    tol = 10e-16
    return near(x[0], x1, tol) or near(x[0], x2, tol) or near(x[1], y1, tol) or near(x[1], y2, tol)


class MySubDomain(fcs.SubDomain):
    def __init__(self):
        self.patch = 0
        super().__init__()

    def inside(self, x, on_boundary):
        return on_boundary


class OuterBoundary(MySubDomain):
    def __init__(self, **kwargs):
        super().__init__()


class OuterSphere(OuterBoundary):
    center = [0, 0]
    radius = 0.5

    def __init__(self, c, r, **kwargs):
        self.center = c
        self.radius = r
        super().__init__(**kwargs)

    def getGeometry(self, dim):
        if dim == 2:
            return Circle(Point(self.center[0], self.center[1]), self.radius)
        elif dim == 3:
            return Sphere(Point(self.center[0], self.center[1], self.center[2]), self.radius)

    def inside(self, x, on_boundary):
        tol = 10e-2
        return on_boundary and near(np.linalg.norm(x - self.center), self.radius, tol)


class OuterCube(OuterBoundary):

    def __init__(self, p1, p2, **kwargs):
        self.p1 = p1
        self.p2 = p2
        super().__init__(**kwargs)

    def getGeometry(self, dim):
        if dim == 2:
            return Rectangle(Point(self.p1[0], self.p1[1]), Point(self.p2[0], self.p2[1]))
        elif dim == 3:
            return Box(Point(self.p1[0], self.p1[1], self.p1[2]), Point(self.p2[0], self.p2[1], self.p2[2]))

    def inside(self, x, on_boundary):

        if len(x) == 2:
            return on_boundary and insideRectangle(self.p1, self.p2, x)
        elif len(x) == 3:
            def near(a, b):
                return (a - b) < 10e-3

            p1 = self.p2
            p0 = self.p1

            return on_boundary and near(x[0], p1[0]) or near(x[0], p0[0]) or near(x[1], p1[1]) or near(x[1],
                                                                                                       p0[1]) or near(
                x[2], p1[2]) or near(x[2], p0[2])



class CellSubDomain(MySubDomain):

    def __init__(self, c, r):
        self.radius = r
        self.center = c
        super().__init__()

    def getGeometry(self, dim):
        if dim == 2:
            return Circle(Point(self.center[0], self.center[1]), self.radius)
        elif dim == 3:
            return Sphere(Point(self.center[0], self.center[1], self.center[2]), self.radius)

    def inside(self, x, on_boundary):
        tol = 10e-2
        #        message(str(x))
        return on_boundary and near(np.linalg.norm(x - self.center), self.radius, tol)
