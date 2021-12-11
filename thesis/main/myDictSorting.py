#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 09:52:45 2019

@author: kiwitz
"""

import itertools
import logging

module_logger = logging.getLogger(__name__)

def sortDict(x, keys):
    res = x
    for i in keys:
        res = res[i]
    return res


def groupByKey(ls, keys):
    ls.sort(key=lambda x: sortDict(x, keys))
    l = []
    for i, e in itertools.groupby(ls, key=lambda x: sortDict(x, keys)):
        l.append(list(e))
    return l
