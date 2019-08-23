#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 12:18:02 2019

@author: lukas
"""

import matplotlib.pyplot as plt
import numpy as np
import json
import os
import pandas as pd
import decimal




data = []
loc = "./logs/"
for i,o in enumerate(os.listdir(loc)):
        with open(loc+str(o),"r") as file:
            d = file.read()
            data.append((float(o),json.loads(d)))
data.sort()

