#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 15:34:07 2022

@author: bennett
"""
from numba import njit, jit, int64, float64
from numba.experimental import jitclass
import numpy as np
###############################################################################
data = [('N_ang', int64), 
        ('N_space', int64),
        ('M', int64),
        ('tfinal', float64),
        ('IC', float64[:,:,:]),
        ('x0', float64),
        ('source', int64[:]),
        ("source_type", int64[:]),
        ("uncollided", int64),
        ]
@jitclass(data)
class IC_func(object):
    def __init__(self, source, uncollided):
        self.source = source
        self.uncollided = uncollided
    def function(self, x):
        # if self.uncollided == True:
        return 0