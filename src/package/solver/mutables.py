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
        ('x', float64[:])
        ]
@jitclass(data)
class IC_func(object):
    def __init__(self, source_type, uncollided, x0):
        self.source_type = source_type
        self.uncollided = uncollided
        self.x0 = x0
        
    def function(self, x):
        if self.uncollided == True:
            return np.zeros(x.size)
        elif self.uncollided == False and self.source_type[1] == 1:
            return self.square_IC(x)
        
    def square_IC(self, x):
        temp = np.greater(x, -self.x0)*1.0 - np.greater(x, self.x0)*1.0
        return temp/2.0
        
