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
        ('x', float64[:]),
        ('source_strength', float64),
        ('sigma', float64)
        ]
@jitclass(data)
class IC_func(object):
    def __init__(self, source_type, uncollided, x0, source_strength, sigma):
        self.source_type = source_type
        self.uncollided = uncollided
        self.x0 = x0
        self.source_strength = source_strength
        self.sigma = sigma


    def function(self, x):
        if self.uncollided == True:
            return np.zeros(x.size)
        elif self.uncollided == False and self.source_type[0] == 1:
            return self.plane_and_square_IC(x)/self.x0/2.0
        elif self.uncollided == False and self.source_type[1] == 1:
            return self.plane_and_square_IC(x)
        elif self.uncollided == False and self.source_type[2] == 1:
            return np.zeros(x.size)
        elif self.uncollided == False and self.source_type[3] == 1:
            return self.gaussian_IC(x)
        elif self.source_type[4] == 1:
            return self.MMS_IC(x)
        else:
            return np.zeros(x.size)
        
    def plane_and_square_IC(self, x):
        temp = (np.greater(x, -self.x0) - np.greater(x, self.x0))*self.source_strength
            # temp = x/x
        return temp/2.0
    
    def gaussian_IC(self, x):
        temp = np.exp(-x*x/self.sigma**2)*self.source_strength
        return temp/2.0
    
    def MMS_IC(self, x):
        # temp = np.greater(x, -self.x0)*1.0 - np.greater(x, self.x0)*1.0 * np.exp(-x*x/2)/(2)
        temp = np.exp(-x*x/2)/(2)
        return temp
        
        
