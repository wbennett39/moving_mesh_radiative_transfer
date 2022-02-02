#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 14:38:34 2022

@author: bennett
"""
import numpy as np
from build_problem import build
import math
###############################################################################
class source_class:
    def __init__(self, build):
        self.S = np.zeros(build.M+1).transpose()
        self.source_type = build.source_type
        self.uncollided = build.uncollided
        self.x0 = build.x0
        
    def square_ic_func(xs, t):
        # this function is broken
        temp = xs*0
        # temp += np.where(np.less_equal(np.abs(xs), t + x0) or np.less_equal(np.abs(xs), t - x0), math.exp(-t)/(2*t + build.x0), 0) 
        # temp += np.where(np.greater_equal(x0 - xs, t) and np.less_equal(x0 + xs, t), math.exp(-t)*(t + xs + x0)/(4*t*t*x0), 0)
        # temp += np.where(np.less_equal(x0 - xs, t) and np.greater_equal(x0 + xs, t), math.exp(-t)*(t - xs + x0)/(4*t*t*x0), 0)
        return temp
        
        
    def __call__(self, t, xR, xL):
        if self.uncollided == True:
            if self.source_type == "plane":
                self.S[0] = (math.exp(-t)/(2*t+self.x0)*math.sqrt(xR-xL))
                return self.S
            elif self.source_type == "square_IC":
                return "placeholder"
                
            