#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 12:49:44 2022

@author: bennett
"""
from .benchmark_functions import F1, F1_spacefirst, find_intervals_time
import scipy.integrate as integrate
import math
import numpy as np
###############################################################################

def opts0(*args, **kwargs):
       return {'limit':10000000,'epsabs': 1.0e-11, 'epsrel':1.0e-11}
   
def opts1(*args, **kwargs):
       return {'limit':10000000,'epsabs': 1.0e-11, 'epsrel':1.0e-11}
def opts2(*args, **kwargs):
       return {'limit':10000000,'epsabs': 1.0e-11, 'epsrel':1.0e-11}

class collided:
    def __init__(self, source_type, x0, t0):
        self.x0 = x0
        self.t0 = t0
        self.source_type = source_type
    
    def plane_IC(self, xs,t):
        temp = xs*0
        for ix in range(xs.size):
            temp[ix] = integrate.nquad(F1, [[0, math.pi]], args =  (0.0, 0.0, xs[ix], t, 0), opts = [opts0])[0]
        return temp
    
    def square_IC(self, xs, t):
        temp = xs*0
        for ix in range(xs.size):
            temp[ix] = integrate.nquad(F1, [[0, math.pi], [-self.x0, self.x0]], args =  (0.0, xs[ix], t, 0), opts = [opts0, opts1])[0]
        return temp
    
    def source_double_integral_time(self, s, x, t, source):
        ab = find_intervals_time(t, x, s)
        solution = integrate.nquad(F1_spacefirst, [[0, math.pi], [ab[0],ab[1]]], args =  (s, x, t, source), opts = [opts0, opts1])[0]
        return solution
        
    def square_source(self, xs, t):
        temp = xs*0
        for ix in range(xs.size):
            temp[ix] = integrate.nquad(self.source_double_integral_time, [[-self.x0, self.x0]], args = (xs[ix], t, 0), opts = [opts2])[0]
        return temp
    
    def gaussian_IC(self, xs, t):
        temp = xs*0
        for ix in range(xs.size):
            temp[ix] = integrate.nquad(F1, [[0, math.pi], [-np.inf, np.inf]], args =  (0.0, xs[ix], t, 1), opts = [opts0, opts1])[0]
        return temp

    def gaussian_source(self, xs, t):
        temp = xs*0
        for ix in range(xs.size):
            temp[ix] = integrate.nquad(self.source_double_integral_time, [[-self.x0, self.x0]], args = (xs[ix], t, 1), opts = [opts2])[0]
        return temp
    
    def __call__(self, xs, t):
        if self.source_type == 'plane_IC':
            return self.plane_IC(xs, t)
        elif self.source_type == 'square_IC':
            return self.square_IC(xs, t)
        elif self.source_type == 'square_source':
            return self.square_source(xs, t)
        elif self.source_type == 'gaussian_IC':                
            return self.gaussian_IC(xs, t)
        elif self.source_type == 'gaussian_source':
            return self.gaussian_source(xs, t)