
                
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  5 10:42:11 2022

@author: bennett
"""

from .build_problem import build
from .functions import normPn

from numba.experimental import jitclass
from numba import int64, float64, deferred_type, prange
import numpy as np
import math

build_type = deferred_type()
build_type.define(build.class_type.instance_type)



data = [('temp_function', int64[:]),
        ('e_vec', float64[:]),
        ('e', float64),
        ('H', float64[:]),
        ('alpha', float64),
        ('a', float64),
        ('M', int64),
        ("xs_quad", float64[:]),
        ("ws_quad", float64[:]),
        ("T", float64[:]),
        ('cv0', float64)
        ]
###############################################################################

@jitclass(data)
class T_function(object):
    def __init__(self, build):
        self.temp_function = build.temp_function
        self.H = np.zeros(build.M+1).transpose()
        self.M = build.M
        
        self.a = 0.0137225 # GJ/cm$^3$/keV$^4
        self.alpha = 4 * self.a
        
        self.xs_quad = build.xs_quad
        self.ws_quad = build.ws_quad
        self.cv0 = build.cv0
        if self.cv0 != 0:
            print('cv0 is ', self.cv0)

        
        
    def make_e(self, xs, a, b):
        temp = xs*0
        for ix in range(xs.size):
            for j in range(self.M+1):
                temp[ix] += normPn(j, xs[ix:ix+1], a, b)[0] * self.e_vec[j] 
        return temp 
    
    def integrate_quad(self, a, b, j):
        argument = (b-a)/2 * self.xs_quad + (a+b)/2
        self.H[j] = (b-a)/2 * np.sum(self.ws_quad * self.T_func(argument, a, b) * normPn(j, argument, a, b))
        
    def T_func(self, argument, a, b):
        if self.temp_function[0] == 1:
            T = self.su_olson_source(argument, a, b)
            return self.a * np.power(T,4)
        elif self.temp_function[1] == 1:
            e = self.make_e(argument,a,b)
            T =  e / self.cv0
            return self.a * np.power(T,4)
        else:
            assert(0)
        
    def su_olson_source(self, x, a, b):
        e = self.make_e(x, a, b)
        for count in range(e.size):
            if math.isnan(e[count]) == True:
                            print("nan")
                            print(e)
                            assert 0       
        t1 = np.abs(4*e/self.alpha)
        return np.power(t1,0.25)
        
        
    def make_H(self, xL, xR, e_vec):
        self.e_vec = e_vec
            
        for j in range(self.M+1):
            self.integrate_quad(xL, xR, j)