
                
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
from numba import types, typed
import numba as nb

build_type = deferred_type()
build_type.define(build.class_type.instance_type)
kv_ty = (types.int64, types.unicode_type)
params_default = nb.typed.Dict.empty(key_type=nb.typeof('par_1'),value_type=nb.typeof(1))




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
        ('cv0', float64),
        ('fudge_factor', float64[:]),
        ('clight', float64),
        ('test_dimensional_rhs', int64),
        ('save_derivative', int64),
        ('xs_points', float64[:]),
        ('e_points', float64[:]),
        ('thermal_couple', nb.typeof(params_default)),


        ]
###############################################################################

@jitclass(data)
class T_function(object):
    def __init__(self, build):
        self.temp_function = np.array(list(build.temp_function), dtype = np.int64) 
        self.H = np.zeros(build.M+1).transpose()
        self.M = build.M
        
        self.a = 0.0137225 # GJ/cm$^3$/keV$^4
        self.alpha = 4 * self.a
        self.clight = 299.98
        
        self.xs_quad = build.xs_quad
        self.ws_quad = build.ws_quad
        self.cv0 = build.cv0 / self.a 
        if (self.cv0) != 0.0:
            print('cv0 is ', self.cv0)
        self.test_dimensional_rhs = False
        self.save_derivative = build.save_wave_loc
        self.thermal_couple = build.thermal_couple

        
    def make_e(self, xs, a, b):
        temp = xs*0
        for ix in range(xs.size):
            for j in range(self.M+1):
                temp[ix] += normPn(j, xs[ix:ix+1], a, b)[0] * self.e_vec[j] 
        return temp 
    
    def integrate_quad(self, a, b, j):
        if self.thermal_couple['none'] != True:
            argument = (b-a)/2 * self.xs_quad + (a+b)/2
            self.H[j] = (b-a)/2 * np.sum(self.ws_quad * self.T_func(argument, a, b) * normPn(j, argument, a, b))
        else:
             self.H[j] = 0
        
    def T_func(self, argument, a, b):
        e = self.make_e(argument, a, b)

        self.xs_points = argument
        self.e_points = e
        if self.temp_function[0] == 1:
            T = self.su_olson_source(e, argument, a, b)
            return self.a * np.power(T,4) * self.fudge_factor

        elif self.temp_function[1] == 1:
            if self.test_dimensional_rhs == True:
                T =  e / self.cv0
                return np.power(T,4) * self.a * self.clight
            else:
                T =  e / self.cv0 
                return np.power(T,4)
        else:
            assert(0)

        

        
    def su_olson_source(self, e, x, a, b):
        
        self.fudge_factor = np.ones(e.size)
    
        for count in range(e.size):
            if math.isnan(e[count]) == True:
                            print("nan")
                            print(e)
                            assert 0     
            elif (e[count]) < 0.:
                self.fudge_factor[count] = -1.


        t1 = np.abs(4*e/self.alpha)
        return np.power(t1,0.25)
        
        
    def make_H(self, xL, xR, e_vec):
        self.e_vec = e_vec
            
        for j in range(self.M+1):
            self.integrate_quad(xL, xR, j)

        
