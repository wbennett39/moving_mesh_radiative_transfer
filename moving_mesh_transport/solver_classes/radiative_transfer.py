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
        ("ws_quad", float64[:])
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
        
        
    def integrate_quad(self, a, b, j, func):
        argument = (b-a)/2 * self.xs_quad + (a+b)/2
        self.H[j] = (b-a)/2 * np.sum(self.ws_quad * self.a * func(argument, a, b)**4)
            
    def make_e(self, xs, a, b):
        temp = xs*0
        for j in range(self.M+1):
            temp += normPn(j, xs, a, b) * self.e_vec[j]
        return temp 
    
    def su_olson_source(self, x, xL, xR):
        temp = 0.0
        e = self.make_e(x, xL, xR)
        temp = 4 * e**0.25 / self.alpha
        return temp
        
    def make_H(self, xL, xR, e_vec):
        self.e_vec = e_vec
        
        if self.temp_function[0] == 1:
            for j in range(self.M+1):
                self.integrate_quad(xL, xR, j, self.su_olson_source)