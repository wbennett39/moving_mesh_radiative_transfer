#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 14:38:34 2022

@author: bennett
"""
import numpy as np
from build_problem import build
import math
from functions import normPn

from numba import float64, int64, deferred_type
from numba.experimental import jitclass
###############################################################################
build_type = deferred_type()
build_type.define(build.class_type.instance_type)

data = [("S", float64[:]),
        ("source_type", int64[:]),
        ("uncollided", int64),
        ("moving", int64),
        ("M", int64),
        ("x0", float64),
        ("t", float64),
        ("xL", float64),
        ("xR", float64),
        ("argument", float64[:]),
        ("source_vector", float64[:]),
        ("temp", float64[:]),
        ("abxx", float64),
        ("xx", float64),
        ("ix", int64),
        ("xs_quad", float64[:]),
        ("ws_quad", float64[:]),
        ("mag", float64),
        ("term1", float64),
        ("term2", float64)
        
        ]
###############################################################################
@jitclass(data)
class source_class(object):
    def __init__(self, build):
        self.S = np.zeros(build.M+1).transpose()
        self.source_type = build.source_type
        self.uncollided = build.uncollided
        self.x0 = build.x0
        self.M = build.M
        self.xs_quad = build.xs_quad
        self.ws_quad = build.ws_quad
        self.moving = build.moving
        
    def integrate_quad(self, t, a, b, j, func):
        argument = (b-a)/2 * self.xs_quad + (a+b)/2
        self.S[j] = (b-a)/2 * np.sum(self.ws_quad * func(argument, t) * normPn(j, argument, a, b))
        
    def square_IC_uncollidided_solution(self, xs, t):
        temp = xs*0
        for ix in range(xs.size):
            xx = xs[ix]
            # abxx = abs(xx)
            if (t <= self.x0) and (xx >= -self.x0 + t) and (xx <= self.x0 - t):
                temp[ix] = math.exp(-t)
            elif t > self.x0  and (-t + self.x0 <=  xx) and (t - self.x0 >= xx):
                temp[ix] = math.exp(-t) * self.x0 / (t + 1e-12)
            elif (xx < t + self.x0) and (xx > -t - self.x0):
                if (self.x0 - xx >= t) and (self.x0 + xx <= t):
                    temp[ix] = math.exp(-t)*(t + xx + self.x0)/(2.0 * t + 1e-12)
                elif (self.x0 - xx <= t) and (self.x0 + xx >= t):
                    temp[ix] = math.exp(-t)*(t - xx + self.x0)/(2.0 * t + 1e-12)
                # else:
                #     temp[ix] = 0.0
            # else: 
            #     temp[ix] = 0.0
        return temp
    def gaussian_IC_uncollided_solution(self, xs, t):
        temp = xs*0
        sqrtpi = math.sqrt(math.pi)
        for ix in range(xs.size):
            xx = xs[ix]
            temp[ix] = math.exp(-t) * sqrtpi * (math.erf(2*t-2*xx) + math.erf(2*t+2*xx))/(8.0 * t + 1e-12)
        return temp
        
    def plane_IC_uncollided_solution(self, xs, t):
        temp = xs*0
        for ix in range(xs.size):
            if (-t <= xs[ix] <= t):
                temp[ix] = math.exp(-t)/(2*t+1e-18)
        return temp
    
    def plane_IC_uncollided_solution_integrated(self, t, xL, xR):
        self.S[0] = (math.exp(-t)/(2*t+self.x0)*math.sqrt(xR-xL))
    
    def make_source(self, t, xL, xR):
        if self.uncollided == True:
            if self.source_type[0] == 1:
                if self.moving == True:
                    self.plane_IC_uncollided_solution_integrated(t, xL, xR)
                elif self.moving == False:
                    for j in range(self.M+1):
                        self.integrate_quad(t, xL, xR, j, self.plane_IC_uncollided_solution)
                        # if xL >= -t and xR <= t:
                        #     if self.S[j] > 1e-14:
                                # print(self.S[j], math.exp(-t)/2/t*math.sqrt(xR-xL),t)
                                   
            elif self.source_type[1] == 1:
                for j in range(self.M+1):
                    self.integrate_quad(t, xL, xR, j, self.square_IC_uncollidided_solution)
                    
            elif self.source_type[3] == 1:
                for j in range(self.M+1):
                    self.integrate_quad(t, xL, xR, j, self.gaussian_IC_uncollided_solution)
                
                    
    def uncollided_solution(self, xs, t):
        if self.uncollided == True:
            if self.source_type[0] == 1:
                return self.plane_IC_uncollided_solution(xs, t)
            elif self.source_type[1] == 1:
                return self.square_IC_uncollidided_solution(xs, t)
            elif self.source_type[3] == 1:
                return self.gaussian_IC_uncollided_solution(xs, t)
        else:
            return xs*0
        
        
    
                

            