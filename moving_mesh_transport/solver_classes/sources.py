#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 14:38:34 2022

@author: bennett
"""
import numpy as np
import math
from numba import float64, int64, deferred_type
from numba.experimental import jitclass

from .build_problem import build
from .functions import normPn
from .functions import numba_expi as expi
from .uncollided_solutions import uncollided_solution
from scipy.special import expi as expi2
###############################################################################
build_type = deferred_type()
build_type.define(build.class_type.instance_type)
uncollided_solution_type = deferred_type()
uncollided_solution_type.define(uncollided_solution.class_type.instance_type)

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
        ("term2", float64),
        ("tfinal", float64),
        ("t0", float64),
        ("t1", float64),
        ("t2", float64), 
        ("t3", float64),
        ("tau", float64)
        
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
        self.tfinal = build.tfinal
        self.t0 = self.tfinal
    
    def integrate_quad(self, t, a, b, j, func):
        argument = (b-a)/2 * self.xs_quad + (a+b)/2
        self.S[j] = (b-a)/2 * np.sum(self.ws_quad * func(argument, t) * normPn(j, argument, a, b))
        
    def integrate_quad_not_isotropic(self, t, a, b, j, mu, func):
        argument = (b-a)/2 * self.xs_quad + (a+b)/2
        self.S[j] = (b-a)/2 * np.sum(self.ws_quad * func(argument, t, mu) * normPn(j, argument, a, b))
    
    def MMS_source(self, xs, t, mu):
        temp = xs*0
        for ix in range(xs.size):
            if -t - self.x0 <= xs[ix] <= t + self.x0:
                # temp[ix] = - math.exp(-xs[ix]*xs[ix]/2)*(1 + (1+t)*xs[ix]*mu)/((1+t)**2)/2
                temp[ix] = -0.5*(1 + (1 + t)*xs[ix]*mu)/(math.exp(xs[ix]**2/2.)*(1 + t)**2)
        return temp*2.0
    
    def square_source(self, xs, t):
        temp = xs*0
        for ix in range(xs.size):
            if abs(xs[ix]) <= self.x0 and t <= self.t0:
                temp[ix] = 1.0
        return temp
            
    def gaussian_source(self, xs, t):
        temp = xs*0
        for ix in range(xs.size):
            x = xs[ix]
            if t <= self.t0:
                temp[ix] = math.exp(-4.0*x*x)
        return temp
        
        
    def make_source(self, t, xL, xR, uncollided_solution):
        if self.uncollided == True:
            if (self.source_type[0] == 1) and  (self.moving == True):
                    self.S[0] = uncollided_solution.plane_IC_uncollided_solution_integrated(t, xL, xR)
            else:
                for j in range(self.M+1):
                    self.integrate_quad(t, xL, xR, j, uncollided_solution.uncollided_solution)
                                  
        elif self.uncollided == False:
            if self.source_type[2] == 1:
                for j in range(self.M+1):
                    self.integrate_quad(t, xL, xR, j, self.square_source)
            elif self.source_type[5] == 1:
                for j in range(self.M+1):
                    self.integrate_quad(t, xL, xR, j, self.gaussian_source)

    def make_source_not_isotropic(self, t, mu, xL, xR):
            if self.source_type[4] ==1:
                for j in range(self.M+1):
                    self.integrate_quad_not_isotropic(t, xL, xR, j, mu, self.MMS_source)
                

        
        

            
