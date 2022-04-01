#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 12:01:21 2022

@author: bennett
"""
import math
import scipy.integrate as integrate 
from .benchmark_functions import uncollided_square_source, uncollided_square_IC, gaussian_source_integrand, uncollided_gauss_2D_integrand

def opts0(*args, **kwargs):
       return {'limit':10000000}
   
###############################################################################
class uncollided_class:
    
    def __init__(self, source_type, x0, t0):
        self.source_type = source_type
        self.x0 = x0
        self.t0 = t0
    
    def plane_IC(self, xs, t):
        """ uncollided scalar flux for 1D plane pulse 
        """
        temp = xs*0
        for ix in range(xs.size):
            if (-t <= xs[ix] <= t):
                temp[ix] = math.exp(-t)/(2*t+1e-12)
        return temp
        
    def square_IC(self, xs, t):
        """ uncollided scalar flux for 1D square pulse 
        """
        temp = xs*0
        for ix in range(xs.size):
            x = xs[ix]
            temp[ix] = uncollided_square_IC(x, t, self.x0)
        return temp
    
    def square_source(self, xs, t):
        """ uncollided scalar flux for 1D square source
        """
        temp = xs*0
        for ix in range(xs.size):
            x = xs[ix]
            temp[ix] = uncollided_square_source(x, t, self.x0, self.t0)
        return temp
    def gaussian_IC(self, xs, t):
        """ uncollided scalar flux for 1D Gaussian pulse with standard deviation x0
        """
        temp = xs*0
        sqrtpi = math.sqrt(math.pi)
        for ix in range(xs.size):
            xx = xs[ix]
            temp[ix] = math.exp(-t) * sqrtpi * (math.erf(2*t-2*xx) + math.erf(2*t+2*xx))/(8.0 * t + 1e-12)
        return temp 
    
    def gaussian_source(self, xs, t):
        """ uncollided scalar flux for 1D Gaussian source with standard deviation x0
        """
        self.t0 = min(self.t0, t)
        temp = xs*0
        sqrtpi = math.sqrt(math.pi)
        for ix in range(xs.size):
            x = xs[ix]
            result = integrate.nquad(gaussian_source_integrand, [[0, self.t0]], args =  (t, x))[0]
            temp[ix] = result
        return temp*sqrtpi/8.0  
    
    def gaussian_IC_2D(self, rhos, t):
        temp = rhos*0
        for ix in range(rhos.size):
            rho = rhos[ix]
            b = rho + t
            a = max(0.0, rho-t)
            temp[ix] = integrate.nquad(uncollided_gauss_2D_integrand, [[a, b]], args = (rho, t, self.x0), opts = [opts0])[0]
        
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
        elif self.source_type == 'gaussian_IC_2D':
            return self.gaussian_IC_2D(xs, t)