#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 16:40:27 2022

@author: bennett
"""

import numpy as np
import scipy.integrate as si
import numba
from numba import cfunc,carray,njit
from numba.types import intc, CPointer, float64
from scipy import LowLevelCallable
import math as math
from timeit import default_timer as timer

def jit_F1(integrand_function):
    jitted_function = numba.jit(integrand_function, nopython=True)
    @cfunc(float64(intc, CPointer(float64)))
    def wrapped(n, xx):
        values = carray(xx,n)
        return jitted_function(values)
    return LowLevelCallable(wrapped.ctypes)

def jit_F(integrand_function):
    jitted_function = numba.jit(integrand_function, nopython=True)
    @cfunc(float64(intc, CPointer(float64)))
    def wrapped(n, xx):
        values = carray(xx,n)
        return jitted_function(values)
    return LowLevelCallable(wrapped.ctypes)

def jit_F_gaussian_source(integrand_function):
    jitted_function = numba.jit(integrand_function, nopython=True)
    @cfunc(float64(intc, CPointer(float64)))
    def wrapped(n, xx):
        values = carray(xx,n)
        return jitted_function(values)
    return LowLevelCallable(wrapped.ctypes)


@njit
def source(s, source_type):
    if source_type == 0:     # square 
        return 1.0
    elif source_type == 1:   # gaussian 
        return np.exp(-4*s*s)
    

@jit_F1
def F1(args):
    """ The integrand for the triple integral args = (u, s, tau, x, t)
    """
    u = args[0]
    s = args[1]
    tau = args[2]
    x = args[3]
    t = args[4]
    source_type = args[5]
    
    ## define new variables  ##
    xp = x-s
    tp = t-tau
    if abs(xp) < tp and tp > 0:
        eta = xp/tp
        
        ## find xi ##
        q = (1+eta)/(1-eta)
        zz = np.tan(u/2)
        xi = (np.log(q) + u*1j)/(eta + zz*1j)
        
        complex_term = np.exp(tp*((1 - eta**2)*xi/2.))*xi**2
        
        return (1/np.cos(u/2.0))**2*complex_term.real * (tp/4/math.pi) * (1 - eta**2) * math.exp(-tp)/2/tp * source(s, source_type)
    else:
        return 0.0
@jit_F
def F(args):
    """ integrand for the double integral. ags = s, tau, t, x
    the  sqrt(pi)/8 is left out 
    """
    s = args[0]
    tau = args[1]
    t = args[2]
    x = args[3]
    source_type = args[4]
    ## define new variables
    xp = x - s
    tp = t - tau
    ###
    if 1 - abs(xp/tp) > 0.0 :  
        return math.exp(-tp)/2/tp * source(s, source_type)
    else:
        return 0.0
@jit_F_gaussian_source
def F_gaussian_source(args):
    tau = args[0]
    t = args[1]
    x = args[2]
    
    abx = abs(x)
    tp = t - tau
    
    if tp != 0:
        erf1 = math.erf(2*(tp - abx)) 
        erf2 = math.erf(2*(tp + abx))
        return math.exp(-tp)* (erf1 + erf2) / tp
    else:
        return 0.0
        
    
    
    

# parms = np.array([0,1])
# start = timer()
# integral_2 = si.nquad(F1, [[0, math.pi], [-0.5, 0.5], [0, 1.0]], parms)[0]
# end = timer()
# print(end-start)