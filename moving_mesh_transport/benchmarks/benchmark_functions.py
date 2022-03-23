#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 12 11:20:05 2022

@author: bennett
"""
import numpy as np
import math
import numba 
from scipy.special import expi
from numba import njit
from numba import cfunc,carray
from numba.types import intc, CPointer, float64
from scipy import LowLevelCallable


###############################################################################
def xi(u, eta):
    q = (1+eta)/(1-eta)
    if q <= 0:
        print(eta)
    zz = np.tan(u/2)
    
    # return (np.log(q) + complex(0, u))/(eta + complex(0, zz))
    return (np.log(q) + u*1j)/(eta + zz*1j)

# @cfunc("float64(float64, float64, float64, float64, float64)")
# @jit
def pyF1(u, s, tau, x, t):
    xp = x-s
    tp = t-tau
    if abs(xp) <= tp:
        eta = xp/tp
        eval_xi = xi(u, eta)
        complex_term = np.exp(tp*((1 - eta**2)*eval_xi/2.))*eval_xi**2
        return (1/np.cos(u/2.0))**2*complex_term.real * (4/math.pi) * (1 - eta**2) * math.exp(-tp)/2 # this seems wrong
    else:
        return 0.0

def get_intervals(x, x0, t, t0):
    intervals = np.zeros(4)
    intervals[0] = 0.0 
    intervals[3] = min(t, t0, t - abs(x) + x0)
    intervals[2] = min(intervals[3], t + abs(x) - x0)
    intervals[1] = min(intervals[3], t - abs(x) - x0)
    for i in range(4):
        if intervals[i] < 0:
            intervals[i] = 0
    return intervals    
    
# @cfunc("float64(float64, float64, float64, float64)")
# @jit
def pyF(s, tau, t, x):
    xp = x-s
    tp = t - tau
    if tp != 0 and abs(xp) <= tp:     
        return math.exp(-tp)/2/tp
    else:
        return 0
    
def f1(t, tau, x0):
    return -x0 * expi(tau-t)
    
def f2(t, tau, x0, x):
    if tau != t:
        return 0.5*((-x0 + abs(x)) * expi(tau-t) + math.exp(tau - t))
    else:
        return 0.5 * math.exp(0.0)

def f3(t, tau, x0, x):
    return math.exp(tau-t)

def uncollided_square_s2(x, t, x0, t0):
    tau_1 = 0.0
    end = min(t0, t - abs(x) + x0)
    if end <= 0.0:
        return 0.0
    tau_2 = min(end, t - x0 - abs(x))
    if tau_2 < 0.0:
        tau_2 = 0.0
    tau_3 = min(end, t - x0 + abs(x))
    if tau_3 < 0.0:
        tau_3 = 0.0
    tau_4 = end
    
    t1 = f1(t, tau_2, x0) - f1(t, tau_1, x0)
    t2 = f2(t, tau_3, x0, x) - f2(t, tau_2, x0, x)
    t3 = f3(t, tau_4, x0, x) - f3(t, tau_3, x0, x)
    
    return t1 + t2 + t3

def F1_integrand(tau, t, x, x0):
    tp = t - tau
    abstp = abs(tp)

    if abstp <= abs(x) - x0 or abstp == 0:
        return 0.0
    else:
        if abstp == 0 or t < tau:
            return 0.0
        else:
            if abstp >= abs(x) + x0:  
                return   x0 * math.exp(-tp)/(tp)
            
            elif abstp < x0 + abs(x) and abstp >= x0 - abs(x):
               return (x0 - abs(x) + abstp) * math.exp(-tp)/(2*tp)
           
            # elif abstp < x0 + abs(x) and abstp >= x0 + x:
            #     return  (abstp - abs(x) + x0) * math.exp(-tp)/(2*tp)*0
            
            elif abstp < x0 + abs(x) and abstp < x0 - abs(x):
                return  abstp * math.exp(-tp)/(tp)
            
            else:
                print("exception", tau, t, x)
                return 0.0



############# low level callable functions ####################################
def jit_F1(integrand_function):
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
    
@njit 
def heaviside(arg):
    if arg >= 0.0:
        return 1.0
    else:
        return 0.0

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
        if abs(xi.real) < 1e-15:
            xi = 0.0 + xi.imag
        if abs(xi.imag) < 1e-15:
            xi = xi.real + 0.0*1j
        
        complex_term = np.exp(tp*((1 - eta**2)*xi/2.))*xi**2

        res = (1/np.cos(u/2.0))**2*complex_term.real * (1/math.pi/8) * (1 - eta**2) * math.exp(-tp) * source(s, source_type)
    
        return res
    
    else:
        return 0.0
    
@jit_F1
def F1_spacefirst(args):
    """ The integrand for the triple integral args = (u, s, tau, x, t)
    """
    u = args[0]
    s = args[2]
    tau = args[1]
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
        # if abs(xi.real) < 1e-15:
        #     xi = 0.0 + xi.imag
        # if abs(xi.imag) < 1e-15:
        #     xi = xi.real + 0.0*1j
        
        complex_term = np.exp(tp*((1 - eta**2)*xi/2.))*xi**2

        res = (1/np.cos(u/2.0))**2*complex_term.real * (1/math.pi/8) * (1 - eta**2) * math.exp(-tp) * source(s, source_type)
    
        return res
    
    else:
        return 0.0
    
@jit_F1
def F1_c2(args):
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
    eta = xp/(tp + 1e-10)
    
    heaviside_arg = x - (abs(x) - t - tau)
    
    ## find xi ##
    q = (1+eta)/(1-eta)
    zz = np.tan(u/2)
    xi = (np.log(q) + u*1j)/(eta + zz*1j)
    
    complex_term = np.exp(tp*((1 - eta**2)*xi/2.))*xi**2
    
    res = (1/np.cos(u/2.0))**2*complex_term.real * (tp/4/math.pi) * (1 - eta**2) * math.exp(-tp)/2/tp * source(s, source_type)*heaviside(heaviside_arg)
    return res

@jit_F1
def F1_c3(args):
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
    eta = xp/(tp + 1e-10)
    
    heaviside_arg_1 = x - (t-tau-s)
    heaviside_arg_2 = x - (s + (t-tau))
    
    ## find xi ##
    q = (1+eta)/(1-eta)
    zz = np.tan(u/2)
    xi = (np.log(q) + u*1j)/(eta + zz*1j)
    
    complex_term = np.exp(tp*((1 - eta**2)*xi/2.))*xi**2
    
    res = (1/np.cos(u/2.0))**2*complex_term.real * (tp/4/math.pi) * (1 - eta**2) * math.exp(-tp)/2/tp * source(s, source_type)*heaviside(heaviside_arg_1)*heaviside(heaviside_arg_2)
    return res

@jit_F1
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
    
@jit_F1
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
