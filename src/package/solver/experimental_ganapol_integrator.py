#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 26 11:43:33 2022

@author: bennett
"""


import scipy.integrate as integrate 
import numpy as np
import math as math

x0 = 0.5

def Q(s, tau):
    return np.heaviside(x0 - np.abs(s), 1)

def xi(u, eta):
    q = (1+eta)/(1-eta)
    zz = np.tan(u/2)
    return (np.log(q) + complex(0, u))/(eta + complex(0, zz))
    
def F1(u, tnew, xnew, eta):
    eval_xi = xi(u, eta)
    # complex_term = (eval_xi**2 * math.exp(t * (1-eta**2) * eval_xi / 2))
    complex_term = np.exp(tnew*((1 - eta**2)*eval_xi/2.))*eval_xi**2
    return (1/np.cos(u/2.0))**2*complex_term.real

            
def ival(u, t, tau, x, s):
    tnew = t - tau
    xnew = x - s 
    eta = xnew/tnew
    if abs(xnew) < tnew and tnew != 0:
        t1 = math.exp(-tnew)/2/tnew
        t2 = (tnew/4/math.pi) * (1 - eta**2)
        term =  t1 * (1 + t2 * F1(u, tnew, xnew, eta))
    elif abs(xnew) == (tnew) and tnew != 0:
        # term = math.exp(-t)/2/t
        term = 0
    elif (tnew == 0):
        term = 0
    else:
        term = 0
    return term

def integrand(u, s, tau, x, t):
    res = Q(s, tau) * ival(u, t, tau, x, s)
    return res

def opts0(*args, **kwargs):
       return {'limit':5}

# result = integrate.nquad(integrand, [[0, math.pi], [-x0,x0], [0.0,1.0]], args = (0, 1.0), opts = [opts0, opts0, opts0])[0]




