#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 26 17:40:52 2022

@author: bennett
"""
import numpy as np
import math as math
import scipy.integrate as integrate 
from numba import njit, cfunc

import matplotlib.pyplot as plt


x0 = 0.5
tfinal = 1.0 

# @cfunc("float64[:](float64[:], float64)")
def Q(s, tau):
    return np.heaviside(x0 - np.abs(s), 1)

@cfunc("complex128(float64, float64)")
def xi(u, eta):
    q = (1+eta)/(1-eta)
    zz = np.tan(u/2)
    return (np.log(q) + complex(0, u))/(eta + complex(0, zz))

@cfunc("float64(float64, float64, float64)")
def F1(u, t, eta):
    eval_xi = xi(u, eta)
    # complex_term = (eval_xi**2 * math.exp(t * (1-eta**2) * eval_xi / 2))
    complex_term = np.exp(t*((1 - eta**2)*eval_xi/2.))*eval_xi**2
    return (1/np.cos(u/2.0))**2*complex_term.real * (t/4/math.pi) * (1 - eta**2) * math.exp(-t)/2/t

def ival(u, x, t):
    eta = x/t
    if abs(x) < t and t != 0:
        t1 = math.exp(-t)/2/t
        t2 = (t/4/math.pi) * (1 - eta**2)
        term =  t1 * (1 + t2 * F1(u, t, eta))
    elif abs(x) == t and t != 0:
        term = math.exp(-t)/2/t
        # term = 0
    elif t == 0:
        term = 0
    else:
        term = 0
    return term

def integrand(u, s, tau, x, t):
    res = Q(s, tau) * ival(u, x-s, t-tau)
    return res
def opts0(*args, **kwargs):
       return {'limit':100}

# def make_integral(x,t):
    




# result = integrate.nquad(integrand, [[0, math.pi], [-x0, x0], [0, tfinal]], args = (0.0, 1.0), opts = [opts0, opts0, opts0])[0]

# result = integrate.nquad(integrand, [[0, math.pi], [-x0, x0]], args = (0.0, 0.0, 1.0), opts = [opts0, opts0, opts0])[0]


npnts = 25
xlist = np.linspace(0, 1.5, npnts)
phi = np.zeros(npnts)
t = 1.0

for i in range(npnts):
    eta = xlist[i]/t
    # phi[i] = integrate.nquad(integrand, [[0, math.pi]], args = (0.0,0.0, xlist[i], 1.0), opts = [opts0, opts0, opts0])[0]
    FF1 = integrate.nquad(F1, [[0, math.pi]], args =  (1.0, xlist[i]/1.0), opts = [opts0, opts0, opts0])[0]
    t1 = math.exp(-t)/2/t
    # t2 = (t/4/math.pi) * (1 - eta**2)
    term =   (t1 +  FF1)
    phi[i] =  term
plt.figure(10)
plt.plot(xlist, phi, "-x")



