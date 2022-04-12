#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 20:32:27 2022

@author: bennett
"""

import numpy as np
import math as math
import scipy.integrate as integrate 
import matplotlib.pyplot as plt
import numba 
from numba import njit
from numba import cfunc,carray
from numba.types import intc, CPointer, float64
from scipy import LowLevelCallable

def jit_F1(integrand_function):
    jitted_function = numba.jit(integrand_function, nopython=True)
    @cfunc(float64(intc, CPointer(float64)))
    def wrapped(n, xx):
        values = carray(xx,n)
        return jitted_function(values)
    return LowLevelCallable(wrapped.ctypes)

def opts0(*args, **kwargs):
       return {'limit':2}
@njit 
def point_collided(u, r, t):
    eta = r/t
    if eta >= 1:
        return 0.0
    else:
        
        first = math.exp(-t)/4/math.pi/r/t * math.log((1 + eta)/(1-eta))
        q = (1+eta)/(1-eta)
        zz = np.tan(u/2)
        xi = (np.log(q) + u*1j)/(eta + zz*1j)
        exp_arg = t/2 * (1- eta**2) * xi
        complex_term = (eta + 1j * zz)*xi**3 * np.exp(exp_arg)
        
        result = first + (1/2/math.pi) * math.exp(-t)/4/math.pi/r/t * (1/4) * (1- eta**2) * (1/math.cos(u/2))**2 * complex_term.real
        
    return result 

@njit
def eta_func_2d_gauss_cartesian(x, s, y, v):
    res = (x-s)**2 + (y-v)**2
    return math.sqrt(res)

@jit_F1
def F_2D_gaussian_pulse(args):
    u = args[0]
    omega = args[1]
    s = args[2]
    v = args[3]
    x = args[4]
    y = args[5]
    t = args[6]
    x0 = args[7]
    
    eta = eta = eta_func_2d_gauss_cartesian(x, s, y, v)
    
    if eta < 1:
        r_arg = t * math.sqrt(eta**2 + omega**2)
        return  2 * t * point_collided(u, r_arg, t) * math.exp(-(s**2 + v**2)/x0**2)
    else: 
        return 0 


######################uncollided_functions######################################


def integrand(s, v, x, y, t):
    res = 0.0
    eta = eta_func_2d_gauss_cartesian(x, s, y, v)
    if eta < 1:
        ft = math.exp(-t)/2/math.pi/t/t
        garg = -(s**2 + v**2)/0.5**2
        gt = math.exp(garg)
        res = ft / math.sqrt(1-eta**2) * gt  
    return res
    
def greens_solution_first_integral(v, x, y, t):
    sqrt_term = t**2 - v**2 + 2*v*y - y**2
    res = 0.0
    if sqrt_term >= 0:
        a = x - math.sqrt(sqrt_term)
        b = x + math.sqrt(sqrt_term)
        
        res = integrate.nquad(integrand, [[a, b]], args = (v, x, y, t), opts = [opts0])[0]
    return res

def greens_solution_second_integral(x, y, t):
    res = integrate.nquad(greens_solution_first_integral, [[-np.inf, np.inf]], args = (x, y, t), opts = [opts0])[0]
    return res
    
######################collided functions#######################################
def gaussian_pulse_2D_double_integral(s, v, x, y, t, x0):
    """ integrates over omega, and u
    """
    res = 0.0
    eta = eta = eta_func_2d_gauss_cartesian(x, s, y, v)
    omega_a = 0
    if eta < 1:
        omega_b = math.sqrt(1-eta**2)
        res = integrate.nquad(F_2D_gaussian_pulse, [[0, math.pi], [omega_a, omega_b]], args = (s, v, x, y, t, x0), opts = [opts0, opts0])[0]
    return res

def gaussian_pulse_2D_s_integral(v, x, y, t, x0):
    """ integrates over s and sets bounds based on the heaviside
    """
    sqrt_term = t**2 - v**2 + 2*v*y - y**2
    res = 0.0
    if sqrt_term >= 0:
        a = x - math.sqrt(sqrt_term)
        b = x + math.sqrt(sqrt_term)
        res = integrate.nquad(gaussian_pulse_2D_double_integral, [[a,b]], args = [v, x, y, t, x0])[0]
    return res
    
def collided_gauss_2D(x, y, t, x0):
    """ calls the integrator over s which, in turn calls the integrator over omega 
    and u. Basically we integrate the imaginary variable u, then z, then integrate
    over s with the heaviside absorbed into the interval, then finally over v
    """
    b = np.inf
    # b = rho + t
    # a = max(0.0, rho-t)
    a = -np.inf 
    res = integrate.nquad(gaussian_pulse_2D_s_integral, [[a, b]], args = (x, y, t, x0), opts = [opts0])[0]
    
    return res


################################################################################
def solve(tfinal, pnts):
    y = 0
    xs = np.linspace(0, tfinal + 1  , pnts)
    phi_u = xs*0
    phi_c = xs*0
    for ix in range(xs.size):
        print(ix)
        print("x = ", xs[ix])
        phi_u[ix] = greens_solution_second_integral(xs[ix], y, tfinal)
        phi_c[ix] = collided_gauss_2D(xs[ix], y, tfinal, 0.5)
    plt.plot(xs, phi_u, "-o")
    plt.plot(xs, phi_c, "-s")
    plt.plot(xs, phi_u + phi_c, "-^")
    plt.show()
    
solve(1, 25)




