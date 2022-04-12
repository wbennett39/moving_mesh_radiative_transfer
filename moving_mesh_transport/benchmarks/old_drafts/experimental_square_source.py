#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 07:05:26 2022

@author: bennett
"""

import numpy as np
import math as math
import scipy.integrate as integrate 
# from numba import njit, cfunc, jit
import matplotlib.pyplot as plt
from timeit import default_timer as timer
import h5py
from .benchmark_functions import F, F1, F1_spacefirst, F_gaussian_source, uncollided_square_source, pyF, get_intervals
from .make_benchmarks import write_to_file
from pathlib import Path


###############################################################################
    
def opts0(*args, **kwargs):
       return {'limit':10000000}
   
def opts1(*args, **kwargs):
       return {'limit':10000000}
def opts2(*args, **kwargs):
       return {'limit':10000000}
###############################################################################

def find_intervals(t, tau, x, x0):
    if (t - tau >= abs(x) + x0):
        a = -x0
        b = x0
    elif (t - tau < x0 + abs(x)) and (t - tau) >= x0 - x:
        a = abs(x) - t - tau
        b = x0
    elif t - tau < x0 + abs(x) and t - tau > abs(x) - x0:
        a = -(t-tau)
        b = (t-tau)
    elif t-tau < abs(x) - x0:
        a = 0 
        b = 0 
    return [a, b]
 
def find_intervals_time(t, x, s):
    a = 0 
    b = min(t, t - abs(x-s))
    if b < 0:
        b = 0
    return [a,b]
    

def double_integral(tau, x, t):
    x0 = 0.5
    ab = find_intervals(t, tau, x, x0)
    # ab[0] = - x0
    # ab[1] = x0
    collided_solution = integrate.nquad(F1, [[0, math.pi], [ab[0], ab[1]]], args =  (tau, x, t, 0), opts = [opts0, opts1])[0]
    
    return collided_solution

def time_integral(x, t):
    x0 = 0.5
    
    collided_solution = integrate.nquad(double_integral, [[0, t]], args = (x, t), opts = [opts2])[0]
    uncollided_solution = uncollided_square_source(x, t, x0, t)
    
    # collided_solution = integrate.nquad(F1, [[0, math.pi], [-x0, x0], [0,tfinal]], args =  (x, tfinal, 0), opts = [opts0, opts1, opts2])[0]

    
    return uncollided_solution + collided_solution

def double_integral_time(s, x, t):
    ab = find_intervals_time(t, x, s)
    collided_solution = integrate.nquad(F1_spacefirst, [[0, math.pi], [ab[0],ab[1]]], args =  (s, x, t, 0), opts = [opts0, opts1])[0]
    
    return collided_solution

def space_integral(x,t):
    x0 = 0.5
    
    collided_solution = integrate.nquad(double_integral_time, [[-x0, x0]], args = (x, t), opts = [opts2])[0]
    uncollided_solution = uncollided_square_source(x, t, x0, t)
    
    return uncollided_solution + collided_solution
    
def make_square_source(tfinal, npnts = [100]):
    x0 = 0.5
    print("t = ", tfinal)
    xs1 = np.linspace(0, x0, npnts[0])
    # xs3 = np.array([0.176,1.5])
    times = np.zeros(1)
    phi_sqs = xs1*0

    start = timer()
    for k in range(npnts[0]):
        phi_sqs[k] = time_integral(xs1[k], tfinal)
    times[0] = timer()  - start
    print("square source finished")
    
    plt.figure(-1)
    plt.ion()
    # plt.plot(xs1, phi_pl, "-.",label = "plane")
    # plt.plot(xs2, phi_sq, "--", label = "square IC")
    plt.plot(xs1, phi_sqs, "-", label = "square source")
    plt.show()
    
    write_to_file(xs1, phi_sqs, tfinal, 'square_source', npnts[0])
    
def make_square_source_spacefirst(tfinal, npnts = [100]):
    x0 = 0.5
    print("t = ", tfinal)
    xs1 = np.linspace(0, tfinal + x0, npnts[0])
    # xs3 = np.array([0.176,1.5])
    times = np.zeros(1)
    phi_sqs = xs1*0

    start = timer()
    for k in range(npnts[0]):
        phi_sqs[k] = space_integral(xs1[k], tfinal)
    times[0] = timer()  - start
    print("square source finished")
    
    plt.figure(-1)
    plt.ion()
    # plt.plot(xs1, phi_pl, "-.",label = "plane")
    # plt.plot(xs2, phi_sq, "--", label = "square IC")
    plt.plot(xs1, phi_sqs, "-", label = "square source")
    plt.show()
    
    write_to_file(xs1, phi_sqs, tfinal, 'square_source', npnts[0])