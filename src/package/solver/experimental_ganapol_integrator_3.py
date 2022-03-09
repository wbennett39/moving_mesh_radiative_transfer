#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 26 20:51:32 2022

@author: bennett
"""

"""
[] add gaussian IC, source
[] clean up scripts/drafts

"""

import numpy as np
import math as math
import scipy.integrate as integrate 
from numba import njit, cfunc, jit
import matplotlib.pyplot as plt
from timeit import default_timer as timer
import h5py
from low_level_ganapol_integrator import F, F1
from experimental_F1_integrator import integrand


# @cfunc("complex128(float64, float64)")
# @jit
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
    if abs(xp) < tp:
        eta = xp/tp
        eval_xi = xi(u, eta)
        complex_term = np.exp(tp*((1 - eta**2)*eval_xi/2.))*eval_xi**2
        return (1/np.cos(u/2.0))**2*complex_term.real * (tp/4/math.pi) * (1 - eta**2) * math.exp(-tp)/2/tp
    else:
        return 0.0


# @cfunc("float64(float64, float64, float64, float64)")
# @jit
def pyF(s, tau, t, x):
    xp = x-s
    tp = t - tau
    if tp != 0 and abs(xp) <= tp:     
        return math.exp(-tp)/2/tp
    else:
        return 0
    
def opts0(*args, **kwargs):
       return {'limit':1000}
   
def opts1(*args, **kwargs):
       return {'limit':1000}
   
def opts2(*args, **kwargs):
       return {'limit':1000}

def do_ganapol(x, tfinal, x0):
    integral_1 = pyF(0.0, 0.0, tfinal, x)
    integral_2 = integrate.nquad(F1, [[0, math.pi]], args =  (0.0, 0.0, x, tfinal), opts = [opts0])[0]
    return integral_1 + integral_2

def do_square_ic(x, tfinal, x0):
    integral_1 = integrate.nquad(F, [[-x0, x0]], args =  ([0.0, tfinal, x]), opts = [opts0])[0]
    integral_2 = integrate.nquad(F1, [[0, math.pi], [-x0, x0]], args =  (0.0,x, tfinal), opts = [opts0, opts0, opts0])[0]
    return integral_1 + integral_2

def do_square_source(x, tfinal, x0):
    # integral_1 = 0*integrate.nquad(F, [[-x0, x0], [0, tfinal]], args =  (tfinal, x), opts = [opts0, opts1, opts2])[0]
    integral_1 = integrate.nquad(integrand, [[0.0, tfinal]], args =  (1.0, x, 0.5))[0]
    integral_2 = integrate.nquad(F1, [[0, math.pi], [-x0, x0], [0, tfinal]], args =  (x, tfinal), opts = [opts0, opts1, opts2])[0]
    # integral_1 = 0.0
    return integral_1 + integral_2

def plotter(tfinal, x0, npnts = [1000, 500, 500]):
    xs1 = np.linspace(0, tfinal, npnts[0])
    xs2 = np.linspace(0, tfinal + x0, npnts[1])
    xs3 = np.linspace(0, tfinal + x0, npnts[2])
    # xs3 = np.array([0.176,1.5])
    phi_pl = xs1*0
    phi_sq = xs2*0
    phi_sqs = xs3*0
    plt.figure(10)
    start = timer()
    for i in range(npnts[0]):
        phi_pl[i] = do_ganapol(xs1[i], tfinal, 0.0)
    print("plane finished")
    for j in range(npnts[1]):
        phi_sq[j] = do_square_ic(xs2[j], tfinal, x0)
    print("square IC finished")
    for k in range(npnts[2]):
        phi_sqs[k] = do_square_source(xs3[k], tfinal, x0)
    end = timer()
    print("square source finished")
    
    plt.plot(xs1, phi_pl, "-.",label = "plane")
    plt.plot(xs2, phi_sq, "--", label = "square IC")
    plt.plot(xs3, phi_sqs, ":", label = "square source")
    plt.xlabel("x")
    plt.ylabel("scalar flux")
    plt.legend()
    print("-   -   -   -   -   -   -   -   -")
    print("time elapsed", end-start)
    print(phi_sqs)
    print("-   -   -   -   -   -   -   -   -")
    
    f = h5py.File("benchmarks_3_4.hdf5", "a")
    
    plane_IC = f.create_group("plane_IC")
    square_IC = f.create_group("square_IC")
    square_source = f.create_group("square_source")
    
    plane_IC.create_dataset("t = 1", (2, npnts[0]), dtype = "f", data = (xs1, phi_pl))
    square_IC.create_dataset("t = 1", (2, npnts[1]), dtype = "f", data = (xs2, phi_sq))
    square_source.create_dataset("t = 1", (2, npnts[2]), dtype = "f", data = (xs3, phi_sqs))

    f.close()
    
    
    
# result = do_integral(0.0, 1.0, 0.5)
# result = do_heaviside(0.0, 1.0, 0.5)
plotter(1.0, 0.5)






