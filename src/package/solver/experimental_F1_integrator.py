#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 13:31:27 2022

@author: bennett
"""

import scipy.integrate as integrate
import numpy as np
import math 
import matplotlib.pyplot as plt
from scipy.special import expi
from numba import njit


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


    
    


    
    

# def uncollided_square_s(x, t, x0 ,t0):
#     taus = find_intervals(x, t, x0, t0)
#     tau_1 = taus[0]
#     tau_2 = taus[1]
#     tau_3 = taus[2]
#     tau_4 = taus[3]
#     t1 = f1(t, tau_2, x0) - f1(t, tau_1, x0)
#     t2 = f2(t, tau_3, x0, x) - f2(t, tau_2, x0, x)
#     t3 = f3(t, tau_4, x0, x) - f2(t, tau_3, x0, x)
#     return t1 + t2 + t3
    

def integrand(tau, t, x, x0):
    tp = t - tau
    abstp = abs(tp)

    if abstp <= abs(x) - x0:
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



# xs  = np.linspace(0.0,1.5, 100)
# phi = xs*0
# phi_u = xs*0
# tf = 0.005

# for i in range(len(xs)):
#     integral = integrate.nquad(integrand, [[0.0, tf]], args =  (tf, xs[i], 0.5))[0]
#     phi[i] = integral
#     phi_u[i] = uncollided_square_s2(xs[i], tf, 0.5, tf)
# plt.figure(11)
# plt.plot(xs, phi, "-o")
# plt.plot(xs, phi_u, "x")

