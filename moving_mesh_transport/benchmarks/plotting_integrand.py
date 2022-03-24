#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 17:21:00 2022

@author: bennett
"""

import numpy as np
import math
import matplotlib.pyplot as plt

def uncol_integrand(x,s,t,tau,u):
    xp = x-s
    tp = t-tau
    
    
    if abs(xp) < (tp):
        eta = xp/(tp)
        ## find xi ##
        q = (1.0+eta)/(1.0-eta)
        zz = 1j*np.tan(u/2.0)
        xi = (np.log(q) + u*1j)/(eta + zz)
        if abs(xi.real) < 1e-12:
            xi = 0.0 + xi.imag
        if abs(xi.imag) < 1e-12:
            xi = xi.real + 0.0*1j
        # if (xi**2).real < 0:
        #     print("negative xi", x)
        complex_term = np.exp(tp*((1 - eta**2)*xi/2.))*xi**2
        
        return (1/np.cos(u/2.0))**2*complex_term.real * (math.exp(-tp)/8/math.pi) * (1 - eta**2)
    else:
        return 0.0
    

N = 50000
ss = np.linspace(0, 0.5, N)
x = 0
ts = np.linspace(0,1, N)
sol= ss*0
sol_t = ts*0
u = 3.1
t = 1
tau = 0.99
s = 0.0
x = 0

for i in range(ss.size):
    sol[i] = uncol_integrand(x,ss[i],t,tau,u)
    sol_t[i] = uncol_integrand(x,s,ts[i],tau,u)
    
plt.figure(1)
plt.plot(ss, sol, label = f"tau = {tau} ")
plt.legend()
plt.show()
# plt.figure(2)
# plt.plot(ts, sol_t)
# plt.show()