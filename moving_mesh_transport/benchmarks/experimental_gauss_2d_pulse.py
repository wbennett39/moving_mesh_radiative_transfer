#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 10:52:51 2022

@author: bennett
"""

import scipy.integrate as integrate
import math
import numpy as np
import matplotlib.pyplot as plt
from benchmark_functions import F_2D_gaussian_pulse

def uncollided_gauss_2D_integrand(s, rho, t, x0):
    
    if rho**2 + s**2 -2*rho*s > 0:
        eta = math.sqrt(rho**2 + s**2 - 2*rho*s)
    
        if abs(eta) < 1 and eta > 0:
        
            # res = s*0 * math.exp(-s**2/x0**2) / math.sqrt(1-eta**2) * math.exp(-t)/t/t
            res =  math.exp(-s**2/x0**2) / math.sqrt(1-eta**2) * math.exp(-t)/t/t
        else:
            res = 0
    else:
        res = 0
    return res

def opts0(*args, **kwargs):
       return {'limit':10000000}
   
""" uncollided """

def uncollided_gauss_2D(rho, t, x0):
    b = rho + t
    a = max(0.0, rho-t)
    res = integrate.nquad(uncollided_gauss_2D_integrand, [[a, b]], args = (rho, t, x0), opts = [opts0])[0]
    
    return res

""" collided """



def gaussian_pulse_2D_double_integral(s, rho, t, x0):
    eta = (rho-s)/t
    omega_a = 0.0
    omega_b = math.sqrt(1-eta**2)
    res = integrate.nquad(F_2D_gaussian_pulse, [[0, math.pi],[omega_a, omega_b]], args = (s, rho, t, x0), opts = [opts0, opts0])[0]
    return res


def collided_gauss_2D(rho, t, x0):
    b = rho + t
    a = max(0.0, rho-t)
    res = integrate.nquad(gaussian_pulse_2D_double_integral, [[a, b]], args = (rho, t, x0), opts = [opts0])[0]
    
    return res

###############################################################################
###############################################################################
###############################################################################


x0 = 0.5
tfinal = 1
pnts = 500

rhos = np.linspace(0, tfinal + 1/x0, pnts)
phi_u = rhos*0
phi_c = rhos*0

for ix in range(rhos.size):
    phi_u[ix] = uncollided_gauss_2D(rhos[ix], tfinal, x0)
    phi_c[ix] = collided_gauss_2D(rhos[ix], tfinal, x0)

def source(rho, x0):
    return np.exp(-rho**2/x0**2)

plt.figure(1)
plt.plot(rhos, phi_u, "-o", mfc = 'none')
plt.plot(rhos, phi_c + phi_u, "-s", mfc = 'none')


# plt.plot(rhos, source(rhos,x0))
plt.show()

# plt.figure(2)
# phi_l = rhos*0

# for ix in range(rhos.size):
#     phi_l[ix] = uncollided_gauss_2D_integrand(0, rhos[ix], tfinal, x0)
    

# plt.plot(rhos, phi_l)
# plt.show()