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
from benchmark_functions import F1_2D_gaussian_pulse, F2_2D_gaussian_pulse
from scipy import LowLevelCallable
import numba 
from numba import cfunc,carray
from numba.types import intc, CPointer, float64
from scipy.interpolate import interp1d as interp


def jit_F1(integrand_function):
    jitted_function = numba.jit(integrand_function, nopython=True)
    @cfunc(float64(intc, CPointer(float64)))
    def wrapped(n, xx):
        values = carray(xx,n)
        return jitted_function(values)
    return LowLevelCallable(wrapped.ctypes)

#############################################################################
def arcsec(x):
    return 1/math.cos(x)

def find_intervals_uncollided_theta(r, t, theta, rho):
    interval_final = [0,0]
    sqrt_term = -(r**2*(r - t - rho)*(r + t - rho)*rho**2*(r - t + rho)*(r + t + rho)*math.sin(theta)**2)
    if (sqrt_term > 0):
        t1 = -arcsec((2*r**2*rho**2)/(r*rho*(r**2 - t**2 + rho**2)*math.cos(theta) - math.sqrt(sqrt_term)))
        t2 = -t1
        t3 = -arcsec((2*r**2*rho**2)/(r*rho*(r**2 - t**2 + rho**2)*math.cos(theta) + math.sqrt(-(r**2*(r - t - rho)*(r + t - rho)*rho**2*(r - t + rho)*(r + t + rho)*math.sin(theta)**2))))
        t4 = -t3
        interval = np.array([t1, t2, t3, t4])
        for count, val in enumerate(interval):
            if val < 0:
                interval[count] = 0 
        interval = np.sort(interval)
        interval_final = [interval[0], interval[-1]]
    return interval_final
        
        
    return interval 
def find_intervals_uncollided_s(r, t, theta, thetap):
    sqrt_term = -r**2 + 2*t**2 + r**2*math.cos(2*theta-2*thetap)
    a = 0.0
    b = 0.0
    if sqrt_term >=0:
        denominator = 2*(math.cos(thetap)**2 + math.sin(thetap)**2)
        t2 = 2 * r * (math.cos(theta) * math.cos(thetap) + math.sin(theta) * math.sin(thetap))
        a = (-math.sqrt(2) * math.sqrt(sqrt_term) + t2)/denominator
        b = (math.sqrt(2) * math.sqrt(sqrt_term) + t2)/denominator
        if a < 0:
            a = 0
        if b < 0:
            b = 0
    return [a,b]
        
@jit_F1
def uncollided_gauss_2D_integrand(args):
    
    thetap = args[1]
    s = args[0]
    rho = args[2]
    t = args[3]
    theta = args[4]
    x0 = args[5]
    
    x = rho * math.cos(theta)
    y = rho * math.sin(theta)
    q = s * math.cos(thetap)
    v = s * math.sin(thetap)
    new_r = math.sqrt((x-q)**2 + (y-v)**2)
    eta = new_r/t
    res = 0.0
    if eta < 1:
        res =   s * math.exp(-s**2/x0**2) / math.sqrt(1-eta**2) * math.exp(-t)/t/t/2/math.pi
    return res



def opts0(*args, **kwargs):
       return {'limit':50}


def opts1(*args, **kwargs):
       return {'limit':51}
   
""" uncollided """

def uncollided_gauss_2D_s(thetap, theta, rho, t, x0):
    interval = find_intervals_uncollided_s(rho, t, theta, thetap)
       
    res = integrate.nquad(uncollided_gauss_2D_integrand,  [interval], args = (thetap, rho, t, theta, x0), opts = [opts1])[0]
    
    return res

def uncollided_gauss_2D_theta(theta, rho, t, x0):
    
    res = integrate.nquad(uncollided_gauss_2D_s,  [[0, 2*math.pi]], args = (theta, rho, t, x0), opts = [opts1])[0]
    
    return res


""" collided """



def gaussian_pulse_2D_double_integral(s, thetap, rho, t, theta, x0):
    """ integrates over u, omega
    """
    x = rho * math.cos(theta)
    y = rho * math.sin(theta)
    q = s * math.cos(thetap)
    v = s * math.sin(thetap)
    new_r = math.sqrt((x-q)**2 + (y-v)**2)
    eta = new_r/t
    omega_a = 0.0
    res = 0.0
    
    if eta < 1:
        omega_b = math.sqrt(1-eta**2)
        rest_collided = integrate.nquad(F2_2D_gaussian_pulse, [[0, math.pi], [omega_a, omega_b]], args = (thetap, s, rho, theta, t,  x0), opts = [opts0, opts0])[0]
        first_collided = integrate.nquad(F1_2D_gaussian_pulse, [[omega_a, omega_b]], args = (thetap, s, rho, theta, t,  x0), opts = [opts0])[0]
        res = rest_collided + first_collided
    return res

def collided_gauss_2D_s(thetap, rho, t, x0):
    """ integrates over s
    """
    theta = 0
    # b = np.inf
    # b = rho + t
    # interval = [a, b]
    interval = find_intervals_uncollided_s(rho, t, theta, thetap)
    res = integrate.nquad(gaussian_pulse_2D_double_integral, [interval], args = (thetap, rho, t, theta, x0), opts = [opts0])[0]
    
    return res

def collided_gauss_2D_theta(rho, t, x0):
    """ integrates over thetap
    """

    res = integrate.nquad(collided_gauss_2D_s, [[0, math.pi*2]], args = (rho, t, x0), opts = [opts0])[0]
    
    return res

###############################################################################
###############################################################################
###############################################################################


x0 = 0.5
tfinal = 1
pnts = 100

rhos = np.linspace(0, tfinal + 1, pnts)
# rhos = np.linspace(0, tfinal, pnts)
xs = np.linspace(0, tfinal, pnts)
ys = np.linspace(0, tfinal, pnts)

phi_u = rhos*0
phi_c = rhos*0

phi_uxy = np.zeros((xs.size, ys.size))

# for ix in range(rhos.size):
#     phi_u[ix] = uncollided_gauss_2D(rhos[ix], tfinal, x0)
    # phi_c[ix] = collided_gauss_2D(rhos[ix], tfinal, x0)
    
    
for ix in range(rhos.size):
    print(ix)
    print(rhos[ix])
    phi_u[ix] = uncollided_gauss_2D_theta(math.pi, rhos[ix], tfinal, x0)
    phi_c[ix] = collided_gauss_2D_theta(rhos[ix], tfinal, x0)
    
    

plt.figure(3)
# plt.plot(rhos, phi_u, "-o", mfc = 'none')
# plt.plot(rhos, phi_c, "-^", mfc = 'none')
plt.plot(rhos, phi_c + phi_u, "-s", mfc = 'none')


# plt.plot(rhos, source(rhos,x0))
plt.show()


results = np.zeros((rhos.size, 2))
results[:, 0] = rhos
results[:, 1] = phi_u + phi_c


with open('2d_gaussian_IC_scalar_flux.txt', 'w') as f:
    for row in results:
        f.write("%s\n" % str(row))
        
 
"""Minwoo's code"""
phi_min = np.loadtxt('new_rho.dat', unpack = True)
rhos_min = np.loadtxt('new_x.dat', unpack = True)
plt.plot(rhos_min, phi_min, "--k")

minwoos = interp(rhos_min, phi_min, kind = "cubic")

RMS = np.sqrt(np.mean(((phi_c + phi_u) - minwoos(rhos))**2))
print("RMS = ", RMS)
