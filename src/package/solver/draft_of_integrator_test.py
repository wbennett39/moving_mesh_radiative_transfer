#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 11:57:16 2022

@author: bennett
"""
import numpy as np
from functions import normPn
import quadpy
import scipy.integrate as integrate
import math

M = 6
xs_quad = quadpy.c1.gauss_legendre(M+2).points
ws_quad = quadpy.c1.gauss_legendre(M+2).weights
t = 1/2

def integrate_quad(a, b, i, j):
    argument = (b-a)/2 * xs_quad + (a+b)/2
    result =  (b-a)/2 * np.sum(ws_quad * normPn(i, argument, a, b)* normPn(j, argument, a, b))
    print("j = ", j, "i = ", i, result)\
        
def integrate_quad2(a, b):
    # t = 1/2
    argument = (b-a)/2 * xs_quad + (a+b)/2
    result =  (b-a)/2 * np.sum(ws_quad * math.exp(-t)/2/t* np.heaviside(t - abs(argument),1) * normPn(0, argument, a, b))
    return result

def norm_func(x, i, j, a, b):
    return normPn(i, x, a, b)* normPn(j, x, a, b)

def test_func_1(x):
    # t = 1/2
    return math.exp(-t)/2/t*normPn(0, x, 0, 1)*np.heaviside(t - abs(x),1)
    
    
i = 0
# for j in range(M+1):
#     integrate_quad(-1, 1, i, j)
#     print("scipy", integrate.fixed_quad(norm_func, -1,1, args = (i,j,-1,1))[0])
    
res = integrate.fixed_quad(test_func_1, 0, 1)[0]
res2 = integrate_quad2(0,1)
print(res)
print(res2)
# print(math.exp(-t)/2/t*math.sqrt(1),"answer")
