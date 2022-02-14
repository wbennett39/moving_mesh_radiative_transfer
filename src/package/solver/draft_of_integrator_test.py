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

M = 6
xs_quad = quadpy.c1.gauss_lobatto(M+2).points
ws_quad = quadpy.c1.gauss_lobatto(M+2).weights


def integrate_quad(a, b, i, j):
    argument = (b-a)/2 * xs_quad + (a+b)/2
    result =  (b-a)/2 * np.sum(ws_quad * normPn(i, argument, a, b)* normPn(j, argument, a, b))
    print("j = ", j, "i = ", i, result)
def norm_func(x, i, j, a, b):
    return normPn(i, x, a, b)* normPn(j, x, a, b)
    
    
i = 0
for j in range(M+1):
    integrate_quad(-1, 1, i, j)
    print("scipy", integrate.fixed_quad(norm_func, -1,1, args = (i,j,-1,1))[0])