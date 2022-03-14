#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 10:05:00 2022

@author: bennett
"""

from ..solver_classes.functions import *
from ..solver_classes.build_problem import build
from ..solver_classes.sources import source_class
from ..solver_classes.uncollided_solutions import uncollided_solution
# from functions import *
import scipy.integrate as integrate
from numpy.testing import assert_allclose as npassert
from scipy.special import expi
import numpy as np
import quadpy
import math
from numba import jit

def test_normPns():
    # normalized
    for i in range(10):
         res = integrate.fixed_quad(lambda x: normPn(i,x,0.,2.)**2,0,2,n=i+1)[0]
         npassert(res, 1.0)
         # orthogonal 
         for j in range(10):
             if j  !=  i:
                 res = integrate.fixed_quad(lambda x: normPn(i,x,0.,2.) * normPn(j,x,0.,2.),0,2,n= i + j + 2)[0]
                 npassert(abs(res), 0.0, rtol = 1e-11, atol = 1e-11)
def test_expi():
    xs_neg = np.linspace(-25, -1)
    xs_pos = np.linspace(1, 25)
    for i in range(xs_neg.size):
        xn = xs_neg[i]
        xp = xs_pos[i]
        res = (numba_expi(xn) - expi(xn)) + (numba_expi(xp) - expi(xp)) 
        npassert(abs(res), 0.0, rtol = 1e-11, atol = 1e-11)
        
def test_node_finder():
    edges = np.array([-1, 1])
    Ms = [1,2,3,4,5,6,7,8,9,10]
    for M in Ms:
        nodes = find_nodes(edges, M)
        true_nodes = quadpy.c1.gauss_legendre(M+1)
        npassert(nodes, true_nodes.points, rtol = 1e-11, atol = 1e-11)


def test_integrate_quad():
    t = 0.5
    
    edges = np.array([-t - 0.5, -0.5, 0.0, 0.5, t + 0.5])
    
    N_ang = 128
    N_space = 4
    M = 5
    tfinal = 5
    x0 = 0.5
    mus = quadpy.c1.gauss_lobatto(N_ang).points
    ws = quadpy.c1.gauss_lobatto(N_ang).weights
    xs_quad = quadpy.c1.gauss_legendre(M+2).points
    ws_quad = quadpy.c1.gauss_legendre(M+2).weights
    t_quad = quadpy.c1.gauss_legendre(100).points
    t_ws = quadpy.c1.gauss_legendre(100).weights
    sigma_t = np.ones(N_space)
    sigma_s = np.ones(N_space)
    move_type = np.array([1,0,0])
    source_type = np.array([0,1,0,0,0,0])
    uncollided = True
    moving = True
    
    initialize = build(N_ang, N_space, M, tfinal, x0, mus, ws, xs_quad, ws_quad, sigma_t, sigma_s, source_type,
                 uncollided, moving, move_type, t_quad, t_ws)
    
    source = source_class(initialize)
    uncollided = uncollided_solution(initialize)
    @jit
    def func(x,t):
        return np.cos(x)**2 + x**3 + 100*x
    for ix in range(N_space):
        for j in range(M+1):
            xL = edges[ix]
            xR = edges[ix+1]
            source.integrate_quad(t, xL, xR, j, func)
            normed_uncollided = lambda x: normPn(j, x, edges[ix], edges[ix+1]) * (np.cos(x)**2 + x**3 + 100*x)
            res_2 = integrate.fixed_quad(normed_uncollided, edges[ix], edges[ix+1], n = M+2)[0]
            npassert(source.S[j],res_2, rtol = 1e-11, atol = 1e-11)
    
    
    

