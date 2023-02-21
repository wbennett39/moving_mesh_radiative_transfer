#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 18:17:08 2022

@author: bennett
"""
import numpy as np
from .functions import normPn, dx_normPn
from numba.experimental import jitclass
from numba import int64, float64, deferred_type
from .uncollided_solutions import uncollided_solution

uncollided_solution_type = deferred_type()
uncollided_solution_type.define(uncollided_solution.class_type.instance_type)

data = [('N_ang', int64), 
        ('t', float64),
        ('M', int64),
        ('ws', float64[:]),
        ('xs', float64[:]),
        ('u', float64[:,:,:]),
        ('edges', float64[:]),
        ('uncollided', int64),
        ('dx_e', float64[:]),
        ('psi_out', float64[:,:]),
        ('phi_out', float64[:]), 
        ('e_out', float64[:])
        ]
@jitclass(data)
class make_output:
    def __init__(self, t, N_ang, ws, xs, u, M, edges, uncollided):
        self.N_ang = N_ang
        self.ws = ws
        self.xs = xs
        self.u = u 
        self.M = M
        self.edges = edges
        self.uncollided = uncollided
        self.t = t
    def make_phi(self, uncollided_solution):
        output = self.xs*0
        psi = np.zeros((self.N_ang, self.xs.size))
        for ang in range(self.N_ang):
            for count in range(self.xs.size):
                idx = np.searchsorted(self.edges[:], self.xs[count])
                if (idx == 0):
                    idx = 1
                if (idx >= self.edges.size):
                    idx = self.edges.size - 1
                if self.edges[0] <= self.xs[count] <= self.edges[-1]:
                    for i in range(self.M+1):
                        psi[ang, count] += self.u[ang,idx-1,i] * normPn(i,self.xs[count:count+1],float(self.edges[idx-1]),float(self.edges[idx]))[0]
        output = np.sum(np.multiply(psi.transpose(), self.ws), axis = 1)
        if self.uncollided == True:
            uncol = uncollided_solution.uncollided_solution(self.xs, self.t)
            output += uncol
        self.psi_out = psi
        self.phi_out = output
        return output
    
    def make_e(self):
        e = np.zeros((self.xs.size))
        self.dx_e = np.zeros((self.xs.size))
        for count in range(self.xs.size):
            idx = np.searchsorted(self.edges[:], self.xs[count])
            if (idx == 0):
                idx = 1
            if (idx >= self.edges.size):
                idx = self.edges.size - 1
            for i in range(self.M+1):
                if self.edges[0] <= self.xs[count] <= self.edges[-1]:
                    e[count] += self.u[self.N_ang,idx-1,i] * normPn(i,self.xs[count:count+1],float(self.edges[idx-1]),float(self.edges[idx]))[0]
                    if self.M <=11:
                        self.dx_e[count] += self.u[self.N_ang,idx-1,i] * dx_normPn(i,self.xs[count:count+1],float(self.edges[idx-1]),float(self.edges[idx]))[0]
        self.e_out = e
        return e

    
