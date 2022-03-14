#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 18:17:08 2022

@author: bennett
"""
import numpy as np
from .functions import normPn

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
                for i in range(self.M+1):
                    psi[ang, count] += self.u[ang,idx-1,i] * normPn(i,self.xs[count:count+1],float(self.edges[idx-1]),float(self.edges[idx]))[0]
        output = np.sum(np.multiply(psi.transpose(), self.ws), axis = 1)
        if self.uncollided == True:
            output += uncollided_solution.uncollided_solution(self.xs, self.t)
        return output
