#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 11:23:59 2022

@author: bennett
"""
import numpy as np
from build_problem import build
import math
###############################################################################

class G_L(build):
    def __init__(self, build):
        self.M = build.M
        self.L = np.zeros((self.M+1, self.M+1))
        self.G = np.zeros((self.M+1, self.M+1))
        for i in range(0,self.M+1):
            for j in range(0,self.M+1):
                if i > j and (i+j) % 2 !=0: 
                    self.L[i,j] = 2 * math.sqrt(2*i +1 )* math.sqrt(2*j+1)
                else:
                    self.L[i,j] = 0
    def __call__(self, xL, xR, dxL, dxR, mat):
        h = xR - xL
        ih = 1/h
        b = dxR
        a = dxL
        if mat == "L":
            return self.L/(xR-xL)
        if mat == "G":
            for i in range(0,self.M+1):
                for j in range(0,self.M+1):
                    if i==j:
                        self.G[i,j] = -0.5*(2*i+1)*ih*(b-a)
                    elif i>j:
                        if (i+j)%2 ==0:
                            self.G[i,j] = -math.sqrt(2*j+1)*math.sqrt(2*i+1)*ih*(b-a)
                        else:
                            self.G[i,j] = -math.sqrt(2*j+1)*math.sqrt(2*i+1)*ih*(b+a)
            return self.G