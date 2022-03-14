#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 11:23:59 2022

@author: bennett
"""
import numpy as np
from .build_problem import build
import math

from numba import int64, float64, deferred_type
from numba.experimental import jitclass
###############################################################################
build_type = deferred_type()
build_type.define(build.class_type.instance_type)

data = [("M", int64),
        ("L", float64[:,:]),
        ("L_const", float64[:,:]),
        ("G", float64[:,:]),
        ("xL", float64),
        ("xR", float64),
        ("dxL", float64),
        ("dxR", float64),
        ("mat", int64)
        ]
@jitclass(data)
class G_L:
    def __init__(self, build):
        self.M = build.M
        self.L = np.zeros((self.M+1, self.M+1))
        self.L_const = np.zeros((self.M+1, self.M+1))
        self.G = np.zeros((self.M+1, self.M+1))
        for i in range(0,self.M+1):
            for j in range(0,self.M+1):
                if i > j and (i+j) % 2 !=0: 
                    self.L_const[i,j] = 2 * math.sqrt(2*i +1 )* math.sqrt(2*j+1)
                else:
                    self.L_const[i,j] = 0
                    
    def make_L(self, xL, xR):
        self.L = self.L_const/(xR-xL)
        
    def make_G(self, xL, xR, dxL, dxR):
        h = xR - xL
        ih = 1/h
        b = dxR
        a = dxL
        for i in range(0,self.M+1):
            for j in range(0,self.M+1):
                if i==j:
                    self.G[i,j] = -0.5*(2*i+1)*ih*(b-a)
                elif i>j:
                    if (i+j)%2 ==0:
                        self.G[i,j] = -math.sqrt(2*j+1)*math.sqrt(2*i+1)*ih*(b-a)
                    else:
                        self.G[i,j] = -math.sqrt(2*j+1)*math.sqrt(2*i+1)*ih*(b+a)