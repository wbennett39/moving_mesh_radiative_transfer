#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 12:14:04 2022

@author: bennett
"""
import numpy as np
import math
from numba import float64, int64, deferred_type
from numba.experimental import jitclass

from .mesh import mesh_class
from .functions import normPn
from .build_problem import build
###############################################################################
mesh_class_type = deferred_type()
mesh_class_type.define(mesh_class.class_type.instance_type)

data = [("M", int64),
        ("N_space", int64),
        ("source_type", int64[:]),
        ("edges", float64[:]),
        ("Dedges", float64[:]),
        ("i", int64),
        ("h", float64),
        ("u", float64[:,:]),
        ("space", int64),
        ("mul", float64),
        ("ws_quad", float64[:]),
        ("xs_quad", float64[:]),
        ("LU", float64[:]),
        ("v0", float64),
        ("v1", float64),
        ("v2", float64),
        ("v3", float64),
        ("t", float64),
        ("argument", float64[:]),
        ("h", float64),
        ("hp", float64),
        ("hm", float64),
        ("xL_minus", float64),
        ("xR_plus", float64),
        
        
        
        ]
build_type = deferred_type()
build_type.define(build.class_type.instance_type)
###############################################################################
@jitclass(data)
class LU_surf(object):
    def __init__(self, build):
        self.LU = np.zeros(build.M+1).transpose()
        self.M = build.M
        self.source_type = build.source_type
        self.ws_quad = build.ws_quad
        self.xs_quad = build.xs_quad
        self.N_space = build.N_space
        self.v0 = 0.0 
        self.v1 = 0.0
        self.v2 = 0.0
        self.v3 = 0.0
        self.xL_minus = 0.0
        self.xR_plus = 0.0

    def integrate_quad(self, t, a, b, j, side):
        argument = (b-a)/2 * self.xs_quad + (a+b)/2
        res = (b-a)/2 * np.sum(self.ws_quad * self.BC_func(argument, t) * normPn(j, argument, a, b))
        return res
        
    def B_LR_func(self, i, h):
        if i == 0:
            B_right = 1/h
            B_left = 1/h
        elif i>0:
            B_right = math.sqrt(2*i+1)/h
            if i%2 ==0:
                B_left = math.sqrt(2*i+1)/h
            else: 
                B_left = -math.sqrt(2*i+1)/h
        return B_left, B_right
    
    def BC_func(self, xs, t):
        temp = xs*0
        if self.source_type[4] == 1:
            temp = np.exp(-xs*xs/2)/(t + 1)/2.0
        return temp
    
    def make_h(self, space):
        
        xR = self.edges[space+1]
        xL = self.edges[space]
        
        dx = xR - xL
        
        self.h = math.sqrt(xR-xL)
        
        if space != self.N_space - 1:
            self.hp = math.sqrt(self.edges[space+2] - self.edges[space+1])
        else:
            self.hm = math.sqrt(dx)
        
        if space != 0:
            self.hm = math.sqrt(self.edges[space]-self.edges[space-1])
        else:
            self.hm = math.sqrt(dx)
            
    def extend_mesh(self, space):
        xR = self.edges[space+1]
        xL = self.edges[space]
        dx = xR - xL
            
        self.xL_minus = self.edges[0] - dx
        self.xR_plus = self.edges[-1] + dx
    
    def make_sol(self, space, u, t):
        for j in range(self.M+1):
                if space != 0:
                    self.v0 += self.B_LR_func(j, self.hm)[1]*(u[space-1,j])
                    
                elif space == 0 and self.source_type[4] == 1:
                    self.v0 += self.integrate_quad(t, self.xL_minus, self.edges[space], j, "l") * self.B_LR_func(j, self.h)[1] 
                    
                self.v1 += self.B_LR_func(j, self.h)[0]*(u[space, j])
                self.v2 += self.B_LR_func(j, self.h)[1]*(u[space, j])
                
                if space != self.N_space - 1:
                    self.v3 += self.B_LR_func(j, self.hp)[0]*(u[space+1,j])
                    
                elif space == self.N_space - 1 and self.source_type[4] == 1:
                    self.v3 += self.integrate_quad(t, self.edges[space+1], self.xR_plus, j, "r") * self.B_LR_func(j, self.h)[0] 
            
    
    def make_LU(self, t, mesh_class, u, space, mul):
        self.v0 = 0 
        self.v1 = 0
        self.v2 = 0
        self.v3 = 0
        self.edges = mesh_class.edges
        self.Dedges = mesh_class.Dedges
        
        leftspeed = mul - self.Dedges[space]
        rightspeed = mul - self.Dedges[space+1]
        
        self.make_h(space)
        self.extend_mesh(space)
        self.make_sol(space, u, t)

        if leftspeed >= 0: 
            psi_minus = self.v0
        elif leftspeed < 0: 
            psi_minus = self.v1
        if rightspeed >= 0: 
            psi_plus = self.v2
        elif rightspeed < 0:
            psi_plus = self.v3
        for i in range(0,self.M+1):
            B_left, B_right = self.B_LR_func(i, self.h)
            self.LU[i] = (B_right*rightspeed*psi_plus - B_left*leftspeed*psi_minus)
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            