#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 12:14:04 2022

@author: bennett
"""
import numpy as np
import math
from mesh import mesh_class
###############################################################################
class LU_surf:
    def __init__(self, M):
        self.LU = np.zeros(M+1).transpose()
        self.M = M
        # self.edges = mesh.edges
        # self.Dedges = mesh.Dedges
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
    def __call__(self, mesh_class, u, space, mul):
        self.v0 = 0 
        self.v1 = 0
        self.v2 = 0
        self.v3 = 0
        leftspeed = mul - mesh_class.Dedges[space]
        rightspeed = mul - mesh_class.Dedges[space+1]
        h = math.sqrt(mesh_class.edges[space+1]-mesh_class.edges[space])
        if space != mesh_class.N_space - 1:
            hp = math.sqrt(mesh_class.edges[space+2]-mesh_class.edges[space+1])
        if space != 0:
            hm = math.sqrt(mesh_class.edges[space]-mesh_class.edges[space-1])
        for j in range(self.M+1):
            if space != 0:
                self.v0 += self.B_LR_func(j, hm)[1]*(u[space-1,j])
            self.v1 += self.B_LR_func(j, h)[0]*(u[space, j])
            self.v2 += self.B_LR_func(j, h)[1]*(u[space, j])
            if space != mesh_class.N_space - 1:
                self.v3 += self.B_LR_func(j, hp)[0]*(u[space+1,j])
        if leftspeed >= 0: 
            psi_minus = self.v0
        elif leftspeed < 0: 
            psi_minus = self.v1
        if rightspeed >= 0: 
            psi_plus = self.v2
        elif rightspeed < 0:
            psi_plus = self.v3
        for i in range(0,self.M+1):
            B_left, B_right = self.B_LR_func(i, h)
            self.LU[i] = (B_right*rightspeed*psi_plus - B_left*leftspeed*psi_minus)
        return self.LU