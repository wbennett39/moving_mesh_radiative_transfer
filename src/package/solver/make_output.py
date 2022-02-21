#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 18:17:08 2022

@author: bennett
"""
import numpy as np
import quadpy
import math

from functions import normPn
from build_problem import build

class make_output:
    def __init__(self, u, build, mesh, tfinal):
        # div = build.M+1
        div = 100
        self.M = build.M
        self.N_space = build.N_space
        self.N_ang = build.N_ang
        self.xs_list = np.zeros((self.N_space*div))
        self.phi_list = np.zeros((self.N_space*div))
        self.ws = build.ws
        
        mesh.move(tfinal)
        edges = mesh.edges
        for k in range(self.N_space):
            xL = edges[k]
            xR = edges[k+1]
            scheme = quadpy.c1.gauss_legendre(self.M+1)
            # xs = scheme.points
            xs = np.linspace(xL, xR, div)
            self.xs_list[k*div:(k+1)*div] = xL + (xs+1)/(2)*(xR-xL)
            sol_vec = np.zeros((self.N_ang,self.M+1,div))
            for ang in range(self.N_ang):
                for j in range(0,self.M+1):
                    sol_vec[ang,j,:] = math.sqrt(2*j+1)*normPn(j, xs, xL, xR)*u[ang,k,j]/math.sqrt(xR-xL)
            psi = np.zeros((self.N_ang,div)) 
            for ang2 in range(0,self.N_ang):
                psi[ang2,:] = np.sum(sol_vec[ang2,:,:],axis=0)
            phi = np.sum(np.multiply(psi.transpose(),self.ws),axis=1)*0
            # self.phi_list[k*div:(k+1)*div] += phi      
