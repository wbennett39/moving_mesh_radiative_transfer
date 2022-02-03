#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 11:25:35 2022

@author: bennett
"""
import numpy as np
from build_problem import build
from functions import LU_surf_func
###############################################################################

class rhs_class():
    def __init__(self, build):
        self.N_ang = build.N_ang 
        self.N_space = build.N_space
        self.M = build.M
        self.mus = build.mus
        self.ws = build.ws
    def __call__(self,t, V, mesh, matrices, num_flux, source, flux):
        
        V_new = V.copy().reshape((self.N_ang, self.N_space, self.M+1))
        V_old = V_new.copy()
        
        mesh.move(t)

        for space in range(self.N_space):
            
            xR = mesh.edges[space+1]
            xL = mesh.edges[space]
            dxR = mesh.Dedges[space+1]
            dxL = mesh.Dedges[space]
            matrices.make_L(xL, xR)
            matrices.make_G(xL, xR, dxL, dxR)
            L = matrices.L
            G = matrices.G
            flux.make_P(V_old[:,space,:])
            P = flux.P
            source.make_source(t, xL, xR)
            S = source.S
            for angle in range(self.N_ang):
                mul = self.mus[angle]
                num_flux.make_LU(mesh, V_old[angle,:,:], space, mul)
                # LU = num_flux.LU
                LU = np.zeros(self.M+1).transpose()
                LU = LU_surf_func(V_old[angle,:,:],space,self.N_space,mul,self.M,xL,xR,dxL,dxR)
                U = np.zeros(self.M+1).transpose()
                U[:] = V_old[angle,space,:]
                RHS = np.dot(G,U) + -LU + mul*np.dot(L,U) - U + P + S/2
                V_new[angle,space,:] = RHS
        return V_new.reshape(self.N_ang*self.N_space*(self.M+1))
            
       