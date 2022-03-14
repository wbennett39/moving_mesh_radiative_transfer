#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 11:25:35 2022

@author: bennett
"""
import numpy as np

from .build_problem import build
from .matrices import G_L
from .sources import source_class
from .phi_class import scalar_flux
from .uncollided_solutions import uncollided_solution
from .numerical_flux import LU_surf

from numba.experimental import jitclass
from numba import int64, float64, deferred_type, prange

build_type = deferred_type()
build_type.define(build.class_type.instance_type)
matrices_type = deferred_type()
matrices_type.define(G_L.class_type.instance_type)
num_flux_type = deferred_type()
num_flux_type.define(LU_surf.class_type.instance_type)
source_type = deferred_type()
source_type.define(source_class.class_type.instance_type)
flux_type = deferred_type()
flux_type.define(scalar_flux.class_type.instance_type)
uncollided_solution_type = deferred_type()
uncollided_solution_type.define(uncollided_solution.class_type.instance_type)


data = [('N_ang', int64), 
        ('N_space', int64),
        ('M', int64),
        ('source_type', int64[:]),
        ('t', float64),
        ('sigma_t', float64[:]),
        ('sigma_s', float64[:]),
        ('IC', float64[:,:,:]),
        ('mus', float64[:]),
        ('ws', float64[:]),
        ('x0', float64),
        ("xL", float64),
        ("xR", float64),
        ("dxL", float64),
        ("dxR", float64),
        ("L", float64[:,:]),
        ("G", float64[:,:]),
        ("P", float64[:]),
        ("S", float64[:]),
        ("LU", float64[:]),
        ("U", float64[:]),
        ("V_new", float64[:,:,:]),
        ("V", float64[:,:,:]),
        ("V_old", float64[:,:,:]),
        
        
        ]
##############################################################################
@jitclass(data)
class rhs_class():
    def __init__(self, build):
        self.N_ang = build.N_ang 
        self.N_space = build.N_space
        self.M = build.M
        self.mus = build.mus
        self.ws = build.ws
        self.source_type = build.source_type
    # def __call__(self,t, V, mesh, matrices, num_flux, source, flux):
    def call(self,t, V, mesh, matrices, num_flux, source, uncollided_sol, flux):
        
        V_new = V.copy().reshape((self.N_ang, self.N_space, self.M+1))
        V_old = V_new.copy()
        
        mesh.move(t)

        for space in prange(self.N_space):
            
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
            if self.source_type[4] != 1:
                source.make_source(t, xL, xR, uncollided_sol)
            S = source.S
            for angle in range(self.N_ang):
                mul = self.mus[angle]
                if self.source_type[4] == 1:
                    source.make_source_not_isotropic(t, mul, xL, xR)
                num_flux.make_LU(t, mesh, V_old[angle,:,:], space, mul)
                LU = num_flux.LU
                # LU = np.zeros(self.M+1).transpose()
                # LU = LU_surf_func(V_old[angle,:,:],space,self.N_space,mul,self.M,xL,xR,dxL,dxR)
                U = np.zeros(self.M+1).transpose()
                U[:] = V_old[angle,space,:]
                RHS = np.dot(G,U) + -LU + mul*np.dot(L,U) - U + P + S/2
                V_new[angle,space,:] = RHS
        return V_new.reshape(self.N_ang*self.N_space*(self.M+1))
            
       
