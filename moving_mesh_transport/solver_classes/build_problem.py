#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 07:24:05 2022

@author: William Bennett
"""
import numpy as np
from numba import int64, float64, jit, njit, deferred_type
from numba.experimental import jitclass
# from main import IC_func 
from .mesh import mesh_class
from .functions import normPn
from .mutables import IC_func

import yaml
from pathlib import Path

###############################################################################
mesh_class_type = deferred_type()
mesh_class_type.define(mesh_class.class_type.instance_type)
IC_func_type = deferred_type()
IC_func_type.define(IC_func.class_type.instance_type)

data = [('N_ang', int64), 
        ('N_space', int64),
        ('M', int64),
        ('tfinal', float64),
        ('sigma_t', float64[:]),
        ('sigma_s', float64[:]),
        ('IC', float64[:,:,:]),
        ('mus', float64[:]),
        ('ws', float64[:]),
        ('xs_quad', float64[:]),
        ('ws_quad', float64[:]),
        ('x0', float64),
        ('t0', float64),
        ("source_type", int64[:]),
        ("uncollided", int64),
        ("moving", int64),
        ("move_type", int64[:]),
        ("argument", float64[:]),
        ("temp", float64),
        ("t_quad", float64[:]),
        ("t_ws", float64[:]),
        ('scattering_ratio', float64),
        ('thermal_couple', int64),
        ('temp_function', int64[:]),
        ('e_init', float64),
        ('sigma', float64),
        ('particle_v', float64),
        ('edge_v', float64)
        ]
###############################################################################

@jitclass(data)
class build(object):
    def __init__(self, N_ang, N_space, M, tfinal, x0, t0, scattering_ratio, mus, ws, xs_quad, ws_quad, sigma_t, sigma_s, source_type,
                 uncollided, moving, move_type, t_quad, t_ws, thermal_couple, temp_function, e_initial, sigma, particle_v, edge_v):
        self.N_ang = N_ang
        self.N_space = N_space
        self.M = M
        self.tfinal = tfinal
        self.sigma_t = sigma_t
        self.sigma_s = sigma_s
        
        self.mus = mus
        self.ws = ws/np.sum(ws)
        self.xs_quad = xs_quad
        self.ws_quad = ws_quad
        self.x0 = x0
        self.source_type = source_type
        self.uncollided = uncollided 
        self.moving = moving
        self.move_type = move_type
        self.t_quad = t_quad
        self.t_ws = t_ws
        self.t0 = t0
        self.scattering_ratio = scattering_ratio
        self.thermal_couple = thermal_couple
        self.temp_function = temp_function
        self.sigma = sigma
        self.particle_v = particle_v
        self.edge_v = edge_v
        
        
        if self.thermal_couple == 0:
            self.IC = np.zeros((N_ang, N_space, M+1))
        elif self.thermal_couple == 1:
            self.IC = np.zeros((N_ang + 1, N_space, M+1))
            
       
        self.e_init = e_initial
        # self.e_initial = 1e-4
        
        
    def integrate_quad(self, a, b, ang, space, j, ic):
        argument = (b-a)/2*self.xs_quad + (a+b)/2
        self.IC[ang,space,j] = (b-a)/2 * np.sum(self.ws_quad * ic.function(argument) * normPn(j, argument, a, b))
        
        
    def integrate_e(self, a, b, space, j):
        argument = (b-a)/2*self.xs_quad + (a+b)/2
        self.IC[self.N_ang,space,j] = (b-a)/2 * np.sum(self.ws_quad * self.IC_e_func(argument) * normPn(j, argument, a, b))
    
    def IC_e_func(self,x):
        return np.ones(x.size) * self.e_init
                
    def make_IC(self):
        edges = mesh_class(self.N_space, self.x0, self.tfinal, self.moving, self.move_type, self.source_type, self.edge_v)
        edges_init = edges.edges
        
        if self.moving == False and self.source_type[0] == 1 and self.uncollided == False and self.N_space%2 == 0:
            # print(edges[self.N_space/2 + 1] - edges[self.N_space/2 - 1]) 
            right_edge_index = int(self.N_space/2 + 1)
            left_edge_index = int(self.N_space/2 - 1)
            self.x0 = edges_init[right_edge_index] - edges_init[left_edge_index]
            # temp = (edges_init[self.N_space/2 + 1] - edges_init[self.N_space/2 - 1]) 
            
        if self.thermal_couple == 1:
            
            for space in range(self.N_space):
                for j in range(self.M + 1):
                    self.integrate_e(edges_init[space], edges_init[space+1], space, j)
            
            
        ic = IC_func(self.source_type, self.uncollided, self.x0)
        for ang in range(self.N_ang):
            for space in range(self.N_space):
                for j in range(self.M + 1):
                    self.integrate_quad(edges_init[space], edges_init[space+1], ang, space, j, ic)
    
        
        
        
        
        
