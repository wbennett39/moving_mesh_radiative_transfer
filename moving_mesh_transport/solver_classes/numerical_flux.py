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
from numba import types, typed
import numba as nb

from .mesh import mesh_class
from .functions import normPn
from .build_problem import build
###############################################################################
mesh_class_type = deferred_type()
mesh_class_type.define(mesh_class.class_type.instance_type)
params_default = nb.typed.Dict.empty(key_type=nb.typeof('par_1'),value_type=nb.typeof(1))

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
        ('thermal_couple', nb.typeof(params_default)),
        ('moving', int64),
        ('e_init', float64),
        ('speed', float64),
        ('test_dimensional_rhs', int64),
        ('boundary_on', int64[:]), 
        ('boundary_source_strength', float64),
        ('boundary_source', int64),
        ('uncollided', int64),
        ('t0', float64),
        ('geometry', nb.typeof(params_default))
        ]
build_type = deferred_type()
build_type.define(build.class_type.instance_type)
###############################################################################
@jitclass(data)
class LU_surf(object):
    def __init__(self, build):
        self.test_dimensional_rhs = build.test_dimensional_rhs
        self.LU = np.zeros(build.M+1).transpose()
        self.M = build.M
        self.source_type = np.array(list(build.source_type), dtype = np.int64)
        self.ws_quad = build.ws_quad
        self.xs_quad = build.xs_quad
        self.N_space = build.N_space
        self.t0 = build.t0
        self.v0 = 0.0 
        self.v1 = 0.0
        self.v2 = 0.0
        self.v3 = 0.0
        self.xL_minus = 0.0
        self.xR_plus = 0.0
        self.thermal_couple = build.thermal_couple
        self.moving = build.moving
        self.e_init = build.e_init
        self.speed = 1.0
        self.boundary_source = build.boundary_source
        self.boundary_on = np.array(list(build.boundary_on), dtype = np.int64) # [left, right]
        self.boundary_source_strength = build.boundary_source_strength
        self.uncollided = build.uncollided
        self.geometry = build.geometry
        # if build.test_dimensional_rhs == True:
        #     self.speed = build.particle_v

 

    def integrate_quad(self, t, a, b, j, side):
        argument = (b-a)/2 * self.xs_quad + (a+b)/2
        res = (b-a)/2 * np.sum(self.ws_quad * self.BC_func(argument, t) * normPn(j, argument, a, b))
        return res
        
    def B_LR_func(self, i, h):
        if self.geometry['slab'] == True:
            if i == 0:
                B_right = 1/h
                B_left = 1/h
            elif i>0:
                B_right = math.sqrt(2*i+1)/h
                if i%2 ==0:
                    B_left = math.sqrt(2*i+1)/h
                else: 
                    B_left = -math.sqrt(2*i+1)/h
            
        elif self.geometry['sphere'] == True:
            edgeval = 1 / h / math.sqrt(math.pi)
            if i == 0:
                B_left = edgeval
                B_right = edgeval
            else:
                B_right = math.sqrt(2) * edgeval
                if i%2 == 0: 
                    B_left = B_right
                else:
                    B_left = -B_right
        return B_left, B_right
            

    
    def BC_func(self, xs, t):
        temp = xs*0
        if self.source_type[4] == 1:
            temp = np.exp(-xs*xs/2)/(t + 1)/2.0
        elif (self.thermal_couple == 1) and (self.moving == 1):
            temp = np.ones(xs.size) * self.e_init 
        elif self.boundary_source == True:
            temp = np.ones(xs.size) * self.boundary_source_strength

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
        
    def is_boundary_source_on(self, space, t):
        returnval = False
        if self.source_type[4] == 1:
            returnval = True
        elif self.thermal_couple == 1:
            returnval = True
        elif self.boundary_source == True:
            if self.uncollided == False:
                if t <= self.t0:
                    if space == 0:
                        if self.boundary_on[0] == 1:
                            returnval = True
                    elif space == self.N_space - 1:
                        if self.boundary_on[1] == 1:
                            returnval = True
        return returnval


    def make_sol(self, space, u, t, u_refl):
        for j in range(self.M+1):
                if space != 0:
                    self.v0 += self.B_LR_func(j, self.hm)[1]*(u[space-1,j])
                    
                elif space == 0 and self.is_boundary_source_on(space, t): # special MMS case
                    self.v0 += self.integrate_quad(t, self.xL_minus, self.edges[space], j, "l") * self.B_LR_func(j, self.h)[1] 
                
                elif space == 0 and self.geometry['sphere'] == True: #reflecing BC for sphere
                    self.v0 += self.B_LR_func(j, self.h)[0]*(u_refl[j])

                    
                self.v1 += self.B_LR_func(j, self.h)[0]*(u[space, j])
                self.v2 += self.B_LR_func(j, self.h)[1]*(u[space, j])

                
                if space != self.N_space - 1:
                    self.v3 += self.B_LR_func(j, self.hp)[0]*(u[space+1,j])
                    
                elif space == self.N_space - 1 and self.is_boundary_source_on(space, t):
                    self.v3 += self.integrate_quad(t, self.edges[space+1], self.xR_plus, j, "r") * self.B_LR_func(j, self.h)[0] 
            
    
    def make_LU(self, t, mesh_class, u, space, mul, u_refl):
        self.v0 = 0 
        self.v1 = 0
        self.v2 = 0
        self.v3 = 0
        self.edges = mesh_class.edges
        self.Dedges = mesh_class.Dedges

        xL = self.edges[space]
        xR = self.edges[space + 1]
        
        leftspeed = self.speed * mul -  self.Dedges[space] 
        rightspeed = self.speed * mul -  self.Dedges[space+1] 
        
        self.make_h(space)
        self.extend_mesh(space)
        self.make_sol(space, u, t, u_refl)



        if leftspeed >= 0: 
            psi_minus = self.v0
            # if space == 0:
            #     assert (psi_minus == 1 / self.h / math.sqrt(math.pi) * u[space, 0]) 
        elif leftspeed < 0: 
            psi_minus = self.v1
        if rightspeed >= 0: 
            psi_plus = self.v2
        elif rightspeed < 0:
            psi_plus = self.v3
        

        for i in range(0,self.M+1):
            B_left, B_right = self.B_LR_func(i, self.h)
            if self.geometry['slab'] == True:
                self.LU[i] = (B_right*rightspeed*psi_plus - B_left*leftspeed*psi_minus)
            elif self.geometry['sphere'] == True:
                self.LU[i] = (xR**2*B_right*rightspeed*psi_plus - xL**2*B_left*leftspeed*psi_minus)

                # if space == 0:
                #     if self.geometry['sphere'] == True:
                #         if mul > 0:
                #             LUanalytic = mul * math.sqrt(1/math.pi) * math.sqrt(1/(xR-xL)) * (xR**2* psi_plus - xL**2*psi_minus)
                #             if abs(LUanalytic- self.LU[0]) >= 1e-4:
                #                 print('error', abs(LUanalytic- self.LU[0]) )
                #                 assert(0)
                # elif space == self.N_space - 1:
                #     if rightspeed < 0:
                #         assert(abs(psi_plus)<=1e-10) 






 

            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            