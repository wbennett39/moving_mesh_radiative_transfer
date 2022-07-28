#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 07:24:05 2022

@author: William Bennett
"""
import numpy as np
from numba import int64, float64
from numba.experimental import jitclass
import math


###############################################################################
data = [('N_ang', int64), 
        ('N_space', int64),
        ('M', int64),
        ('tfinal', float64),
        ('mus', float64[:]),
        ('ws', float64[:]),
        ('x0', float64),
        ("moving", int64),
        ("move_type", int64[:]),
        ("edges", float64[:]),
        ("edges0", float64[:]),
        ("Dedges", float64[:]),
        ("N_space", int64),
        ('middlebin', int64),
        ('sidebin', int64),
        ('speed', float64),
        ('Dedges_const', float64[:]),
        ('source_type', int64[:])
        ]
###############################################################################

    
@jitclass(data)
class mesh_class(object):
    def __init__(self, N_space, x0, tfinal, moving, move_type, source_type, edge_v):
        self.tfinal = tfinal
        self.N_space = N_space
        self.x0 = x0
        self.moving = moving
        self.move_type = move_type
        self.edges = np.zeros(N_space+1)
        self.edges0 = np.zeros(N_space+1)
        self.Dedges = np.zeros(N_space+1)
        self.N_space = N_space
        self.speed = edge_v
        self.source_type = source_type

        self.initialize_mesh()
    
    def mesh_sharing_square_source(self):
        if self.N_space == 4:
            self.middlebin = 2
            self.sidebin = 1
        elif self.N_space == 8:
            self.middlebin = 4
            self.sidebin = 2
        elif self.N_space == 16:
            self.middlebin = 4
            self.sidebin = 6
    
    def square_source_init_func(self):
        self.middlebin = int(self.N_space/2)
        self.sidebin = int(self.middlebin/2)
        # self.mesh_sharing_square_source()
        dx = 1e-11
        left = np.linspace(-self.x0-dx, -self.x0, self.sidebin + 1)
        right = np.linspace(self.x0, self.x0 + dx, self.sidebin + 1)
        middle = np.linspace(-self.x0, self.x0, self.middlebin + 1)
        self.edges = np.concatenate((left[:-1], middle[:-1], right[:]))
        self.edges0 = self.edges
        # initialize derivatives
        self.Dedges[0:self.sidebin] = (self.edges[0:self.sidebin] + self.x0 )/(self.edges[-1] - self.x0)
        self.Dedges[self.sidebin:self.sidebin+self.middlebin] = 0
        self.Dedges[self.middlebin+self.sidebin + 1:] = (self.edges[self.middlebin+self.sidebin + 1:] - self.x0)/(self.edges[-1] - self.x0)
        self.Dedges = self.Dedges * self.speed
        
    # def square_IC_init_func(self):
    #     self.middlebin = 2
    #     self.sidebin = int(self.N_space-2)/2
    #     dx = 1e-9
    #     left = np.linspace(-self.x0-dx, -self.x0, self.sidebin+1)
    #     right = np.linspace(self.x0, self.x0+dx, self.sidebin+1)
    #     middle = np.linspace(-self.x0, self.x0, 3)
    #     self.edges = np.concatenate((left[:-1], middle[:], right[1:]))
    #     self.edges0 = self.edges
        
    #     self.Dedges[0:self.sidebin] = (self.edges[0:self.sidebin] + self.x0 )/(self.edges[-1] - self.x0)
    #     self.Dedges[self.sidebin] = 1
    #     self.Dedges[self.sidebin + 2] = -1
    #     self.Dedges[self.middlebin+self.sidebin + 1:] = (self.edges[self.middlebin+self.sidebin + 1:] - self.x0)/(self.edges[-1] - self.x0)
    #     print("initial edges")
    #     print(self.edges0)
    #     print("############")
    #     print("initial speed")
    #     print(self.Dedges)
   
        
    # def square_IC_pass_x0(self):
    #     self.Dedges[1:-1] = self.edges0[1:-1]/self.edges0[-2]
        
    # def square_IC_move_func(self, t):
    #     dx = 1e-8
    #     dangerzone = self.x0 - dx
    #     outofthewoods = self.x0 + dx
        
        
    #     if t < dangerzone:
    #         self.edges = self.edges0 +  self.Dedges*t
            
    #     elif t >= dangerzone and t <= outofthewoods:
    #         self.edges[0:self.sidebin] = self.edges0[0:self.sidebin] +  self.Dedges[0:self.sidebin]*t
    #         self.edges[self.sidebin+2:] = self.edges0[self.sidebin+2:] +  self.Dedges[self.sidebin+2:]*t
    #         self.Dedges[self.sidebin:self.sidebin+3] = 0
    #         temp_array = np.array([1,0,-1])
    #         self.edges[self.sidebin:self.sidebin+3] = self.edges0[self.sidebin:self.sidebin+3] + temp_array * (dangerzone)
    #         # self.edges0[self.sidebin:self.sidebin+3] = self.edges[self.sidebin:self.sidebin+3]
            
    #     elif t > outofthewoods and t <= 2 * self.x0:
    #         self.edges[0:self.sidebin] = self.edges0[0:self.sidebin] +  self.Dedges[0:self.sidebin]*t
    #         self.edges[self.sidebin+2:] = self.edges0[self.sidebin+2:] +  self.Dedges[self.sidebin+2:]*t
    #         self.Dedges[self.sidebin:self.sidebin+3] = np.array([-1,0,1])
    #         self.edges[self.sidebin:self.sidebin+3] = np.array([-dx,0, dx]) + self.Dedges[self.sidebin:self.sidebin+3] * (t - outofthewoods)
        
                
    def move(self, t):
        if self.moving == True:
             if self.move_type[4] != 1:
                self.edges = self.edges0 + self.Dedges*t
            
             elif self.move_type[4] == 1:
                if t > 1e-10:
                    sqrt_t = math.sqrt(t)
                else:
                    sqrt_t = math.sqrt(1e-10)
                assert not math.isnan(sqrt_t)
                self.Dedges = self.Dedges_const * 8.0 * 0.5 / sqrt_t
                self.edges = self.edges0 + self.Dedges_const * 8 * sqrt_t


            # elif self.move_type[3] == 1 and self.edges.size !=3:
            #     if t < 2*self.x0:
            #         self.square_source_move_func(t)
            #     else:
            #         print(self.edges, "edges")
            #         print(self.Dedges, "Dedges")
            #         self.square_source_move_func(t)
            #         self.Dedges[self.sidebin] = -0.999
            #         self.Dedges[self.sidebin+self.middlebin] = 0.999
                    
            #         #stopping protocol
            #         if abs(self.edges[self.sidebin] - self.edges[self.sidebin-1]) < 1e-3:
            #             self.Dedges[self.sidebin] = 0
            #             self.Dedges[self.sidebin+self.middlebin] = 0
            #             print(self.edges)
            #             print("stopping")
            #             print("##   ##   ##   ")
                        
                    
                    # self.edges[self.sidebin] = self.edges0[self.sidebin] + (t-self.x0)*self.Dedges[self.sidebin]
                    # self.edges[self.sidebin+self.middlebin] =  self.edges0[self.sidebin+self.middlebin] + (t-self.x0)*self.Dedges[self.sidebin+self.middlebin]
                    
                    # self.square_source_move_func(t)
                    
            
        
    def initialize_mesh(self):
        if self.moving == True:
            if self.move_type[0] == 1 or self.move_type[2] == 1 or self.edges.size == 3:
                self.edges = np.linspace(-self.x0, self.x0, self.N_space+1)
                self.edges0 = np.linspace(-self.x0, self.x0, self.N_space+1)
                self.Dedges = self.edges/self.edges[-1] * self.speed
                self.Dedges_const = self.Dedges
                
            elif self.move_type[1] == 1 and self.edges.size !=3: # square source function
                self.square_source_init_func()
                self.Dedges_const = self.Dedges
            
            # elif self.move_type[3] == 1 and self.edges.size!=3:
            #     self.square_source_init_func()

        else:
            # static meshes

            if self.move_type[4] != 1:
                self.edges = np.linspace((-self.tfinal - self.x0) * self.speed, (self.tfinal + self.x0) * self.speed, self.N_space+1)
            elif self.move_type[4] == 1:
                if self.source_type[5] == 1: # Gaussian source 
                    self.edges = np.linspace((-4 * math.sqrt(self.tfinal) * self.speed - self.x0) , (4 * math.sqrt(self.tfinal) * self.speed + self.x0)  , self.N_space+1)
                elif self.source_type[2] == 1 or self.source_type[1] == 1: # Square source, IC
                    if self.N_space == 4:
                        half = self.N_space/2
                        left_edge = -8 * math.sqrt(self.tfinal) * self.speed - self.x0
                        right_edge = -1*left_edge 
                        edges_left = np.linspace(left_edge, -self.x0, int(half/2 + 1))
                        edges_right = np.linspace(self.x0, right_edge, int(half/2 + 1))
                        middle = np.linspace(-self.x0, self.x0, int(half + 1))
                        self.edges = np.concatenate((edges_left[:-1], middle, edges_right[1:]))
                    
                    elif self.N_space > 4:
                        half = self.N_space/2
                        left_edge = -8 * math.sqrt(self.tfinal) * self.speed - self.x0
                        right_edge = -1*left_edge 
                        edges_left = np.linspace(left_edge, -self.x0, int(half/2 + 1))
                        edges_right = np.linspace(self.x0, right_edge, int(half/2 + 1))
                        middle_1 = np.linspace(-self.x0 + 8 * math.sqrt(self.tfinal) * self.speed, self.x0 - 6 * math.sqrt(self.tfinal) * self.speed, int(3))
                        middle_2 = np.linspace(-self.x0, -self.x0 + 8 * math.sqrt(self.tfinal) * self.speed, int((half-2)/2 + 1))
                        middle_3 = np.linspace(self.x0 - 8 * math.sqrt(self.tfinal) * self.speed, self.x0, int((half-2)/2 + 1))
                        self.edges = np.concatenate((edges_left[:-1], middle_2[:-1], middle_1[:-1], middle_3[:], edges_right[1:]))

                
