#!/usr/bin/env python3# -*- coding: utf-8 -*-"""Created on Wed Jan 26 07:24:05 2022@author: William Bennett"""import numpy as npfrom numba import int64, float64from numba.experimental import jitclassimport math###############################################################################data = [('N_ang', int64),         ('N_space', int64),        ('M', int64),        ('tfinal', float64),        ('mus', float64[:]),        ('ws', float64[:]),        ('x0', float64),        ("moving", int64),        ("move_type", int64[:]),        ("edges", float64[:]),        ("edges0", float64[:]),        ("Dedges", float64[:]),        ("N_space", int64),        ('middlebin', int64),        ('sidebin', int64),        ('speed', float64)        ]###############################################################################    @jitclass(data)class mesh_class(object):    def __init__(self, N_space, x0, tfinal, moving, move_type):        self.tfinal = tfinal        self.N_space = N_space        self.x0 = x0        self.moving = moving        self.move_type = move_type        self.edges = np.zeros(N_space+1)        self.edges0 = np.zeros(N_space+1)        self.Dedges = np.zeros(N_space+1)        self.N_space = N_space                if self.move_type[2] == 1:            self.speed = 1/math.sqrt(3)        else:            self.speed = 1                    self.initialize_mesh()        def mesh_sharing_square_source(self):        if self.N_space == 4:            self.middlebin = 2            self.sidebin = 1        elif self.N_space == 8:            self.middlebin = 4            self.sidebin = 2        elif self.N_space == 16:            self.middlebin = 4            self.sidebin = 6        def square_source_init_func(self):        self.middlebin = int(self.N_space/2)        self.sidebin = int(self.middlebin/2)        # self.mesh_sharing_square_source()        dx = 1e-11        left = np.linspace(-self.x0-dx, -self.x0, self.sidebin + 1)        right = np.linspace(self.x0, self.x0 + dx, self.sidebin + 1)        middle = np.linspace(-self.x0, self.x0, self.middlebin + 1)        self.edges = np.concatenate((left[:-1], middle[:-1], right[:]))        self.edges0 = self.edges        # initialize derivatives        self.Dedges[0:self.sidebin] = (self.edges[0:self.sidebin] + self.x0 )/(self.edges[-1] - self.x0)        self.Dedges[self.sidebin:self.sidebin+self.middlebin] = 0        self.Dedges[self.middlebin+self.sidebin + 1:] = (self.edges[self.middlebin+self.sidebin + 1:] - self.x0)/(self.edges[-1] - self.x0)        self.Dedges = self.Dedges * self.speed    def square_source_move_func(self, t):        self.edges = self.edges0 +  self.Dedges*t        def move(self, t):        if self.moving == True:            if self.move_type[0] == 1 or self.edges.size == 3:                self.edges = self.edges0 + self.Dedges*t            elif self.move_type[1] == 1 and self.edges.size != 3:                self.square_source_move_func(t)                        def initialize_mesh(self):        if self.moving != False:            if self.move_type[0] == 1 or self.edges.size == 3:                self.edges = np.linspace(-self.x0, self.x0, self.N_space+1)                self.edges0 = np.linspace(-self.x0, self.x0, self.N_space+1)                self.Dedges = self.edges/self.edges[-1] * self.speed                            elif self.move_type[1] == 1 and self.edges.size !=3:                self.square_source_init_func()        else:            self.edges = np.linspace(-self.tfinal - self.x0, self.tfinal + self.x0, self.N_space+1)        