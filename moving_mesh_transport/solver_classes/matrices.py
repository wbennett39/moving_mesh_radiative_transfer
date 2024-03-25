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
import numba as nb

###############################################################################
build_type = deferred_type()
build_type.define(build.class_type.instance_type)
params_default = nb.typed.Dict.empty(key_type=nb.typeof('par_1'),value_type=nb.typeof(1))

data = [("M", int64),
        ('N', int64),
        ("L", float64[:,:]),
        ("L_const", float64[:,:]),
        ("G", float64[:,:]),
        ("xL", float64),
        ("xR", float64),
        ("dxL", float64),
        ("dxR", float64),
        ("mat", int64), 
        ('Mass_denom', int64[:,:]),
        ('J_denom', int64[:,:]),
        ('G_denom', int64[:,:]),
        ('L_denom', int64[:,:]),
        ('VV_denom', int64[:,:,:]),
        ('Mass_coeff_even', float64[:,:,:]),
        ('Mass_coeff_odd', float64[:,:,:]),
        ('J_coeff_even', float64[:,:,:]),
        ('J_coeff_odd', float64[:,:,:]),
        ('G_coeff', float64[:,:,:]),
        ('L_coeff', float64[:,:,:]),
        ('VV_coeff_even', float64[:,:,:,:]),
        ('VV_coeff_odd', float64[:,:,:,:]),
        ('Mass', float64[:,:]),
        ('J', float64[:,:]),
        ('VV', float64[:,:,:]),
        ('geometry', nb.typeof(params_default)),

        ]
@jitclass(data)
class G_L:
    def __init__(self, build, Mass_denom, J_denom, G_denom, L_denom, VV_denom, Mass_coeff_even, Mass_coeff_odd, 
                J_coeff_even, J_coeff_odd, G_coeff, L_coeff, VV_coeff_even, VV_coeff_odd):
        self.M = build.M
        self.L = np.zeros((self.M+1, self.M+1))
        self.L_const = np.zeros((self.M+1, self.M+1))
        self.G = np.zeros((self.M+1, self.M+1))
        self.Mass = np.zeros((self.M+1, self.M+1))
        self.J = np.zeros((self.M+1, self.M+1))
        self.VV = np.zeros((self.N+1, self.M+1, self.M+1))

        self.Mass_denom = Mass_denom[:]
        self.J_denom = J_denom[:]
        self.G_denom = G_denom[:]
        self.L_denom = L_denom[:]
        self.VV_denom = VV_denom[:]
        self.Mass_coeff_even = Mass_coeff_even[:]
        self.Mass_coeff_odd = Mass_coeff_odd[:]
        self.J_coeff_even = J_coeff_even[:]
        self.J_coeff_odd = J_coeff_odd[:]
        self.G_coeff = G_coeff
        self.L_coeff = L_coeff
        self.VV_coeff_even = VV_coeff_even
        self.VV_coeff_odd = VV_coeff_odd
        self.geometry = build.geometry
        self.N = build.Msigma

        for i in range(0,self.M+1):
            for j in range(0,self.M+1):
                if i > j and (i+j) % 2 !=0: 
                    self.L_const[i,j] = 2 * math.sqrt(2*i +1 )* math.sqrt(2*j+1)
                else:
                    self.L_const[i,j] = 0
    
    def make_L(self, xL, xR):
        if self.geometry['slab'] == True:
            self.make_L_slab(xL, xR)
        elif self.geometry['sphere'] == True:
            self.make_L_sphere(xL, xR)
    
    def make_G(self, xL, xR, dxL, dxR):
        if self.geometry['slab'] == True:
            self.make_G_slab(xL, xR, dxL, dxR)
        elif self.geometry['sphere'] == True:
            self.make_G_sphere(xL, xR, dxL, dxR)
    
    def make_all_matrices(self, xL, xR, dxL, dxR):
        self.make_L(xL, xR)
        self.make_G(xL, xR, dxL, dxR)
        if self.geometry['sphere'] == True:
            self.make_mass_sphere(xL, xR)
            self.make_J_sphere(xL, xR)
            # self.make_VV_sphere(xL, xR)

        


    def make_L_slab(self, xL, xR):
        self.L = self.L_const/(xR-xL)
        
    def make_G_slab(self, xL, xR, dxL, dxR):
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

    def make_mass_sphere(self, rL, rR):
        """This function builds the mass matrix for the spherical case"""
        rL2 = rL**2
        rR2 = rR**2
        rL3 = rL**3
        rR3 = rR**3
        rLrR = rL*rR
        pi = math.pi
        rttwo = math.sqrt(2)
        if self.M ==0:
            self.Mass[0,0] = (rL2 + rLrR + rR2)/3/pi
        elif self.M == 1:
            self.Mass[0,0] = (rL2 + rLrR + rR2)/3/pi
            self.Mass[1,0] = (rR2-rL2)/ 3/rttwo / pi 
            self.Mass[0,1] = (rR2-rL2)/ 3/rttwo / pi
            self.Mass[1,1] = 2*(2*rL2 + rLrR + 2*rR2)/15/pi
        # self.Mass[0,0] = 1/math.pi

    def make_J_sphere(self, rL, rR):
        """This function builds the J matrix for the spherical case"""
        pi = math.pi
        rL2 = rL**2
        rR2 = rR**2
        rttwo = math.sqrt(2)
        if self.M == 0:
            self.J[0,0] = 0.5 * (rR+rL) / pi
        elif self.M == 1:
            self.J[0,0] = 0.5 * (rR+rL) /  pi
            self.J[1,0] = (rR-rL) /3 /rttwo/pi
            self.J[0,1] = (rR-rL) /3 /rttwo/pi
            self.J[1,1] =  (rR+rL) / pi / 3
    
    


        # self.J[0,0] = math.log(rR/rL)/pi/(rR-rL)



    def make_G_sphere(self, rL, rR, rLp, rRp):
        """This function builds the G matrix for the spherical case"""
        rL2 = rL**2
        rR2 = rR**2
        rLrR = rL*rR
        pi = math.pi

        self.G[0,0] = - (rL2 + rLrR + rR2) * (rLp - rRp) / 6 / pi / (rL-rR)

    def make_L_sphere(self, rL, rR):
        """This function builds the L matrix for the spherical case"""
        pi = math.pi
        rL2 = rL**2
        rR2 = rR**2
        rLrR = rL*rR
        rtwo = math.sqrt(2)
        if self.M == 0:
            self.L[0,0] = 0

        elif self.M == 1:
            self.L[0,0] = 0
            self.L[0,1] = 0
            self.L[1,0] = -2 * rtwo * (rL2 + rLrR + rR2) / 3/ pi /(rL-rR)
            self.L[1,1] = 2*(rL + rR)/3/pi

        # self.L[0,0]  = 2*math.log(rR/rL)/pi/(rR-rL)
    
    def make_VV_sphere(self, rL, rR):
        """This function builds the VV matrix for the spherical case"""
        self.VV[0,0] = 0
        