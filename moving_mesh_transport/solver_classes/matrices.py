#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 11:23:59 2022

@author: bennett
"""
import numpy as np
from .build_problem import build
# from .chebyshev_matrix_builder import matrix_builder
import math
from .functions import sqrt_two_mass_func as rtf
from .functions import rttwo_mistake_undoer as rund
from .GMAT_sphere import GMatrix, MPRIME

from numba import int64, float64, deferred_type
from numba.experimental import jitclass
import numba as nb

###############################################################################
build_type = deferred_type()
build_type.define(build.class_type.instance_type)
params_default = nb.typed.Dict.empty(key_type=nb.typeof('par_1'),value_type=nb.typeof(1))

data = [("M", int64),
        ('MPRIME', float64[:, :]),
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
        ('J_coeff', float64[:,:,:]),
        ('G_coeff', float64[:,:,:]),
        ('L_coeff_even', float64[:,:,:]),
        ('L_coeff_odd', float64[:,:,:]),
        ('VV_coeff_even', float64[:,:,:,:]),
        ('VV_coeff_odd', float64[:,:,:,:]),
        ('Mass', float64[:,:]),
        ('J', float64[:,:]),
        ('VV', float64[:,:,:]),
        ('geometry', nb.typeof(params_default)),
        ('testing', int64)

        ]
@jitclass(data)
class G_L:
    def __init__(self, build, Mass_denom, J_denom, G_denom, L_denom, VV_denom, Mass_coeff_even, Mass_coeff_odd, 
                J_coeff, G_coeff, L_coeff_even, L_coeff_odd, VV_coeff_even, VV_coeff_odd):
        self.M = build.M
        self.L = np.zeros((self.M+1, self.M+1))
        self.L_const = np.zeros((self.M+1, self.M+1))
        self.G = np.zeros((self.M+1, self.M+1))
        self.Mass = np.zeros((self.M+1, self.M+1))
        self.J = np.zeros((self.M+1, self.M+1))
        self.VV = np.zeros((self.N+1, self.M+1, self.M+1))
        self.MPRIME = np.zeros((self.M+1, self.M+1))

        self.Mass_denom = Mass_denom[:]
        self.J_denom = J_denom[:]
        self.G_denom = G_denom[:]
        self.L_denom = L_denom[:]
        self.VV_denom = VV_denom[:]
        self.Mass_coeff_even = Mass_coeff_even[:]
        self.Mass_coeff_odd = Mass_coeff_odd[:]
        self.J_coeff = J_coeff[:]
        self.G_coeff = G_coeff
        self.L_coeff_even = L_coeff_even[:]
        self.L_coeff_odd = L_coeff_odd[:]
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
        
        # self.matrix_test()
        self.testing = False 
    
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
            self.make_MPRIME(xL, xR, dxL, dxR)
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
        """This function builds the mass matrix for spherical geometry"""
        rL2 = rL**2
        rR2 = rR**2
        rL3 = rL**3
        rR3 = rR**3
        rLrR = rL*rR
        pi = math.pi
        rttwo = math.sqrt(2)

        # if self.M ==0:
        #     self.Mass[0,0] = (rL2 + rLrR + rR2)/3/pi
        # elif self.M == 1:
        #     self.Mass[0,0] = (rL2 + rLrR + rR2)/3/pi
        #     self.Mass[1,0] = (rR2-rL2)/ 3/rttwo / pi 
        #     self.Mass[0,1] = (rR2-rL2)/ 3/rttwo / pi
        #     self.Mass[1,1] = 2*(2*rL2 + rLrR + 2*rR2)/15/pi
        # self.Mass[0,0] = 1/math.pi

        for ii in range(self.M+1):
            for jj in range(self.M+1):
                if (ii + jj + 2) % 2 == 0 or (ii == 0 and jj == 0):
                    self.Mass[ii][jj] = (self.Mass_coeff_even[ii, jj, 0] * rL2 + self.Mass_coeff_even[ii, jj, 1] * rLrR + self.Mass_coeff_even[ii, jj, 2] * rR2) * rtf(ii, jj) / pi 
                else:
                    self.Mass[ii][jj] = (self.Mass_coeff_odd[ii, jj, 0] * rL2 + self.Mass_coeff_odd[ii, jj, 1] * rR2 ) * rtf(ii, jj) / pi 
                

        self.Mass = np.multiply(self.Mass, 1/self.Mass_denom[0:self.M+1, 0:self.M+1])

        # testing mass matrix
        if self.testing == True:
            assert(abs(self.Mass[0,0] - (rL2 + rLrR + rR2)/3/pi <= 1e-10))
            if self.M>0:
                assert(abs(self.Mass[1,0] - (rR2-rL2)/ 3/rttwo / pi  <= 1e-10))
                assert(abs(self.Mass[0,1] - (rR2-rL2)/ 3/rttwo / pi  <= 1e-10))
                assert(abs(self.Mass[1,1] - 2*(2*rL2 + rLrR + 2*rR2)/15/pi  <= 1e-10))
            if self.M>1:
                assert(abs(self.Mass[2,2] - (30 * rL2 + 38 * rLrR + 30 * rR2)/105/pi)<=1e-10)

               
            

    def make_J_sphere(self, rL, rR):
        """This function builds the J matrix for the spherical case"""
        pi = math.pi
        rL2 = rL**2
        rR2 = rR**2
        rttwo = math.sqrt(2)
        # if self.M == 0:
        #     self.J[0,0] = 0.5 * (rR+rL) / pi
        # elif self.M == 1:
        #     self.J[0,0] = 0.5 * (rR+rL) /  pi
        #     self.J[1,0] = (rR-rL) /3 /rttwo/pi
        #     self.J[0,1] = (rR-rL) /3 /rttwo/pi
        #     self.J[1,1] =  (rR+rL) / pi / 3
        # elif self.M > 1:
        for ii in range(self.M+1):
            for jj in range(self.M+1):
                self.J[ii, jj] = (self.J_coeff[ii, jj, 0] * rL + self.J_coeff[ii, jj, 1] * rR)  / pi
        
        self.J[1:,0] = self.J[1:,0] / rttwo
        self.J[0,1:] = self.J[0,1:] / rttwo

        self.J = np.multiply(self.J, 1/self.J_denom[0:self.M+1, 0:self.M+1])
        if self.testing == True:
            assert(abs(self.J[0,0] - 0.5 * (rR+rL) /  pi)<=1e-10)
            
            if self.M >0:
                assert(abs(self.J[1,0] - (rR-rL) /3 /rttwo/pi)<=1e-10)
                assert(abs(self.J[0,1] - (rR-rL) /3 /rttwo/pi)<=1e-10)
                assert(abs(self.J[1,1] -  (rR+rL) / pi / 3)<=1e-10)
            if self.M > 1:
                assert(abs(self.J[2,2] - 7 * (rL+rR) /15 /pi)<=1e-10)

    
    


        # self.J[0,0] = math.log(rR/rL)/pi/(rR-rL)



    def make_G_sphere(self, rL, rR, rLp, rRp):
        """This function builds the G matrix for the spherical case"""
        rL2 = rL**2
        rR2 = rR**2
        rLrR = rL*rR
        pi = math.pi
        rttwo = math.sqrt(2)

        for ii in range(self.M+1):
            for jj in range(self.M+1):
                self.G[ii, jj] = GMatrix(ii, jj, rL, rR, rLp, rRp) / pi 
        
        # self.G[1:,0] = self.G[1:,0] / rttwo
        # self.G[0,1:] = self.G[0,1:] / rttwo

        # self.G = np.multiply(self.G, 1/self.G_denom[0:self.M+1, 0:self.M+1])


        if self.testing == True:
            a = rL
            b = rR
            ap = rLp
            bp = rRp
            assert(abs(self.G[0,0]  + 0.16666666666666666*((a**2 + a*b + b**2)*(ap - bp))/((a - b)*math.pi))<=1e-10)
            if self.M > 0:
                assert(abs(self.G[0,1] - ((a + b)*(ap - bp))/(6.*math.sqrt(2)*math.pi)) <=1e-10)
                assert(abs(self.G[1, 0] - (ap*(7*a**2 + 4*a*b + b**2) + (a**2 + 4*a*b + 7*b**2)*bp)/(6.*math.sqrt(2)*(a - b)*math.pi))<=1e-10)
                assert(abs(self.G[1,1] - (-(ap*(11*a**2 + 3*a*b + b**2)) + (a**2 + 3*a*b + 11*b**2)*bp)/(15.*(a - b)*math.pi)) <=1e-10)
            elif self.M > 1:
                assert(abs(self.G[2, 0] - (ap*(-7*a**2 - a*b + b**2) + (-a**2 + a*b + 7*b**2)*bp)/(3.*math.sqrt(2)*(a - b)*math.pi))<=1e-10)
                assert(abs(self.G[0,2] - ((a**2 + 3*a*b + b**2)*(ap - bp))/(15.*math.sqrt(2)*(a - b)*math.pi)))
        

    def make_L_sphere(self, rL, rR):
        """This function builds the L matrix for the spherical case"""
        pi = math.pi
        rL2 = rL**2
        rR2 = rR**2
        rLrR = rL*rR
        rtwo = math.sqrt(2)
        # if self.M == 0:
        #     self.L[0,0] = 0

        # if self.M == 1:
        #     self.L[0,0] = 0
        #     self.L[0,1] = 0
        #     self.L[1,0] = -2 * rtwo * (rL2 + rLrR + rR2) / 3/ pi /(rL-rR)
        #     self.L[1,1] = 2*(rL + rR)/3/pi

        for ii in range(1, self.M+1):
            for jj in range(self.M+1):
                if (ii + jj + 2) % 2 != 0:
                    self.L[ii, jj] = (self.L_coeff_odd[ii, jj, 0] * rL2 + self.L_coeff_odd[ii, jj, 1] * rLrR + self.L_coeff_odd[ii, jj, 2] * rR2) / (rL-rR) / pi
                else:
                    self.L[ii, jj] = (self.L_coeff_even[ii, jj, 0] * rL + self.L_coeff_even[ii, jj, 1] * rR)  / pi

        self.L[1:,0] = self.L[1:,0] * rtwo

        self.L[1:,:] = np.multiply(self.L[1:,:], 1/self.L_denom[1:self.M+1, 0:self.M+1])
        #testing L
        if self.testing == True:
            assert(self.L[0,0] == 0)
            if self.M >0:
                assert(self.L[0,1] == 0)
                assert(abs(self.L[1,0] + 2 * rtwo * (rL2 + rLrR + rR2) / 3/ pi /(rL-rR))<=1e-10)
                assert(abs(self.L[1,1] - 2*(rL + rR)/3/pi)<=1e-10)
            if self.M > 1:
                L2ac = -16 * (2* rL2 + rLrR + 2 * rR2) /15/pi/(rL-rR)
                assert(abs(self.L[2,1] - L2ac ) <= 1e-10)
            if self.M > 2:
                L33ac = 18 * (rR+rL) / 35/ pi
                assert(abs(self.L[3,3] - L33ac ) <= 1e-10)

    
    def make_MPRIME(self, a, b, ap, bp):
        pi = math.pi
        for ii in range(self.M+1):
            for jj in range(self.M+1):
                self.MPRIME[ii, jj] = MPRIME(ii, jj, a, b, ap, bp) / pi

        # self.L[0,0]  = 2*math.log(rR/rL)/pi/(rR-rL)
    
    def make_VV_sphere(self, rL, rR):
        """This function builds the VV matrix for the spherical case"""
        self.VV[0,0] = 0

    def matrix_test(self, test):
        self.testing = test
        if self.testing == True:
            a = 0.2
            b = 0.9
            ap = 0.9
            bp = 1.2
            self.make_all_matrices(a, b, ap, bp)

            if self.M == 3:
                L_bench = np.array([[0.0, 0.0, 0.0, 0.0],
                                    [0.441584, 0.233427, -0.168553, -0.140056],
                                    [0.660232, 0.911882, 0.186742, -0.465642],
                                    [0.609643, 0.980394, 1.17502, 0.180072]
                                    ])
                J_bench = np.array([[0.17507, 0.0525185, -0.082529, -0.0315111],
                                    [0.0525185, 0.116714, 0.0148545, -0.0700282],
                                    [-0.082529, 0.0148545, 0.163399, 0.031831],
                                    [-0.0315111, -0.0700282, 0.031831, 0.170068]
                                    ])
                Mass_bench = np.array([[0.109286, 0.0577703, -0.0417147, -0.0346622],
                                    [0.0577703, 0.0797897, 0.0163399, -0.0407437],
                                    [-0.0417147, 0.0163399, 0.0980394, 0.0350141],
                                    [-0.0346622, -0.0407437, 0.0350141, 0.105174]
                                    ])
                G_bench = np.array([[-0.0234185, -0.0123793, 0.00893885, 
                                    0.00742761], [-0.500801, -0.296392, 0.166476, 
                                    0.173252], [-0.781024, -1.04501, -0.250555, 
                                    0.493426], [-0.736684, -1.17347, -1.32831, -0.242067]])
                MPRIME_bench = np.array([
                                [0.378789 , 0.135047, -0.17016, -0.0810285],
                                [0.135047, 0.258468, 0.0381972, -0.148969],
                                [-0.17016, 0.0381972, 0.350141, 0.0818511],
                                [-0.0810285, -0.148969, 0.0818511, 0.367117]
                                ])
                
                if (np.abs(L_bench - self.L)>=1e-5).any():
                    print("L fail")
                    print(np.abs(L_bench - self.L))
                    assert(0)
                if (np.abs(J_bench - self.J)>=1e-5).any():
                    print("J fail")
                    assert(0)
                if (np.abs(Mass_bench - self.Mass)>=1e-5).any():
                    print("Mass fail")
                    assert(0)
                if (np.abs(G_bench - self.G)>=1e-5).any():
                    print("G fail")
                    print(np.abs(G_bench - self.G))
                    assert(0)
                if (np.abs(MPRIME_bench - self.MPRIME)>=1e-5).any():
                    print("M_prime fail")
                    print(np.abs(G_bench - self.G))
                    assert(0)
            self.testing = False