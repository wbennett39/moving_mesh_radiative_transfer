#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 11:25:35 2022

@author: bennett
"""
import numpy as np
import math

from .build_problem import build
from .matrices import G_L
from .sources import source_class
from .phi_class import scalar_flux
from .uncollided_solutions import uncollided_solution
from .numerical_flux import LU_surf
from .radiative_transfer import T_function
from .opacity import sigma_integrator
from .functions import shaper
from .functions import angular_deriv, finite_diff_uneven
import numba as nb

from numba.experimental import jitclass
from numba import int64, float64, deferred_type, prange
from numba import types, typed

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
transfer_class_type = deferred_type()
transfer_class_type.define(T_function.class_type.instance_type)
sigma_class_type = deferred_type()
sigma_class_type.define(sigma_integrator.class_type.instance_type)

# kv_ty = (types.int64, types.unicode_type)
params_default = nb.typed.Dict.empty(key_type=nb.typeof('par_1'),value_type=nb.typeof(1))


data = [('N_ang', int64), 
        ('N_space', int64),
        ('M', int64),
        ('source_type', int64[:]),
        ('t', float64),
        ('sigma_t', float64),
        ('sigma_s', float64),
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
        ("PV", float64[:]),
        ("S", float64[:]),
        ("LU", float64[:]),
        ("U", float64[:]),
        ("H", float64[:]),
        ("V_new", float64[:,:,:]),
        ("V", float64[:,:,:]),
        ("V_old", float64[:,:,:]),
        ('c', float64),
        ('uncollided', int64),
        ('thermal_couple', nb.typeof(params_default)),
        ('test_dimensional_rhs', int64),
        ('told', float64),
        ('division', float64),
        ('c_a', float64),
        ('sigma_a', float64),
        ('mean_free_time', float64),
        ('counter', int64),
        ('delta_tavg', float64),
        ('l', float64),
        ('times_list', float64[:]),
        ('save_derivative', int64),
        ('e_list', float64[:]),
        ('e_xs_list', float64[:]),
        ('wave_loc_list', float64[:]),
        ('sigma_func', nb.typeof(params_default)),
        ('particle_v', float64),
        ('epsilon', float64),
        ('deg_freedom', int64[:]),
        ('geometry', nb.typeof(params_default)),


        ]
##############################################################################
#thermal couple, l scaling


@jitclass(data)
class rhs_class():
    def __init__(self, build):
        self.N_ang = build.N_ang 
        self.N_space = build.N_space
        self.M = build.M
        self.mus = build.mus
        self.ws = build.ws
        self.source_type = np.array(list(build.source_type), dtype = np.int64) 
        self.c = build.sigma_s 
        self.thermal_couple = build.thermal_couple
        self.uncollided = build.uncollided
        self.test_dimensional_rhs = build.test_dimensional_rhs
        self.told = 0.0
        self.sigma_s = build.sigma_s
        self.sigma_a = build.sigma_a
        self.particle_v = build.particle_v
        self.geometry = build.geometry
       
        self.c_a = build.sigma_a / build.sigma_t
        self.mean_free_time = 1/build.sigma_t
        self.division = 1000
        self.counter = 0
        self.delta_tavg = 0.0
        self.l = build.l
        self.times_list = np.array([0.0])
        self.e_list = np.array([0.0])
        self.e_xs_list = np.array([0.0])
        self.wave_loc_list = np.array([0.0])
        self.save_derivative = build.save_wave_loc
        self.sigma_func = build.sigma_func
        self.deg_freedom = shaper(self.N_ang, self.N_space, self.M + 1, self.thermal_couple)

    
    def time_step_counter(self, t, mesh):
        delta_t = abs(self.told - t)
        self.delta_tavg += delta_t / self.division
        if self.counter == self.division:
            print('t = ', t, '|', 'delta_t average= ', self.delta_tavg)
            if self.N_space <= 32:
                if self.geometry['sphere'] == True:
                    print(mesh.edges[int(self.N_space/2):])
                else:
                    print(mesh.edges)
            print('--- --- --- --- --- --- --- --- --- --- --- --- --- ---')
            self.delta_tavg = 0.0
            self.counter = 0
        else:
            self.counter += 1
        self.told = t

        
    def derivative_saver(self, t,  space, transfer_class):
        if self.save_derivative == True:
            self.e_list = np.append(self.e_list, transfer_class.e_points)
            self.e_xs_list = np.append(self.e_xs_list, transfer_class.xs_points)

        if space == self.N_space - 1:
            deriv = np.copy(self.e_list)*0
            for ix in range(1,self.e_list.size-1):
                dx = self.e_xs_list[ix+1] - self.e_xs_list[ix]
                deriv[ix] = (self.e_list[ix+1] - self.e_list[ix])/dx

            max_deriv = max(np.abs(deriv))
            max_deriv_loc = np.argmin(np.abs(np.abs(self.e_list) - max_deriv))
            heat_wave_loc = self.e_xs_list[max_deriv_loc]
            self.wave_loc_list = np.append(self.wave_loc_list, abs(heat_wave_loc)) 
            self.times_list = np.append(self.times_list,t)
            # print(heat_wave_loc, 'wave x')
        
    def call(self, t, V, mesh, matrices, num_flux, source, uncollided_sol, flux, transfer_class, sigma_class):
        # print out timesteps
        self.time_step_counter(t, mesh) 
        # allocate arrays
        # V_new = V.copy().reshape((self.deg_freedom[0], self.deg_freedom[1], self.deg_freedom[2]))
        V_new = V.copy().reshape((self.N_ang, self.N_space, self.M+1))
        V_old = V_new.copy()
        # move mesh to time t 
        mesh.move(t)
        # represent opacity as a polynomial expansion
        sigma_class.sigma_moments(mesh.edges, t)
        flux.get_coeffs(sigma_class)

        # iterate over all cells
        for space in range(self.N_space):  
            # get mesh edges and derivatives          
            xR = mesh.edges[space+1]
            xL = mesh.edges[space]
            dxR = mesh.Dedges[space+1]
            dxL = mesh.Dedges[space]
            matrices.make_all_matrices(xL, xR, dxL, dxR)
            L = matrices.L
            G = matrices.G
             # special matrices for spherical geometries
            if self.geometry['sphere'] == True:
                Mass = matrices.Mass
                J = matrices.J
                VVs = matrices.VV
            # make P if there is any scattering
            if self.c != 0:
                flux.make_P(V_old[:,space,:], space, xL, xR)
            PV = flux.scalar_flux_term
            # integrate the source
            source.make_source(t, xL, xR, uncollided_sol)
            S = source.S

            # radiative transfer term
            # transfer_class.make_H(xL, xR, V_old[self.N_ang, space, :])
            # H = transfer_class.H
        

            ######### solve thermal couple ############

                
            ########## Loop over angle ############
            for angle in range(self.N_ang):
                mul = self.mus[angle]
                # calculate numerical flux
                refl_index = 0
                if angle <= self.N_ang/2 -1:
                    refl_index = self.N_ang-angle-1
                    assert(abs(self.mus[refl_index] - -self.mus[angle])<=1e-10)
                    
                num_flux.make_LU(t, mesh, V_old[angle,:,:], space, mul, V_old[refl_index, 0])
                
                    # print(self.mus[self.N_ang-angle-1], -self.mus[angle])
                    
                LU = num_flux.LU
                # Get absorption term
                sigma_class.make_vectors(mesh.edges, V_old[angle,space,:], space)
                VV = sigma_class.VV
                # Initialize solution vector, RHS
                U = np.zeros(self.M+1).transpose()
                U[:] = V_old[angle,space,:]
                # RHS = np.zeros_like(V_new[angle,space,:])

                # Spatial derivative term
               
                    
                # RHS = U*0
                # # Drift term
                # RHS += np.dot(G,U)
                # # numerical flux
                # RHS -=  LU
                # # Nspatial derivative
                # RHS += mul*np.dot(L,U)
                # # suource term
                # if self.geometry['sphere'] == True:
                #     RHS += S / 4 / math.pi
                # elif self.geometry['slab'] == True:
                #      RHS += 0.5*S 

                if self.geometry['sphere'] == True:
                    Minv = np.linalg.inv(Mass) 

                    if self.M == 0:
                        a = xL
                        b = xR
                        # print(Minv,  3 * math.pi/ (a*b + b**2 + a**2))
                        assert(np.abs(Minv[0,0] - 3 * math.pi/ (a*b + b**2 + a**2)) <=1e-8)
                    dterm = U*0
                    for j in range(self.M+1):
                        vec = (1-self.mus**2) * V_old[:, space, j]
                        dterm[j] = finite_diff_uneven(self.mus, angle, vec, left = (angle==0), right = (angle == self.N_ang - 1))

                    # Minv = np.copy(M)
                    # Minv[0,0] = 1/ M[0,0]
                #     RHS -= dterm * np.dot(J,U) 
                #     if self.M == 0:
                #         if abs(Minv[0,0] - 3/(xR**3-xL**3)) >1e-10:
                #             print(Minv[0,0])
                #             print(np.array([[3/(xR**3-xL**3)]]),'analytic')
                #             print('error')
                #     # Minv = np.array([[3/(xR**3-xL**3)]])
                #     # Minv = np.identity(self.M+1)
                #     # RHS = np.dot(Minv, RHS)
        
                # # Absorption term
                # RHS -= VV
                # # Scalar flux
                # RHS += PV

                # RHS = np.dot(G,U)  - LU + mul*np.dot(L,U) - VV + PV + 0.5*self.c*S
                # if space != 0:
                #     if abs(LU[0]) >1e-10:
                #         print(LU, 'Lu', space)
                #         print(VV, "VV")
                # RHS =  np.dot(G,U)  -  LU + mul*np.dot(L,U) - np.dot(M, VV) + np.dot(M, PV) +  0.5*self.c*S/4/math.pi - dterm * np.dot(J,U)
                     


                if self.geometry['sphere'] == True:
                    # V_new[angle,space,:] = np.dot(Minv,RHS)
                    # assert(Minv[0,0] - 3 * math.pi / (xL**2 + xL*xR + xR**2)) <= 1e-8
                    
                    # RHS = np.dot(G,U)  -  LU + mul*np.dot(L,U) - np.dot(Mass, U) + self.c / 4 / math.pi * np.dot(Mass, PV) + S/4/math.pi - dterm * np.dot(J,U)  
                    # V_new[angle,space,0] = RHS[0] *  3 * math.pi / (xL**2 + xL*xR + xR**2)
                    # RHS = self.c * PV / 4/ math.pi - U - 3 * (dterm * np.dot(J,U) + LU) * math.pi / (xL**2 + xL*xR + xR**2)
                    # V_new[angle,space,:] = np.dot(Minv, RHS)

                    a = xL
                    b = xR
                    RHS = V_old[angle, space, :]*0
                    
                    RHS -= np.dot(Minv, LU)
                    if self.M ==0:
                        assert(np.abs(RHS[0] + 3 * math.pi  * LU[0]  /  (a*b + b**2 + a**2) )<=1e-8)

                    RHS += mul*np.dot(L,U) * 0
                    RHS -= V_old[angle, space,:]

                    scalar_flux = np.sum(np.multiply(V_old[:, space, 0], self.ws)) * (self.c) 
                    PV = U*0 
                    PV[0] = scalar_flux
                    RHS += PV
                    mu_derivative = np.dot(Minv, np.dot(J, dterm))
                    if self.M == 0:
                        assert(abs(mu_derivative[0] -(3 * (a+b) / 2 / (a*b + b**2 + a**2)) * dterm[0])<=1e-8)
                    RHS -= mu_derivative

                    # RHS2 -= (3 * (a+b) / 2 / (a*b + b**2 + a**2)) * dterm 

                    # if self.M == 0:

                    # #hard code M=0 RHS
                    #     a = xL
                    #     b = xR
                    #     RHS2 = U*0
                    #     RHS2 -= 3 * math.pi  * LU  /  (a*b + b**2 + a**2)
                    #     if space == 0:
                    #         if mul > 0:
                    #             analytic_flux = (3 * math.sqrt(math.pi) / (a*b + b**2 + a**2) / math.sqrt(b-a) * (-mul) * V_old[angle, space, 0] * 1 / math.sqrt(b-a) / math.sqrt(math.pi) * b**2)
                    #             # print(analytic_flux)
                    #             # if abs(RHS[0] - analytic_flux) >= 1e-7:
                    #             #     print(abs(RHS[0] - analytic_flux) )
                    #             #     assert(0)
                    #     RHS2 -= (3 * (a+b) / 2 / (a*b + b**2 + a**2)) * dterm 
                    #     RHS2 -= V_old[angle, space,:] 
                    #     # assert(abs(np.sum(self.ws) - 1) <=1e-8)
                    #     PV = np.sum(np.multiply(V_old[:, space, 0], self.ws)) * (self.c /2 / math.pi) 
                    #     RHS2[0] +=  PV

                    #     assert(np.abs(RHS[0] - RHS2[0]) <= 1e-8)
                    
                    V_new[angle,space,:] = RHS  

        return V_new.reshape((self.N_ang) * self.N_space * (self.M+1))
        
    
       