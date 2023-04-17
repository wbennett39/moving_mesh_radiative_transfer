# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# """
# Created on Mon Jan 31 11:25:35 2022

# @author: bennett
# """
# import numpy as np

# from .build_problem import build
# from .matrices import G_L
# from .sources import source_class
# from .phi_class import scalar_flux
# from .uncollided_solutions import uncollided_solution
# from .numerical_flux import LU_surf
# from .radiative_transfer import T_function

# from numba.experimental import jitclass
# from numba import int64, float64, deferred_type, prange

# build_type = deferred_type()
# build_type.define(build.class_type.instance_type)
# matrices_type = deferred_type()
# matrices_type.define(G_L.class_type.instance_type)
# num_flux_type = deferred_type()
# num_flux_type.define(LU_surf.class_type.instance_type)
# source_type = deferred_type()
# source_type.define(source_class.class_type.instance_type)
# flux_type = deferred_type()
# flux_type.define(scalar_flux.class_type.instance_type)
# uncollided_solution_type = deferred_type()
# uncollided_solution_type.define(uncollided_solution.class_type.instance_type)
# transfer_class_type = deferred_type()
# transfer_class_type.define(T_function.class_type.instance_type)


# data = [('N_ang', int64), 
#         ('N_space', int64),
#         ('M', int64),
#         ('source_type', int64[:]),
#         ('t', float64),
#         ('sigma_t', float64[:]),
#         ('sigma_s', float64[:]),
#         ('IC', float64[:,:,:]),
#         ('mus', float64[:]),
#         ('ws', float64[:]),
#         ('x0', float64),
#         ("xL", float64),
#         ("xR", float64),
#         ("dxL", float64),
#         ("dxR", float64),
#         ("L", float64[:,:]),
#         ("G", float64[:,:]),
#         ("P", float64[:]),
#         ("S", float64[:]),
#         ("LU", float64[:]),
#         ("U", float64[:]),
#         ("H", float64[:]),
#         ("V_new", float64[:,:,:]),
#         ("V", float64[:,:,:]),
#         ("V_old", float64[:,:,:]),
#         ('c', float64),
#         ('uncollided', int64),
#         ('thermal_couple', int64),
#         ]
# ##############################################################################
# @jitclass(data)
# class rhs_class():
#     def __init__(self, build):
#         self.N_ang = build.N_ang 
#         self.N_space = build.N_space
#         self.M = build.M
#         self.mus = build.mus
#         self.ws = build.ws
#         self.source_type = build.source_type
#         self.c = build.scattering_ratio
#         self.thermal_couple = build.thermal_couple
#         # self.temperature_function = build.temperature_function
#         self.uncollided = build.uncollided
        
#     def call(self,t, V, mesh, matrices, num_flux, source, uncollided_sol, flux, transfer_class):
#         if self. thermal_couple == 0:
#             V_new = V.copy().reshape((self.N_ang, self.N_space, self.M+1))
#         elif self.thermal_couple == 1:
#             V_new = V.copy().reshape((self.N_ang + 1, self.N_space, self.M+1))
#         V_old = V_new.copy()
#         mesh.move(t)

#         for space in prange(self.N_space):            
#             xR = mesh.edges[space+1]
#             xL = mesh.edges[space]
#             dxR = mesh.Dedges[space+1]
#             dxL = mesh.Dedges[space]
#             matrices.make_L(xL, xR)
#             matrices.make_G(xL, xR, dxL, dxR)
#             L = matrices.L
#             G = matrices.G
#             flux.make_P(V_old[:,space,:])
#             P = flux.P
#             if self.source_type[4] != 1: # MMS source 
#                 source.make_source(t, xL, xR, uncollided_sol)
#             if self.thermal_couple == 1:
#                 transfer_class.make_H(xL, xR, V_old[self.N_ang, space, :])
#                 H = transfer_class.H
                
#             else: 
#                 H = np.zeros(self.M+1)
#             S = source.S
            
#             if self.uncollided == True:
#                 c2 = self.c 
#             else:
#                 c2 = 1
#             ######### solve thermal couple ############
#             if self.thermal_couple == 1:
#                 sigma_a = 1-self.c
#                 U = np.zeros(self.M+1).transpose()
#                 U[:] = V_old[self.N_ang,space,:]
#                 num_flux.make_LU(t, mesh, V_old[self.N_ang,:,:], space, 0.0)
#                 RU = num_flux.LU

#                 # RHS_energy = U*0
#                 RHS_energy = np.dot(G,U) - RU + sigma_a * (2.0 * P  - H)
#                 if self.uncollided == True:
#                     RHS_energy += sigma_a * source.S 
#                 V_new[self.N_ang ,space,:] = RHS_energy
                
#             ########## Loop over angle ############
#             for angle in range(self.N_ang):
#                 mul = self.mus[angle]
        
#                 if self.source_type[4] == 1: # Make MMS source
#                     source.make_source_not_isotropic(t, mul, xL, xR)
                    
#                 num_flux.make_LU(t, mesh, V_old[angle,:,:], space, mul)
#                 LU = num_flux.LU
                
#                 U = np.zeros(self.M+1).transpose()
#                 U[:] = V_old[angle,space,:]
#                 if self.thermal_couple == 0:
#                     deg_freedom = self.N_ang * self.N_space * (self.M+1)
#                     RHS = np.dot(G,U)  - LU + mul*np.dot(L,U) - U + self.c * P + c2*0.5*S 
#                     V_new[angle,space,:] = RHS
#                 elif self.thermal_couple == 1:
#                     deg_freedom = (self.N_ang + 1) * self.N_space * (self.M+1)
#                     RHS_transport = np.dot(G,U) - LU + mul*np.dot(L,U) - U + self.c*P + c2*0.5*S + sigma_a*0.5*H
#                     V_new[angle,space,:] = RHS_transport 
#         return V_new.reshape(deg_freedom)
           



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
transfer_class_type = deferred_type()
transfer_class_type.define(T_function.class_type.instance_type)
sigma_class_type = deferred_type()
sigma_class_type.define(sigma_integrator.class_type.instance_type)


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
        ('thermal_couple', int64),
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
        ('sigma_func', int64[:])
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
        self.c = build.scattering_ratio
        self.thermal_couple = build.thermal_couple
        self.uncollided = build.uncollided
        self.test_dimensional_rhs = False
        self.told = 0.0
        self.sigma_s = build.sigma_s
       
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

    
    def time_step_counter(self, t, mesh):
        delta_t = abs(self.told - t)
        self.delta_tavg += delta_t / self.division
        if self.counter == self.division:
            print('t = ', t, '|', 'delta_t average= ', self.delta_tavg)
            if self.N_space <= 32:
                print(mesh.edges[int(self.N_space/2):])
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
        
    def call(self,t, V, mesh, matrices, num_flux, source, uncollided_sol, flux, transfer_class, sigma_class):
        self.time_step_counter(t, mesh)
        if self. thermal_couple == 0:
            V_new = V.copy().reshape((self.N_ang, self.N_space, self.M+1))
        elif self.thermal_couple == 1:
            V_new = V.copy().reshape((self.N_ang + 1, self.N_space, self.M+1))
        V_old = V_new.copy()
        mesh.move(t)
        sigma_class.sigma_moments(mesh.edges)
        flux.get_coeffs(sigma_class)

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

            if (self.sigma_func[0] == 1) or (self.c == 0.0):
                P = flux.P
            else:
                flux.make_P_nonconstant_opacity(V_old[:, space, :], space)


            if self.source_type[4] != 1: # MMS source 
                source.make_source(t, xL, xR, uncollided_sol)
            if self.thermal_couple == 1:
                transfer_class.make_H(xL, xR, V_old[self.N_ang, space, :])
                H = transfer_class.H
            else: 
                H = np.zeros(self.M+1)

            self.derivative_saver(t, space, transfer_class)

            S = source.S

            ######### solve thermal couple ############
            if self.thermal_couple == 1:
                U = np.zeros(self.M+1).transpose()
                U[:] = V_old[self.N_ang,space,:]
                num_flux.make_LU(t, mesh, V_old[self.N_ang,:,:], space, 0.0)
                RU = num_flux.LU
                if self.test_dimensional_rhs == True:
                    RHS_energy = np.dot(G,U) - RU + self.c_a * (2.0 * P  - H)
                else:
                    RHS_energy = np.dot(G,U) - RU + self.c_a * (2.0 * P  - H) / self.l
                
                if self.uncollided == True:
                    RHS_energy += self.c_a * source.S /self.l
                V_new[self.N_ang,space,:] = RHS_energy
                
            # elif self.thermal_couple == 1 and self.N_ang == 2:
            #     U = np.zeros(self.M+1).transpose()
            #     U[:] = V_old[self.N_ang,space,:]
            #     num_flux.make_LU(t, mesh, V_old[self.N_ang,:,:], space, 0.0)
            #     RU = num_flux.LU
            #     RHS_energy = np.dot(G,U) - RU + self.c_a * (2.0 * P  - H) / self.l
            #     if self.uncollided == True:
            #         RHS_energy += self.c_a * source.S / self.l
            #     V_new[self.N_ang ,space,:] = RHS_energy
                
            ########## Loop over angle ############
            for angle in range(self.N_ang):
                mul = self.mus[angle]
                if self.source_type[4] == 1: # Make MMS source
                    source.make_source_not_isotropic(t, mul, xL, xR)
                num_flux.make_LU(t, mesh, V_old[angle,:,:], space, mul)
                LU = num_flux.LU
                U = np.zeros(self.M+1).transpose()
                U[:] = V_old[angle,space,:]
                
                if self.thermal_couple == 0:
                    
                    deg_freedom = self.N_ang * self.N_space * (self.M+1)
                    if self.sigma_func[0] == 1:
                        if self.uncollided == False:
                            RHS = np.dot(G,U)  - LU + mul*np.dot(L,U) - U + self.c * P + 0.5*S 
                        elif self.uncollided == True:
                            RHS = np.dot(G,U)  - LU + mul*np.dot(L,U) - U + self.c * (P + 0.5*S)

                    elif self.sigma_func[2]== 1: #siewert problem
                        VV = sigma_class.make_vectors(mesh.edges, V_old[angle,space,:], space)
                        PV = flux.call_P_noncon(xL, xR)
                        # Q = np.zeros(PV.size) # this is for testing the steady state source problem
                        # Q[0] = math.sqrt(xR-xL)
                        RHS = np.dot(G,U)  - LU + mul*np.dot(L,U) - U + PV
                    
                    else:
                        VV = sigma_class.make_vectors(mesh.edges, V_old[angle,space,:], space)
                        PV = flux.call_P_noncon(xL, xR)
                        # PV =  self.sigma_s*flux.P
                        # PV = VV*0
                        # print(np.abs(PV-flux.P))
                        assert(np.abs(flux.cs[space,:] - sigma_class.cs[space,:]).all() <= 1e-10)
                        # if (np.abs(self.sigma_s * flux.P - PV).all() > 1e-6):
                        #     print(flux.P - PV)
                        RHS = np.dot(G,U)  - LU + mul*np.dot(L,U) - VV + PV

                    V_new[angle,space,:] = RHS
                    
                elif self.thermal_couple == 1:

                    deg_freedom = (self.N_ang + 1) * self.N_space * (self.M+1)
                    
                    if self.N_ang == 2:
                        if self.uncollided == True:
                            RHS_transport = np.dot(G,U) - LU + mul*np.dot(L,U) - U/self.l + self.c * (P/self.l + 0.5*S/self.l) + self.c_a*0.5*H/self.l
                        elif self.uncollided == False:
                            RHS_transport = np.dot(G,U) - LU + mul*np.dot(L,U) - U/self.l + self.c*P/self.l + 0.5*S/self.l + self.c_a*0.5*H/self.l
                    elif self.N_ang !=2:
                        if self.uncollided == True:
                            RHS_transport = np.dot(G,U) - LU + mul*np.dot(L,U) - U/self.l + self.c * (P + S*0.5)/self.l + self.c_a*0.5*H/self.l
                        elif self.uncollided == False:
                            if self.test_dimensional_rhs == True:
                                RHS_transport = np.dot(G,U) - LU + 299.98*mul*np.dot(L,U) - 299.98*U + 299.98*self.c * P + 299.98 * S*0.5 + 299.98*self.c_a*0.5*H
                            else:
                                RHS_transport = np.dot(G,U) - LU + mul*np.dot(L,U) - U/self.l + self.c * P /self.l + S*0.5/self.l + self.c_a*0.5*H/self.l
                    V_new[angle,space,:] = RHS_transport 
                    
        return V_new.reshape(deg_freedom)
    
       
