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

from .build_problem import build
from .matrices import G_L
from .sources import source_class
from .phi_class import scalar_flux
from .uncollided_solutions import uncollided_solution
from .numerical_flux import LU_surf
from .radiative_transfer import T_function

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
        ("H", float64[:]),
        ("V_new", float64[:,:,:]),
        ("V", float64[:,:,:]),
        ("V_old", float64[:,:,:]),
        ('c', float64),
        ('uncollided', int64),
        ('thermal_couple', int64),
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
        
    def call(self,t, V, mesh, matrices, num_flux, source, uncollided_sol, flux, transfer_class):
        if self. thermal_couple == 0:
            V_new = V.copy().reshape((self.N_ang, self.N_space, self.M+1))
        elif self.thermal_couple == 1:
            V_new = V.copy().reshape((self.N_ang + 1, self.N_space, self.M+1))
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
            if self.source_type[4] != 1: # MMS source 
                source.make_source(t, xL, xR, uncollided_sol)
            if self.thermal_couple == 1:
                transfer_class.make_H(xL, xR, V_old[self.N_ang, space, :])
                H = transfer_class.H
            else: 
                H = np.zeros(self.M+1)
            S = source.S
            
            sigma_a = 1-self.c
            ######### solve thermal couple ############
            if self.thermal_couple == 1 and self.N_ang !=2:
        
                U = np.zeros(self.M+1).transpose()
                U[:] = V_old[self.N_ang,space,:]
                num_flux.make_LU(t, mesh, V_old[self.N_ang,:,:], space, 0.0)
                RU = num_flux.LU
                RHS_energy = np.dot(G,U) - RU + sigma_a * (2.0 * P  - H)
                
                if self.uncollided == True:
                    RHS_energy += sigma_a * source.S 
                V_new[self.N_ang ,space,:] = RHS_energy
                
            elif self.thermal_couple == 1 and self.N_ang ==2:
                U = np.zeros(self.M+1).transpose()
                U[:] = V_old[self.N_ang,space,:]
                num_flux.make_LU(t, mesh, V_old[self.N_ang,:,:], space, 0.0)
                RU = num_flux.LU
                RHS_energy = np.dot(G,U) - RU + sigma_a * (2.0 * P  - H)
                if self.uncollided == True:
                    RHS_energy += sigma_a * source.S 
                V_new[self.N_ang ,space,:] = RHS_energy
                
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
                    if self.uncollided == False:
                        RHS = np.dot(G,U)  - LU + mul*np.dot(L,U) - U + self.c * P + 0.5*S 
                    elif self.uncollided == True:
                        RHS = np.dot(G,U)  - LU + mul*np.dot(L,U) - U + self.c * (P + 0.5*S)
                    V_new[angle,space,:] = RHS
                    
                elif self.thermal_couple == 1:
                    deg_freedom = (self.N_ang + 1) * self.N_space * (self.M+1)
                    
                    if self.N_ang == 2:
                        if self.uncollided == True:
                            RHS_transport = np.dot(G,U) - LU + mul*np.dot(L,U) - U + self.c * (P + 0.5*S) + sigma_a*0.5*H
                        elif self.uncollided == False:
                            RHS_transport = np.dot(G,U) - LU + mul*np.dot(L,U) - U + self.c*P + 0.5*S + sigma_a*0.5*H
                    elif self.N_ang !=2:
                        if self.uncollided == True:
                            RHS_transport = np.dot(G,U) - LU + mul*np.dot(L,U) - U + self.c * (P + S*0.5) + sigma_a*0.5*H
                        elif self.uncollided == False:
                            RHS_transport = np.dot(G,U) - LU + mul*np.dot(L,U) - U + self.c * P + S*0.5 + sigma_a*0.5*H
                    V_new[angle,space,:] = RHS_transport 
                    
        return V_new.reshape(deg_freedom)
       
