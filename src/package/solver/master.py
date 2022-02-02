import numpy as np
import scipy.integrate as integrate
import quadpy

from build_problem import build
from matrices import G_L
from numerical_flux import LU_surf
from sources import source_class
from phi_class import scalar_flux
from mesh import mesh_class
from rhs_class import rhs_class
###############################################################################
""" 
[] have main take inputs from YAML 
[] make sigma_t and sigma_s dependent on space
[] uncollided source for square, truncated, square IC
[] source for no-uncollided cases
[] find uncollided for square source
[] make benchmarks for all cases 
[] njit all classes 
"""
###############################################################################

    
    

def main():
    
    tfinal = 1.0e-8
    angles = [2]
    Ms = [2]
    N_spaces = [2]
    x0 = 1e-12
    
    source_type = np.array([1,0,0,0])                                                     # ["plane", "square_IC", "square_source", "truncated_gaussian"]
    uncollided = True
    moving = True
    move_type = np.array([1,0,0,0])
    time = True 
    plotting = True
    RMS = True
    
    
    
    for count, N_space in enumerate(N_spaces):
        sigma_t = np.ones(N_space)
        sigma_s = np.ones(N_space)
        M = Ms[0]
        N_ang = angles[0]
        
        mus = quadpy.c1.gauss_lobatto(N_ang).points
        ws = quadpy.c1.gauss_lobatto(N_ang).weights
        xs_quad = quadpy.c1.gauss_lobatto(M).points
        ws_quad = quadpy.c1.gauss_lobatto(M).weights
        
        initialize = build(N_ang, N_space, M, tfinal, x0, mus, ws, xs_quad, ws_quad, sigma_t, sigma_s, source_type, uncollided, moving, move_type, time, plotting, RMS)
        initialize.make_IC()
        IC = initialize.IC
        mesh = mesh_class(N_space, x0, tfinal, moving, move_type) 
        matrices = G_L(initialize)
        num_flux = LU_surf(M)
        source = source_class(initialize)
        flux = scalar_flux(initialize)
        rhs = rhs_class(initialize)
        RHS = lambda t, V: rhs(t, V, mesh, matrices, num_flux, source, flux)
        sol = integrate.solve_ivp(RHS, [0.0,tfinal], IC.reshape(N_ang*N_space*(M+1)), method='DOP853')

        
        
        
main()
        
        
        
    
    
    
    

            





