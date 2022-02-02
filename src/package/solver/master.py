import numpy as np
import scipy.integrate as integrate

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
    
    source_type = "plane"
    uncollided = True
    moving = "linear"
    time = True 
    plotting = True
    
    for count, N_space in enumerate(N_spaces):
        sigma_t = np.ones(N_space)
        sigma_s = np.ones(N_space)
        M = Ms[0]
        N_ang = angles[0]
        
        initialize = build(N_ang, N_space, M, tfinal, x0, sigma_t, sigma_s, source_type, uncollided, moving, time, plotting)
        initialize.make_IC()
        IC = initialize.IC
        mesh = mesh_class(N_space, x0, source_type, tfinal, uncollided, moving) 
        matrices = G_L(initialize)
        num_flux = LU_surf(M)
        source = source_class(initialize)
        flux = scalar_flux(initialize)
        rhs = rhs_class(initialize)
        RHS = lambda t, V: rhs(t, V, mesh, matrices, num_flux, source, flux)
        sol = integrate.solve_ivp(RHS, [0.0,tfinal], IC.reshape(N_ang*N_space*(M+1)), method='DOP853')

        
        
        
main()
        
        
        
    
    
    
    

            





