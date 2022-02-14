import numpy as np
import scipy.integrate as integrate
import quadpy
import matplotlib.pyplot as plt
import math

from build_problem import build
from matrices import G_L
from numerical_flux import LU_surf
from sources import source_class
from phi_class import scalar_flux
from mesh import mesh_class
from rhs_class import rhs_class
from make_phi import make_output
from functions import make_phi, find_nodes
from load_bench import load_bench # fix later 
###############################################################################
""" 

class goals:
[] have main take inputs from YAML 
[] make sigma_t and sigma_s dependent on space
[] jitclass RHS
[] fix numerical flux class
[] put imports __init__ file
[] figure out where that factor of two comes from in the source
[] comments
[] usability - make functions that are used for IC, benchmarking, source, mesh, easy to modify 
[] pytest
[] put benchmarks back in its own module
    ideas for tests:
        source - check if integrator returns analytic integral of plane uncollided for x inside t 
        G, and L, check against mathematica result

paper goals:
[x] static mesh plane uncollided
[] run benchmarks with more evaluation points    
[] mesh function to better capture square, truncated IC, delta function
[] mesh function to move with uncollided sol for square, truncated 
[x] uncollided source for square IC
[x] fix find nodes
[] find uncollided for square source
[] uncollided for truncated
[] uncollided for truncated source
[] MMS
    [x] IC for plane -- static probably needs special mesh
    [x] IC for square 
    [] IC for truncated
[] make benchmarks for 
    [x] plane source
    [x] square IC
    [x] truncated gaussian IC
    [] trunc gauss source
    [] square source
[] timing
[x] something is going on with square IC with uncollided
[x] is the benchmark maker off by a little bit?

"""
###############################################################################

def main():
    
    tfinal = 1.0
    angles = [32]
    Ms = [4]
    N_spaces = [2,4]
    # x0 = 1e-10
    x0 = 1/2
    source_type = np.array([0,1,0,0])                                                     # ["plane", "square_IC", "square_source", "truncated_gaussian"]
    uncollided = True
    moving = False
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
        xs_quad = quadpy.c1.gauss_legendre(M+2).points
        ws_quad = quadpy.c1.gauss_legendre(M+2).weights
        
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
        sol = integrate.solve_ivp(RHS, [0.0,tfinal], IC.reshape(N_ang*N_space*(M+1)), method='DOP853', t_eval = [tfinal])
        sol_last = sol.y[:,-1].reshape((N_ang,N_space,M+1))
        mesh.move(tfinal)
        edges = mesh.edges
        
        xs = find_nodes(edges, M)
        output = make_output(tfinal, N_ang, ws, xs, sol_last, M, edges, uncollided)
        phi = output.make_phi(source)
        
        benchmark = load_bench(source_type, tfinal)
        benchmark_solution = benchmark(xs)
        RMS = np.linalg.norm(phi - benchmark_solution)
        
        # phi = make_phi(N_ang, ws, xs, sol_last, M, edges) 
        plt.plot(xs, phi, "-o")
        plt.plot(xs, benchmark_solution, "k-")
        
        
main()
        
        
        
    
    
    
    

            





