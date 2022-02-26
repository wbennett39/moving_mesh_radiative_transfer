import numpy as np
import scipy.integrate as integrate
import quadpy
import matplotlib.pyplot as plt
import math
from numba import njit
from timeit import default_timer as timer

from build_problem import build
from matrices import G_L
from numerical_flux import LU_surf
from sources import source_class
from phi_class import scalar_flux
from mesh import mesh_class
from rhs_class import rhs_class
from make_phi import make_output
from functions import make_phi, find_nodes, convergence
from load_bench import load_bench # fix later 
from save_output import save_output
###############################################################################
""" 

class goals:
[] have main take inputs from YAML 
[] make sigma_t and sigma_s dependent on space
[] jitclass RHS
[x] fix numerical flux class
[] put imports __init__ file
[] figure out where that factor of two comes from in the source
[] comments
[] usability - make functions that are used for IC, benchmarking, source, mesh, easy to modify 
[] pytest
[] bug in hp hm for MMS
[] put benchmarks back in its own module
    ideas for tests:
        source - check if integrator returns analytic integral of plane uncollided for x inside t 
        G, and L, check against mathematica result


paper goals:
[x] static mesh plane uncollided
[x] run benchmarks with more evaluation points    
[] mesh function to better capture square, truncated IC, delta function
[] mesh function to move with uncollided sol for square, truncated 
[x] uncollided source for square IC
[x] fix find nodes
[] find uncollided for square source
[x] uncollided for gaussian
[] uncollided for gaussian source
[x] MMS
    [x] IC for plane -- static probably needs special mesh
    [x] IC for square 
    [x] IC for gaussian
[] make benchmarks for 
    [x] plane source
    [x] square IC
    [x] gaussian IC
    [] gauss source
    [] square source
[x] timing
[x] something is going on with square IC with uncollided
[x] is the benchmark maker off by a little bit?
[] no-uncollided plane source doesnt converge


"""
###############################################################################

def main(uncollided = True, moving = True):
    N_runs = 1
    rt = 1e-9
    at = 1e-7
    tfinal = 1.0
    angles = [256,256,256]
    # pars for no uncollided moving plane
    #angles [256, 512 ]
    # tols 1e-9, 1e-7 -- 1e-9, 1e-7, -- 
    # angles = [2,2,2,2]
    r_times = np.zeros(len(angles))
    Ms = [4]
    N_spaces = [2,4,8]
    RMS_list = []
    # x0 = 1e-5
    x0 = 0.5
    x0s = np.ones(4)*x0
    source_type = np.array([0,0,1,0,0])                                                     # ["plane", "square_IC", "square_source", "gaussian", "MMS"]
    move_type = np.array([1,0,0,0])
    time = True 
    plotting = True
    RMS = True
    saving = save_output(tfinal, Ms[0], source_type, moving, uncollided)
    benchmark = load_bench(source_type, tfinal, x0)
    
    plt.figure(10)
    if source_type[0] == 1:
        xsb = np.linspace(0, tfinal , 100000)
    else:
        xsb = np.linspace(0, tfinal + x0, 100000)
        
    plt.plot(xsb, benchmark(xsb), "k-")
    plt.plot(-xsb, benchmark(xsb), "k-")
    
    print("uncollided", uncollided)
    print("moving", moving)
    for nr in range(N_runs):
        for count, N_space in enumerate(N_spaces):
            sigma_t = np.ones(N_space)
            sigma_s = np.ones(N_space)
            M = Ms[0]
            N_ang = angles[count]
            if source_type[0] == 1 and uncollided == False and moving == True:
                x0 = x0s[count]/N_space
                print(x0)
            mus = quadpy.c1.gauss_lobatto(N_ang).points
            ws = quadpy.c1.gauss_lobatto(N_ang).weights
            xs_quad = quadpy.c1.gauss_legendre(M+2).points
            ws_quad = quadpy.c1.gauss_legendre(M+2).weights
            
            initialize = build(N_ang, N_space, M, tfinal, x0, mus, ws, xs_quad, ws_quad, sigma_t, sigma_s, source_type, uncollided, moving, move_type, time, plotting, RMS)
            initialize.make_IC()
            IC = initialize.IC
            mesh = mesh_class(N_space, x0, tfinal, moving, move_type) 
            matrices = G_L(initialize)
            num_flux = LU_surf(initialize)
            source = source_class(initialize)
            flux = scalar_flux(initialize)
            rhs = rhs_class(initialize)
            # RHS = lambda t, V: rhs(t, V, mesh, matrices, num_flux, source, flux)
            # @njit
            def RHS(t, V):
                return rhs.call(t, V, mesh, matrices, num_flux, source, flux)
            
            start = timer()
            sol = integrate.solve_ivp(RHS, [0.0,tfinal], IC.reshape(N_ang*N_space*(M+1)), method='DOP853', t_eval = [tfinal], rtol = rt, atol = at)
            end = timer()
            r_times[nr] += (end-start)/N_runs
            sol_last = sol.y[:,-1].reshape((N_ang,N_space,M+1))
            mesh.move(tfinal)
            edges = mesh.edges
            
            xs = find_nodes(edges, M)
            output = make_output(tfinal, N_ang, ws, xs, sol_last, M, edges, uncollided)
            phi = output.make_phi(source)
            
            benchmark_solution = benchmark(np.abs(xs))
            RMS = np.sqrt(np.mean((phi - benchmark_solution)**2))
            RMS_list.append(RMS)
            print(N_space, "spaces", "    ", "%.4f" % (end-start), "time elapsed")
            print("RMSE", RMS)
            if count > 0:
                print("Order", "%.2f" % convergence(RMS_list[count-1], N_spaces[count-1], RMS, N_space))
            
            # phi = make_phi(N_ang, ws, xs, sol_last, M, edges) 
            plt.figure(10)
            plt.plot(xs, phi, "-o")
            # plt.plot(xs, benchmark_solution, "k-")
            
    saving.save_RMS(RMS_list, N_spaces, angles, r_times)
        
# main()
main(uncollided = False)
# main(moving = False)
# main(uncollided = False, moving = False)
        
        
        
    
    
    
    

            





