import numpy as np
import scipy.integrate as integrate
import quadpy
import matplotlib.pyplot as plt
from timeit import default_timer as timer
import yaml
from pathlib import Path

from .solver_classes.build_problem import build
from .solver_classes.matrices import G_L
from .solver_classes.numerical_flux import LU_surf
from .solver_classes.sources import source_class
from .solver_classes.uncollided_solutions import uncollided_solution
from .solver_classes.phi_class import scalar_flux
from .solver_classes.mesh import mesh_class
from .solver_classes.rhs_class import rhs_class
from .solver_classes.make_phi import make_output
from .solver_classes.functions import find_nodes, convergence
from .solver_classes.save_output import save_output
from .load_bench import load_bench
###############################################################################
""" 

class goals:
[x] have main take inputs from YAML 
[] make sigma_t and sigma_s dependent on space
[x] jitclass RHS
        is it any faster?
[x] fix numerical flux class
[] put imports __init__ file
[x] figure out where that factor of two comes from in the source
[] uncollided_solution -- make temp[ix] self
[] comments
[] usability - make functions that are used for IC, benchmarking, source, mesh, easy to modify 
[] pytest
[x] bug in hp hm for MMS
[] clean up file structure
    ideas for tests:
        source - check if integrator returns analytic integral of plane uncollided for x inside t 
        G, and L, check against mathematica result
[x] I should probably make a class that just returns uncollided solutions
[] make new benchmark maker into a class
[] benchmarks for t 1= , 5, 10 not just 1
[] fix h5 algorithm on benchmark
[] incorporate the experimental scripts
[] clean up main with auxilliary functions

paper goals:
[x] static mesh plane uncollided
[x] run benchmarks with more evaluation points    
[] mesh function to better capture square,, delta function
[] mesh function to move with uncollided sol for square, truncated <--------
[x] uncollided source for square IC
[x] fix find nodes
[x] find uncollided for square source
[] fix square source uncollided 
[x] uncollided for gaussian
[x] uncollided for gaussian source
[x] MMS
    [x] IC for plane -- static probably needs special mesh
    [x] IC for square 
    [x] IC for gaussian
[] make benchmarks for 
    [x] plane source
    [x] square IC
    [x] gaussian IC
    [x] gauss source
    [x] square source
[x] timing
[x] something is going on with square IC with uncollided
[x] is the benchmark maker off by a little bit?
[] no-uncollided plane source doesnt converge -- take it out
[] x0 function fort the gaussian and plane 
 
"""
###############################################################################
data_folder = Path("package")
config_file_path = data_folder / "config.yaml"
###############################################################################


def run_plane(uncollided = True, moving = True):
    source_name = "plane_IC"
    print("running plane IC")
    main(source_name, uncollided, moving)
    
def run_all():
    return 0
    
def main(source_name, uncollided, moving):
    
    with open(config_file_path, 'r') as file:
       parameters = yaml.safe_load(file)
    file.close()
    #'C:/users/bennett/Documents/GitHub/MovingMesh/moving_mesh_radiative_transfer/src/package/
    tfinal = float(parameters['all']['tfinal'])
    N_spaces = np.array(parameters['all']['N_spaces'])
    Ms = np.array(parameters['all']['Ms'])
    N_runs = int(parameters['all']['N_runs'])
    t_nodes = int(parameters['all']['tnodes'])
    rt = float(parameters['all']['rt'])
    at = float(parameters['all']['at'])
    
    N_angles = np.array(parameters[source_name]['N_angles'])
    x0 = float(parameters[source_name]['x0'])
    source_type = np.array(parameters[source_name]['source_type'])
    move_type = np.array(parameters[source_name]['move_type'])
    
    
    r_times = np.zeros(len(N_angles))
    RMS_list = []
    x0s = np.ones(4)*x0
    
    saving = save_output(tfinal, Ms[0], source_type, moving, uncollided)
    benchmark = load_bench(source_type, tfinal, x0)
    
    plt.figure(11)
    if source_type[0] == 1:
        xsb2 = np.linspace(0, tfinal , 1000)
        xsb = np.concatenate((xsb2, np.ones(1)*1.0000001))
        bench = np.concatenate((benchmark(xsb2),np.zeros(1)))
    else:
        xsb = np.linspace(0, tfinal + x0, 100000)
        bench = benchmark(xsb)
        
    plt.plot(xsb, bench, "k-")
    plt.plot(-xsb, bench, "k-")
    
    print("uncollided", uncollided)
    print("moving", moving)
    for nr in range(N_runs):
        for count, N_space in enumerate(N_spaces):
            sigma_t = np.ones(N_space)
            sigma_s = np.ones(N_space)
            M = Ms[0]
            N_ang = N_angles[count]
            if source_type[0] == 1 and uncollided == False and moving == True:
                x0 = x0s[count]/N_space
            mus = quadpy.c1.gauss_lobatto(N_ang).points
            ws = quadpy.c1.gauss_lobatto(N_ang).weights
            xs_quad = quadpy.c1.gauss_legendre(M+2).points
            ws_quad = quadpy.c1.gauss_legendre(M+2).weights
            t_quad = quadpy.c1.gauss_legendre(t_nodes).points
            t_ws = quadpy.c1.gauss_legendre(t_nodes).weights
            initialize = build(N_ang, N_space, M, tfinal, x0, mus, ws, xs_quad, ws_quad, sigma_t, sigma_s, source_type, uncollided, moving, move_type, t_quad, t_ws)
            initialize.make_IC()
            IC = initialize.IC
            mesh = mesh_class(N_space, x0, tfinal, moving, move_type) 
            matrices = G_L(initialize)
            num_flux = LU_surf(initialize)
            source = source_class(initialize)
            uncollided_sol = uncollided_solution(initialize)
            flux = scalar_flux(initialize)
            rhs = rhs_class(initialize)
            def RHS(t, V):
                return rhs.call(t, V, mesh, matrices, num_flux, source, uncollided_sol, flux)
            
            start = timer()
            sol = integrate.solve_ivp(RHS, [0.0,tfinal], IC.reshape(N_ang*N_space*(M+1)), method='DOP853', t_eval = [tfinal], rtol = rt, atol = at)
            end = timer()
            r_times[nr] += (end-start)/N_runs
            sol_last = sol.y[:,-1].reshape((N_ang,N_space,M+1))
            mesh.move(tfinal)
            edges = mesh.edges
            
            xs = find_nodes(edges, M)
            output = make_output(tfinal, N_ang, ws, xs, sol_last, M, edges, uncollided)
            phi = output.make_phi(uncollided_sol)
            
            benchmark_solution = benchmark(np.abs(xs))
            RMS = np.sqrt(np.mean((phi - benchmark_solution)**2))
            RMS_list.append(RMS)
            print(N_space, "spaces", "    ", "%.4f" % (end-start), "time elapsed")
            print("RMSE", RMS)
            if count > 0:
                print("Order", "%.2f" % convergence(RMS_list[count-1], N_spaces[count-1], RMS, N_space))
            # phi = make_phi(N_ang, ws, xs, sol_last, M, edges) 
            plt.figure(11)
            # plt.plot(xs, phi, "-o")
            plt.xlabel("x")
            plt.ylabel("scalar flux")
            plt.plot(xs, phi, "-o")
            # plt.plot(xs, benchmark_solution, "k-")

            
    saving.save_RMS(RMS_list, N_spaces, N_angles, r_times)
        
        
# run_plane()       
        
    
    
    
    

            





