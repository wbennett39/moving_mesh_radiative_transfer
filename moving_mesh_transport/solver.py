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
from .save_output import save_output
from .load_bench import load_bench
###############################################################################
""" 

class goals:

[x] license
[x] write benchmarking module that evaluates and saves at t = 1, t = 5, and t = 10, renames the experimental scripts
    [x] run the benchmark maker
[x] write plotting function that reads RMS data and makes plots for paper and lives in a \plots folder
    [x] err vs space
    [x] err vs computation time
    [x] err vs best case avg
    [x] benchmark solution 
[x] write plotting function that plots results from run_plane, run_square_IC etc.
[x] find best parameters for each problem type, put into input file 
[] fix time saving 
[x] plot bench last
[x] pytest
[] better convergence triangle
[] Sn labels on plot?
[] update README

[x] no-uncollided plane source doesnt converge -- take it out
[] write report 
[] pull request for report and code

square source:
    -solve uncollided equation vs uncollided bench to see if it converges and confirm the uncollided solution is correct
    -check if the triple integral can be simplified at all
    -maybe <= instead of < in F1
    -cancellation in F1 integrand

ideas for tests:
    source - check if integrator returns analytic integral of plane uncollided for x inside t 
    G, and L, check against integrated normPns -
    check normPns

long term goals:
[] comments for all classes, functions
[] mu mapping for the cases that have smooth but discontinuous $\psi_u$
[] uncollided_solution -- make temp[ix] self
[] x0 function for the gaussian and plane 
[] benchmarking class
[] sigma_t and sigma_s functions of space
[] usability - make functions that are used for IC, benchmarking, source, mesh, easy to modify 
[] mesh function to better capture square,, delta function
[] mesh function to move with uncollided sol for square, truncated <--------
[] ability to automatically add other benchmark times
[] is jitclass RHS any faster?
[] pycuda for integrator
[] finite domain for ganapol? 

"""
###############################################################################
data_folder = Path("moving_mesh_transport")
config_file_path = data_folder / "config.yaml"
###############################################################################


def run_plane_IC(uncollided = True, moving = True):
    plt.ion()
    plt.figure(1)
    source_name = "plane_IC"
    print("---  ---  ---  ---  ---  ---  ---")
    print("running plane IC")
    print("---  ---  ---  ---  ---  ---  ---")
    
    main(source_name, uncollided, moving)
    plt.title("plane IC")
    plt.legend()
    plt.show(block = False)
    
def run_square_IC(uncollided = True, moving = True):
    plt.ion()
    plt.figure(2)
    source_name = "square_IC"
    print("---  ---  ---  ---  ---  ---  ---")
    print("running square IC")
    print("---  ---  ---  ---  ---  ---  ---")
    main(source_name, uncollided, moving)
    plt.title("square IC")
    plt.legend()
    plt.show(block = False)
    
def run_square_source(uncollided = True, moving = True):
    plt.ion()
    plt.figure(3)
    source_name = "square_source"
    print("---  ---  ---  ---  ---  ---  ---")
    print("running square source")
    print("---  ---  ---  ---  ---  ---  ---")
    main(source_name, uncollided, moving)
    plt.title("square source")
    plt.legend()
    plt.show(block = False)
    
def run_gaussian_IC(uncollided = True, moving = True):
    plt.ion()
    plt.figure(4)
    source_name = "gaussian_IC"
    print("---  ---  ---  ---  ---  ---  ---")
    print("running Gaussian IC")
    print("---  ---  ---  ---  ---  ---  ---")
    main(source_name, uncollided, moving)
    plt.title("Gaussian IC")
    plt.legend()
    plt.show(block = False)
    
def run_gaussian_source(uncollided = True, moving = True):
    plt.ion()
    plt.figure(5)
    source_name = "gaussian_source"
    print("---  ---  ---  ---  ---  ---  ---")
    print("running Gaussian source")
    print("---  ---  ---  ---  ---  ---  ---")
    main(source_name, uncollided, moving)
    plt.title("Gaussian source")
    plt.legend()
    plt.show(block = False)
    
def run_MMS(uncollided = False, moving = True):
    plt.ion()
    plt.figure(6)
    source_name = "MMS"
    print("---  ---  ---  ---  ---  ---  ---")
    print("running MMS problem")
    print("---  ---  ---  ---  ---  ---  ---")
    main(source_name, uncollided, moving)
    plt.title("MMS")
    plt.legend()
    plt.show(block = False)
    
def run_all():
    run_plane_IC(True, True)
    run_plane_IC(True, False)
    # run_plane_IC(False, True)        # this doesn't converge
    run_plane_IC(False, False)
    
    run_square_IC(True, True)
    run_square_IC(True, False)
    run_square_IC(False, True)
    run_square_IC(False, False)
    
    run_square_source(True, True)
    run_square_source(True, False)
    run_square_source(False, True)
    run_square_source(False, False)
    
    run_gaussian_IC(True, True)
    run_gaussian_IC(True, False)
    run_gaussian_IC(False, True)
    run_gaussian_IC(False, False)
    
    run_gaussian_source(True, True)
    run_gaussian_source(True, False)
    run_gaussian_source(False, True)
    run_gaussian_source(False, False)
    
    run_MMS(False, True)            # only one case is possible for the MMS
    
     
def main(source_name = "plane_IC", uncollided = True, moving = True):
    
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
    
    if source_type[0] == 1:
        xsb2 = np.linspace(0, tfinal , 1000)
        xsb = np.concatenate((xsb2, np.ones(1)*1.0000001))
        bench = np.concatenate((benchmark(xsb2),np.zeros(1)))
    else:
        xsb = np.linspace(0, tfinal + x0, 100000)
        bench = benchmark(xsb)
        

    
    print("uncollided  = ", uncollided)
    print("moving mesh = ", moving)
    print("---  ---  ---  ---  ---  ---  ---")
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
            print("---  ---  ---  ---  ---  ---  ---")


            plt.xlabel("x")
            plt.ylabel("scalar flux")
            
            plt.plot(xs, phi, "-o", label = f"{N_space} spaces", mfc = "none")
            
            
    saving.save_RMS(RMS_list, N_spaces, N_angles, r_times)
    plt.plot(xsb, bench, "k-", label = "benchmark")
    plt.plot(-xsb, bench, "k-")
# run_plane()       
        
    
    
    
    

            





