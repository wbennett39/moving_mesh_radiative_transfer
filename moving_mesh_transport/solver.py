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
from .solver_classes.radiative_transfer import T_function
from .solver_classes.functions import find_nodes, convergence
from .save_output import save_output
from .load_bench import load_bench
###############################################################################
""" 
to do:

[] update README

[x] no-uncollided plane source doesnt converge -- take it out

[] either s or source, pick one
[] save solution somewhere
[x] variable c 
[] plot all problems on same graph
[] combine Ms and cells with solve function


ideas for tests:
    source - check if integrator returns analytic integral of plane uncollided for x inside t 
    G, and L, check against integrated normPns -
    check normPns

long term goals:
[] comments for all classes, functions
[] fix np dot contiguous warning
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
    # run_plane_IC(True, True)
    # run_plane_IC(True, False)
    # # # run_plane_IC(False, True)        # this doesn't converge
    # run_plane_IC(False, False)
    
    # run_square_IC(True, True)
    # run_square_IC(True, False)
    # run_square_IC(False, True)
    # run_square_IC(False, False)
    
    # run_square_source(True, True)
    # run_square_source(True, False)
    # run_square_source(False, True)
    # run_square_source(False, False)
    
    run_gaussian_IC(True, True)
    run_gaussian_IC(True, False)
    run_gaussian_IC(False, True)
    run_gaussian_IC(False, False)
    
    run_gaussian_source(True, True)
    run_gaussian_source(True, False)
    run_gaussian_source(False, True)
    run_gaussian_source(False, False)
    
    # run_MMS(False, True)            # only one case is possible for the MMS
    
     
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
    t0 = float(parameters['all']['t0'])
    scattering_ratio = float(parameters['all']['c'])
    major = str(parameters['all']['major'])
    thermal_couple = int(parameters['all']['radiative_transfer'])
    temp_function = np.array(parameters['all']['temperature_dependence'])
    
    N_angles = np.array(parameters[source_name]['N_angles'])
    x0 = float(parameters[source_name]['x0'])
    source_type = np.array(parameters[source_name]['source_type'])
    move_type = np.array(parameters[source_name]['move_type'])
    
    
    r_times = np.zeros(len(N_angles))
    RMS_list = np.zeros(len(N_angles))
    x0s = np.ones(4)*x0
    
    saving = save_output(tfinal, N_spaces, Ms, source_type, moving, uncollided, major, thermal_couple, temp_function)
    benchmark = load_bench(source_type, tfinal, x0)
    
    
    print("--- --- --- --- --- --- ")
    print("tfinal = ", tfinal )
    print("--- --- --- --- --- --- ")
    if source_type[0] == 1:
        xsb2 = np.linspace(0, tfinal , 1000)
        xsb = np.concatenate((xsb2, np.ones(1)*1.0000001))
        bench = np.concatenate((benchmark(xsb2)[0],np.zeros(1)))
    else:
        xsb = np.linspace(0, tfinal + x0, 100000)
        bench = benchmark(xsb)[0]
        

    
    print("uncollided  = ", uncollided)
    print("moving mesh = ", moving)
    print("---  ---  ---  ---  ---  ---  ---")
    for nr in range(N_runs):
        if major == 'cells':
            for count, N_space in enumerate(N_spaces):
                if source_type[3] == 1 and N_space >= 32:
                    x0 += 1
                sigma_t = np.ones(N_space)
                sigma_s = np.ones(N_space)
                M = Ms[0]
                print("M = ", M)
                N_ang = N_angles[count]
                xs, phi, time = solve(tfinal, N_space, N_ang, M, x0, t0, sigma_t, sigma_s, t_nodes, scattering_ratio, source_type, 
                          uncollided, moving, move_type, thermal_couple, temp_function, rt, at)
                
                
                r_times[count] += (time)/N_runs
                benchmark_solution = benchmark(np.abs(xs))[0]
                RMS = np.sqrt(np.mean((phi - benchmark_solution)**2))
                RMS_list[count] = RMS
                print(N_space, "spaces", "    ", "%.4f" % (time), "time elapsed")
                print("RMSE", RMS)
                if count > 0:
                    print("Order", "%.2f" % convergence(RMS_list[count-1], N_spaces[count-1], RMS, N_space))
                plt.xlabel("x")
                plt.ylabel("scalar flux")
                
                plt.plot(xs, phi, "-o", label = f"{N_space} spaces", mfc = "none")
                   
            
        elif major == 'Ms':
            for count, M in enumerate(Ms):
                N_space = N_spaces[1]
                print("cells", N_space)
                if source_type[3] == 1 and N_space >= 32:
                    x0 += 1
                sigma_t = np.ones(N_space)
                sigma_s = np.ones(N_space)
                N_ang = N_angles[count]
               
                xs, phi, time = solve(tfinal, N_space, N_ang, M, x0, t0, sigma_t, sigma_s, t_nodes, scattering_ratio, source_type, 
                          uncollided, moving, move_type, thermal_couple, temp_function, rt, at)
               
                r_times[count] += (time)/N_runs

                
                benchmark_solution = benchmark(np.abs(xs))[0]
                RMS = np.sqrt(np.mean((phi - benchmark_solution)**2))
                RMS_list[count] = RMS
                
                if major == 'Ms':
                    print(M, "M", "    ", "%.4f" % (time), "time elapsed")
                    print("RMSE", RMS)
                # if count >0:
                #     print("Order", "%.2f" % convergence(RMS_list[count-1], Ms[count-1], RMS, N_space))
                
                
                
                
                # phi = make_phi(N_ang, ws, xs, sol_last, M, edges) 
                print("---  ---  ---  ---  ---  ---  ---")
    
    
                plt.xlabel("x")
                plt.ylabel("scalar flux")
                
                plt.plot(xs, phi, "-o", label = f"{N_space} spaces", mfc = "none")
            
            
    saving.save_RMS(RMS_list, N_angles, r_times)
    plt.plot(xsb, bench, "k-", label = "benchmark")
    plt.plot(-xsb, bench, "k-")
# run_plane()       
        


def solve(tfinal, N_space, N_ang, M, x0, t0, sigma_t, sigma_s, t_nodes, scattering_ratio, source_type, 
          uncollided, moving, move_type, thermal_couple, temp_function, rt, at):
    mus = quadpy.c1.gauss_lobatto(N_ang).points
    ws = quadpy.c1.gauss_lobatto(N_ang).weights
    xs_quad = quadpy.c1.gauss_legendre(M+2).points
    ws_quad = quadpy.c1.gauss_legendre(M+2).weights
    t_quad = quadpy.c1.gauss_legendre(t_nodes).points
    t_ws = quadpy.c1.gauss_legendre(t_nodes).weights
    initialize = build(N_ang, N_space, M, tfinal, x0, t0, scattering_ratio, mus, ws, xs_quad,
                       ws_quad, sigma_t, sigma_s, source_type, uncollided, moving, move_type, t_quad, t_ws,
                       thermal_couple, temp_function)
    initialize.make_IC()
    IC = initialize.IC
    if thermal_couple == 0:
        deg_freedom = N_ang*N_space*(M+1)
    elif thermal_couple == 1:
        deg_freedom = (N_ang+1)*N_space*(M+1)
    mesh = mesh_class(N_space, x0, tfinal, moving, move_type) 
    matrices = G_L(initialize)
    num_flux = LU_surf(initialize)
    source = source_class(initialize)
    uncollided_sol = uncollided_solution(initialize)
    flux = scalar_flux(initialize)
    rhs = rhs_class(initialize)
    transfer = T_function(initialize)
    
    def RHS(t, V):
        return rhs.call(t, V, mesh, matrices, num_flux, source, uncollided_sol, flux, transfer)
    
    start = timer()
    sol = integrate.solve_ivp(RHS, [0.0,tfinal], IC.reshape(deg_freedom), method='DOP853', t_eval = [tfinal], rtol = rt, atol = at)
    end = timer()
    sol_last = sol.y[:,-1].reshape((N_ang,N_space,M+1))
    mesh.move(tfinal)
    edges = mesh.edges
    
    xs = find_nodes(edges, M)
    output = make_output(tfinal, N_ang, ws, xs, sol_last, M, edges, uncollided)
    phi = output.make_phi(uncollided_sol)
    
    computation_time = start-end
    
    return xs, phi, computation_time
    
    
    
    

            





