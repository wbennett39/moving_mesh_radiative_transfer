import numpy as np
import scipy.integrate as integrate
import quadpy
import matplotlib.pyplot as plt
from timeit import default_timer as timer
import yaml
from pathlib import Path

# from .solver_classes.build_problem import build
# from .solver_classes.matrices import G_L
# from .solver_classes.numerical_flux import LU_surf
# from .solver_classes.sources import source_class
# from .solver_classes.uncollided_solutions import uncollided_solution
# from .solver_classes.phi_class import scalar_flux
# from .solver_classes.mesh import mesh_class
# from .solver_classes.rhs_class import rhs_class
# from .solver_classes.make_phi import make_output
# from .solver_classes.radiative_transfer import T_function
from .solver_classes.functions import find_nodes, convergence
from .save_output import save_output
from .load_bench import load_bench

from .main_functions import parameter_function
# from .main_functions import plot_p1_su_olson_mathematica

from .load_parameters import parameter_load_class

from .main_functions import solve
###############################################################################
""" 
to do:

[] update README


[] either s or source, pick one
[] save solution somewhere
[x] variable c 
[x] plot all problems on same graph
[x] combine Ms and cells with solve function
[x] make a parameter function to merge major ms and major cells 
[x] clean up constants in rhs
[] clean loading and saving scripts
[] fix s2 RMS plots
[] get source type array under control 
[] variable sigma 


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
# ###############################################################################
# data_folder = Path("moving_mesh_transport")
# config_file_path = data_folder / "config.yaml"
# ###############################################################################



class main_class(parameter_load_class):
    
    pass

    
    
#     with open(config_file_path, 'r') as file:
#        parameters = yaml.safe_load(file)
#     file.close()
#     #'C:/users/bennett/Documents/GitHub/MovingMesh/moving_mesh_radiative_transfer/src/package/
#     tfinal = float(parameters['all']['tfinal'])
#     N_spaces = np.array(parameters['all']['N_spaces'])
#     Ms = np.array(parameters['all']['Ms'])
#     N_runs = int(parameters['all']['N_runs'])
#     t_nodes = int(parameters['all']['tnodes'])
#     rt = float(parameters['all']['rt'])
#     at = float(parameters['all']['at'])
#     t0 = float(parameters['all']['t0'])
#     scattering_ratio = float(parameters['all']['c'])
#     major = str(parameters['all']['major'])
#     thermal_couple = int(parameters['all']['radiative_transfer'])
#     temp_function = np.array(parameters['all']['temperature_dependence'])
#     e_initial = float(parameters['all']['e_initial'])
#     weights = str(parameters['all']['weights'])
#     benchmarking = int(parameters['all']['benchmarking'])
#     saving = int(parameters['all']['save_solution'])
    
    
#     N_angles = np.array(parameters[source_name]['N_angles'])
#     x0 = np.array(parameters[source_name]['x0'])
#     source_type = np.array(parameters[source_name]['source_type'])
#     move_type = np.array(parameters[source_name]['move_type'])
    

#     r_times = np.zeros(len(N_angles))
#     RMS_list = np.zeros(len(N_angles))
#     RMS_list_energy = np.zeros(len(N_angles))
    
    
#     if source_name in ["gaussian_IC", "gaussian_source"]:
#         sigma = float(parameters[source_name]['sigma'])
#         x0_or_sigma = sigma
#     else:
#         x0_or_sigma = x0
#         sigma = 0
    def main(self, uncollided = True, moving = True):  
        saving = save_output(self.tfinal, self.N_spaces, self.Ms, self.source_type, 
                             moving, uncollided, self.major, self.thermal_couple, 
                             self.temp_function, self.scattering_ratio, self.sigma)
        
        if (self.weights != "gauss_legendre") and self.sigma != 300:
            benchmark = load_bench(self.source_type, self.tfinal, self.x0_or_sigma)
            
        elif (self.weights != "gauss_legendre") and self.sigma == 300:
            benchmark = load_bench([0,0,0,0,0,0,0,0,0,0,0,0,1,0], self.tfinal, self.x0_or_sigma)
            benchmark_mat = load_bench([0,0,0,0,0,0,0,0,0,0,0,0,0,1], self.tfinal, self.x0_or_sigma)
        
        elif (self.weights == "gauss_legendre") and (self.thermal_couple == 1):
            if self.source_type[2] == 1:
                    benchmark = load_bench([0,0,0,0,0,0,0,0,0,1,0,0,0,0], self.tfinal, self.x0_or_sigma)
                    benchmark_mat = load_bench([0,0,0,0,0,0,0,0,0,0,1,0,0,0], self.tfinal, self.x0_or_sigma)
            elif self.source_type[5] == 1 and self.sigma == 0.5 :
                 benchmark = load_bench([0,0,0,0,0,0,0,0,0,0,0,1,0,0], self.tfinal, self.x0_or_sigma)
                 benchmark_mat = load_bench([0,0,0,0,0,0,0,0,0,0,0,0,1,0], self.tfinal, self.x0_or_sigma)
            elif self.source_type[5] == 1 and self.sigma == 300 :
                 benchmark = load_bench([0,0,0,0,0,0,0,0,0,0,0,0,1,0], self.tfinal, self.x0_or_sigma)
                 benchmark_mat = load_bench([0,0,0,0,0,0,0,0,0,0,0,0,0,1], self.tfinal, self.x0_or_sigma)
                
        
        print("--- --- --- --- --- --- ")
        print("tfinal = ", self.tfinal )
        print("c = ", self.scattering_ratio)
        print("x0s", self.x0)
        print("--- --- --- --- --- --- ")
        if self.benchmarking == True:
            print("verifying with benchmark solution")
        
        #benchmark for plotting
        if self.benchmarking == True:
            if self.source_type[0] == 1:
                xsb2 = np.linspace(0, self.tfinal , 1000)
                xsb = np.concatenate((xsb2, np.ones(1)*1.0000001))
                bench = np.concatenate((benchmark(xsb2)[0],np.zeros(1)))
            else:
                xsb = np.linspace(0, self.tfinal + self.x0[0], 100000)
                bench = benchmark(xsb)[0] 
            
        print("uncollided  = ", uncollided)
        print("moving mesh = ", moving)
        print("---  ---  ---  ---  ---  ---  ---")
    
        
        for nr in range(self.N_runs):
            for count, N_ang in enumerate(self.N_angles):
                N_space, M = parameter_function(self.major, self.N_spaces, self.Ms, count)
                
                if self.source_type[3] or self.source_type[4] == 1:
                    x0_new = self.x0[count]
                    print("x0 = ", x0_new)
                else:
                    x0_new = self.x0[0]
                    
                print("M = ", M)
                print(N_space, "cells")
                print(N_ang, "angles")
                print("---  ---  ---  ---  ---  ---  ---")
    
                sigma_t = np.ones(N_space)
                sigma_s = np.ones(N_space)
                
                if self.benchmarking == True and self.thermal_couple == 1 and self.weights == "gauss_lobatto" and self.sigma != 300:
                    choose_xs = True
                    specified_xs = benchmark(np.abs(np.linspace(0, self.tfinal + self.x0)))[2][:,0]
                else:
                    choose_xs = False
                    specified_xs = 0.0
                    
                xs, phi, e, time = solve(self.tfinal, N_space, N_ang, M, x0_new, 
                                         self.t0, sigma_t, sigma_s, self.t_nodes, 
                                         self.scattering_ratio, self.source_type, 
                                         uncollided, moving, self.move_type, self.thermal_couple,
                                         self.temp_function, self.rt, self.at, self.e_initial, 
                                         choose_xs, specified_xs, self.weights, self.sigma)
                
                
                self.r_times[count] += (time)/self.N_runs
                
                plt.plot(xs, phi, "-o", label = f"{N_space} spaces", mfc = "none")
                
                if self.thermal_couple == 1:
                    plt.plot(xs, e, "-^", label = "energy density", mfc = "none")
                    
                plt.xlabel("x")
                plt.ylabel("scalar flux")
                    
                if self.benchmarking == True:
                    if self.thermal_couple == 0:
                        benchmark_solution = benchmark(np.abs(xs))[0] #benchmark for RMS
                    
                        RMS = np.sqrt(np.mean((phi - benchmark_solution)**2))
                    
                    elif self.thermal_couple == 1:
                        if  self.weights == "gauss_lobatto" and self.sigma != 300:
                            e_xs = benchmark(np.abs(xs))[2][:,0]
                            phi_bench = benchmark(np.abs(xs))[2][:,1]
                            e_bench = benchmark(np.abs(xs))[2][:,2]
        
                            plt.plot(e_xs,phi_bench, "-k")
                            plt.plot(-e_xs,phi_bench, "-k")
                            plt.plot(e_xs,e_bench, "--k")
                            plt.plot(-e_xs,e_bench, "--k")
        
                        elif self.weights == "gauss_legendre" or self.sigma ==300:
                            phi_bench = benchmark(np.abs(xs))[0]
                            e_bench = benchmark_mat(np.abs(xs))[0]
                            
                            if count == len(self.N_angles)-1:
                                plt.plot(xs, phi_bench, "-k")
                                plt.plot(-xs,phi_bench, "-k")
                                plt.plot(xs,e_bench, "--k")
                                plt.plot(-xs,e_bench, "--k")
                                # plot_p1_su_olson_mathematica()
                                            
                        RMS = np.sqrt(np.mean((phi - phi_bench)**2))
                        RMS_energy = np.sqrt(np.mean((e - e_bench)**2))
                        self.RMS_list_energy[count] = RMS_energy
                        
                    self.RMS_list[count] = RMS
                
                print(N_space, "spaces", "    ", "%.4f" % (time), "time elapsed")
                
                if self.benchmarking == True:
                    print("RMSE", RMS)
                    if self.major == 'cells':
                        if self.thermal_couple == 1:
                            print("energy RMSE", RMS_energy)
                            if count > 0:
                                print("material energy order", "%.2f" % convergence(self.RMS_list_energy[count-1], self.N_spaces[count-1], RMS_energy, N_space))
                        if count > 0:
                            print("Order", "%.2f" % convergence(self.RMS_list[count-1], self.N_spaces[count-1], RMS, N_space))
    
                
                print("---  ---  ---  ---  ---  ---  ---")
    
    
    
        if self.benchmarking == True:
            if (self.weights == "gauss_legendre") and (self.thermal_couple == 1):
                saving.save_RMS_P1_su_olson(self.RMS_list, self.RMS_list_energy, self.N_angles, self.r_times)
            else:
                saving.save_RMS(self.RMS_list, self.RMS_list_energy, self.N_angles, self.r_times)
                
            if ((self.tfinal == 1 or self.tfinal == 5 or self.tfinal == 10) and self.thermal_couple == 0):
                plt.plot(xsb, bench, "k-", label = "benchmark")
                plt.plot(-xsb, bench, "k-")

        
        




        
    
    
    
    

            





