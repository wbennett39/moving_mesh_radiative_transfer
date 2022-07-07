import numpy as np

import matplotlib.pyplot as plt

from .solver_classes.functions import  convergence
from .save_output import save_output
from .load_bench import load_bench

from .main_functions import parameter_function
# from .main_functions import plot_p1_su_olson_mathematica

from .load_parameters import parameter_load_class

from .main_functions import solve, s2_source_type_selector
import math
###############################################################################
""" 
to do:

[] update README


[] either s or source, pick one
[x] save solution somewhere
[x] variable c 
[x] plot all problems on same graph
[x] combine Ms and cells with solve function
[x] make a parameter function to merge major ms and major cells 
[x] clean up constants in rhs
[] clean loading and saving scripts
[x] fix s2 RMS plots
[] get source type array under control 
[x] variable sigma 
[] simplify if statements in solver
[] plot final edges


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


    def main(self, uncollided = True, moving = True):  

        
        saving = save_output(self.tfinal, self.N_spaces, self.Ms, self.source_type, 
                             moving, uncollided, self.major, self.thermal_couple, 
                             self.temp_function, self.scattering_ratio, self.sigma,
                             self.x0)
        
        if self.bench_type == 'full':
            benchmark = load_bench(self.source_type, self.tfinal, self.x0_or_sigma)
        
        elif self.bench_type == 'S2':
            s2_source_res = s2_source_type_selector(self.sigma, self.x0[0], self.thermal_couple, self.source_type, self.weights)
            benchmark = load_bench(s2_source_res[0], self.tfinal, self.x0_or_sigma)
            benchmark_mat = load_bench(s2_source_res[1], self.tfinal, self.x0_or_sigma)
                
        
        print("---  ---  ---  ---  ---  ---  ---")
        print("tfinal = ", self.tfinal )
        print("c = ", self.scattering_ratio)
        print("x0s", self.x0)
        print("---  ---  ---  ---  ---  ---  ---")
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
                
                if self.benchmarking == True and self.thermal_couple == 1 and self.weights == "gauss_lobatto" and self.sigma != 300 and self.x0 != 400:
                    choose_xs = True
                    specified_xs = benchmark(np.abs(np.linspace(0, self.tfinal + self.x0)))[2][:,0]
                else:
                    choose_xs = False
                    specified_xs = 0.0
                    
                xs, phi, e, time, sol_matrix, ws = solve(self.tfinal, N_space, N_ang, M, x0_new, 
                                         self.t0, sigma_t, sigma_s, self.t_nodes, 
                                         self.scattering_ratio, self.source_type, 
                                         uncollided, moving, self.move_type, self.thermal_couple,
                                         self.temp_function, self.rt, self.at, self.e_initial, 
                                         choose_xs, specified_xs, self.weights, self.sigma, 
                                         self.particle_v, self.edge_v)
                print(xs[0], xs[-1], 'edges')
                if self.sigma == 0:
                    x0_or_sigma = self.x0
                else:
                    x0_or_sigma = self.sigma
                    
                if self.saving == True:                                         # saves phi and coefficient matrix
                    if N_ang == 2:
                        s2 = True
                    else:
                        s2 = False
                        
                    saving.save_solution(xs, phi, e, sol_matrix, x0_or_sigma, ws, N_space, s2)
                
                
                
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
                        if  (self.weights == "gauss_lobatto" and (self.sigma != 300 and self.x0[0] != 400)):
                            e_xs = benchmark(np.abs(xs))[2][:,0]
                            phi_bench = benchmark(np.abs(xs))[2][:,1]
                            e_bench = benchmark(np.abs(xs))[2][:,2]
        
                            plt.plot(e_xs,phi_bench, "-k")
                            plt.plot(-e_xs,phi_bench, "-k")
                            plt.plot(e_xs,e_bench, "--k")
                            plt.plot(-e_xs,e_bench, "--k")
                            
                            
       
                        elif self.weights == "gauss_legendre" or self.sigma == 300 or self.x0[0] == 400:
                            print("loading s2 bench")
                            phi_bench = benchmark(np.abs(xs))[0]
                            e_bench = benchmark_mat(np.abs(xs))[0]
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
                                print("material energy convergence order", "%.2f" % convergence(self.RMS_list_energy[count-1], self.N_spaces[count-1], RMS_energy, N_space))
                        if count > 0:
                            print("radiation energy density convergence order", "%.2f" % convergence(self.RMS_list[count-1], self.N_spaces[count-1], RMS, N_space))
    
                
                print("---  ---  ---  ---  ---  ---  ---")
    
    
    
        if self.benchmarking == True:
            if self.bench_type == 'S2':
                saving.save_RMS_P1_su_olson(self.RMS_list, self.RMS_list_energy, self.N_angles, self.r_times, self.N_angles[0])
            else:
                saving.save_RMS(self.RMS_list, self.RMS_list_energy, self.N_angles, self.r_times)
                
            if ((self.tfinal == 1 or self.tfinal == 5 or self.tfinal == 10) and self.thermal_couple == 0):
                plt.plot(xsb, bench, "k-", label = "benchmark")
                plt.plot(-xsb, bench, "k-")
                plt.show()
                
            elif self.weights == "gauss_legendre" or self.sigma == 300 or self.x0[0] == 400:
                phi_bench_plot = benchmark(np.abs(xsb))[0]
                e_bench_plot = benchmark_mat(np.abs(xsb))[0]
                
                
                plt.plot(xsb, phi_bench_plot, "-k")
                plt.plot(-xsb, phi_bench_plot, "-k")
                plt.plot(xsb, e_bench_plot, "--k")
                plt.plot(-xsb, e_bench_plot, "--k")
                # if int(self.x0[0]) == 400:
                #     plt.xlim(self.x0[0]-10, self.x0[0] + 4 * math.sqrt(self.tfinal) /math.sqrt(3) + 10)
                    
                


        
        




        
    
    
    
    

            





