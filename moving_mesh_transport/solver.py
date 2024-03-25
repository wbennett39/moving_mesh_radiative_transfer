import numpy as np

import matplotlib.pyplot as plt

from .solver_classes.functions import  convergence
from .loading_and_saving.save_output import save_output
from .loading_and_saving.load_bench import load_bench

from .solver_functions.main_functions import parameter_function
# from .main_functions import plot_p1_su_olson_mathematica

from .loading_and_saving.load_parameters import parameter_load_class
from .loading_and_saving.load_solution import load_sol

from .solver_functions.main_functions import solve, s2_source_type_selector
from .solver_functions.main_functions import plot_edges, x0_function

from .plots.plot_functions.show import show
import math
###############################################################################

"""
NEW TO DO:



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
                             self.x0, self.cv0, self.problem_type, self.N_angles, self.epsilon)
        if self.benchmarking == True:
            if self.bench_type == 'full':
                benchmark = load_bench(self.source_type, self.tfinal, self.x0_or_sigma, self.scattering_ratio, self.c_scaling)
            
            elif self.bench_type == 'S2':
                s2_source_res = s2_source_type_selector(self.sigma, self.x0[0], self.thermal_couple, self.source_type, self.weights)
                benchmark = load_bench(s2_source_res[0], self.tfinal, self.x0_or_sigma, self.scattering_ratio, self.c_scaling)
                benchmark_mat = load_bench(s2_source_res[1], self.tfinal, self.x0_or_sigma, self.scattering_ratio, self.c_scaling)
                
        print("---  ---  ---  ---  ---  ---  ---")
        print("tfinal = ", self.tfinal )
        print("c = ", self.scattering_ratio)
        print('source strength', self.source_strength)
        print("x0s", self.x0)
        print('sigma', self.sigma_t)
        print("---  ---  ---  ---  ---  ---  ---")
        if self.benchmarking == True:
            print("verifying with benchmark solution")

        #benchmark for plotting
        if self.benchmarking == True: # make this a separate function
            if self.source_type[0] == 1:
                xsb2 = np.linspace(0, self.tfinal , 1000)
                xsb = np.concatenate((xsb2, np.ones(1)*1.0000001))
                # bench = np.concatenate((benchmark(xsb2)[0],np.zeros(1)))
                bench = self.scattering_ratio*math.exp(-(1-self.scattering_ratio)*self.tfinal)*benchmark(xsb*self.scattering_ratio)[0] 
            else:
                xsb = np.linspace(0, self.tfinal + self.x0[0], 100000)
                # bench = benchmark(xsb)[0]
                bench = self.scattering_ratio*math.exp(-(1-self.scattering_ratio)*self.tfinal)*benchmark(xsb*self.scattering_ratio)[0] 
                
            
        print("uncollided  = ", uncollided)
        print("moving mesh = ", moving)
        print("---  ---  ---  ---  ---  ---  ---")
        if self.move_type[1] == 1:
            if (self.thick == True and self.source_type[2] == 1 and self.move_type[1] == 1): # make this a separate function
                sol_loader = load_sol(self.problem_type, 'square_s','transfer', self.scattering_ratio, self.N_angles[0]==2, self.cv0)
                sol_loader.call_wavepoints(self.tfinal)
                self.tpnts_wave = sol_loader.tpnts
                self.left_wave = sol_loader.left
                self.right_wave = sol_loader.right
                self.T_wave = sol_loader.T_wave
                self.wave_loc_array = np.array([[(self.tpnts_wave)], [(self.left_wave)], [(self.right_wave)], [(self.T_wave)]])
                print(self.wave_loc_array)

            elif (self.thick == True and self.source_type[5] == 1):
                sol_loader = load_sol(self.problem_type, 'gaussian_s','transfer', self.scattering_ratio, self.N_angles[0]==2, self.cv0)
                sol_loader.call_wavepoints(self.tfinal)
                self.tpnts_wave = sol_loader.tpnts
                self.left_wave = sol_loader.left
                self.right_wave = sol_loader.right
                self.T_wave = sol_loader.T_wave
                self.wave_loc_array = np.array([[(self.tpnts_wave)], [(self.left_wave)], [(self.right_wave)], [(self.T_wave)]])
        else:
            self.wave_loc_array = np.zeros((1,1,1))

    
        
        for nr in range(self.N_runs):
            for count, N_ang in enumerate(self.N_angles):
                N_space, M = parameter_function(self.major, self.N_spaces, self.Ms, count)
                
                x0_new = x0_function(self.x0, self.source_type, count)
                    
                print("M = ", M)
                print(N_space, "cells")
                print(N_ang, "angles")
                print("---  ---  ---  ---  ---  ---  ---")
    
                # sigma_t = np.ones(N_space)*self.sigma_t
                # sigma_s = np.ones(N_space)*self.scattering_ratio*self.sigma_t
                
                if self.choose_xs == True:
                    choose_xs = True
                    specified_xs = self.specified_xs
                else:
                    choose_xs = False
                    specified_xs = 0.0
                print(self.finite_domain, 'finite domain')
                xs, phi, psi, exit_dist, exit_phi, e, time, sol_matrix, angles, ws, edges, wavespeed_array, tpnts, left_edges, right_edges, wave_tpnts, wave_xpnts, T_front_location, mus = solve(self.tfinal,N_space, N_ang, M, x0_new, self.t0, self.sigma_t, 
                self.sigma_s, self.t_nodes, self.source_type, uncollided, moving, self.move_type,
                self.thermal_couple,self.temp_function, self.rt, self.at, self.e_initial, choose_xs, specified_xs, 
                self.weights, self.sigma, self.particle_v, self.edge_v, self.cv0, self.estimate_wavespeed, self.find_wave_loc, 
                self.thick, self.mxstp, self.wave_loc_array, self.find_edges_tol, self.source_strength, self.move_factor, 
                self.integrator, self.l, self.save_wave_loc, self.pad, self.leader_pad, self.xs_quad, self.eval_times, self.eval_array,
                self.boundary_on, self.boundary_source_strength, self.boundary_source, self.sigma_func, self.Msigma, self.finite_domain,
                self.domain_width, self.fake_sedov_v0, self.test_dimensional_rhs, self.epsilon, self.geometry)
                print(edges, 'final edges')
                # print(edges, "edges")
                print(wave_tpnts, wave_xpnts, "wave points")
                
                # self.xs_out = xs
                # self.phi_out = phi
                # self.psi_out = psi
                # print(psi, psi)
                

                # if self.sigma_t == 800:
                #     print('re-scaling thick solution')
                #     xs = xs * self.sigma_t
                # plt.figure(1)
                # plt.plot(xs, phi, 'o', mfc = 'none')
                if self.sigma == 0:
                    x0_or_sigma = self.x0[0]
                else:
                    x0_or_sigma = self.sigma
                    
                if self.saving == True:                                         # saves phi and coefficient matrix
                    if N_ang == 2:
                        s2 = True
                    else:
                        s2 = False
                    
                    if self.eval_times ==False:
                        saving.save_solution(xs, phi, e, sol_matrix, edges, x0_or_sigma, ws, N_space, s2 , psi, self.epsilon, mus )
                    else:
                        for it, tt in enumerate(self.eval_array):
                            saving = save_output(tt, self.N_spaces, self.Ms, self.source_type, 
                            moving, uncollided, self.major, self.thermal_couple, 
                            self.temp_function, self.scattering_ratio, self.sigma,
                            self.x0, self.cv0, self.problem_type, self.N_angles, self.epsilon)

                            saving.save_solution(xs[it], phi[it], e, sol_matrix, edges, x0_or_sigma, ws, N_space, s2, psi[it, :, :], self.epsilon, mus)
                
                
                self.r_times[count] += (time)/self.N_runs
                

                ##################################################################
                plt.figure(1)
                if self.eval_times == False:
                    plt.plot(xs, phi, "-o", label = f"{N_space} spatial cells", mfc = "none")
                    if self.benchmarking == True:
                        plt.plot(xs, benchmark(np.abs(xs))[0], '-k')

                else:
                    plt.plot(xs[-1], phi[-1,:], "-o", label = f"{N_space} spatial cells", mfc = "none")
                    plt.plot(xs[0], phi[0,:], "-o", label = f"{N_space} spatial cells", mfc = "none")
                    if self.benchmarking == True:
                        plt.plot(xs, benchmark(np.abs(xs))[0], '-k')

                plt.xlabel("x")
                plt.ylabel("scalar flux")
                if count == len(self.N_angles)-1:
                    plot_edges(edges, 1)
                # if self.thermal_couple == 1:
                #     plt.plot(xs, e, "-^", label = "energy density", mfc = "none")
                # if self.temp_function[0]==1:
                #         plt.plot(xs, np.power(e,0.25), '-s', mfc = 'none', label = 'T')
                # elif self.temp_function[1] == 1:
                #         plt.plot(xs, e/self.cv0, '-s', mfc = 'none', label = 'T')
                #         plt.plot(xs, (e/self.cv0)**4, '-s', mfc = 'none', label = 'T^4')
                plt.legend()

                plt.show()

                
                # # if count == len(self.N_angles)-1:
                # #     plot_edges(edges, 3)
                # if self.thermal_couple == 1:
                #     plt.plot(xs, e, "-^", label = "energy density", mfc = "none")
                #     if self.temp_function[0]==1:
                #         plt.plot(xs, np.power(e,0.25), '-s', mfc = 'none', label = 'T')
                #     elif self.temp_function[1] == 1:
                #         plt.plot(xs, e/self.cv0* 0.0137225, '-s', mfc = 'none', label = 'T')
                # if self.thick == True and self.sigma_t ==1 and (self.source_type[1] == 1 or self.source_type[2] == 1) :
                #     plt.xlim(self.x0[0] - self.x0[0]/8, edges[-1])
                # plt.legend()
                # plt.show()
                if count == len(self.N_angles)-1:
                    self.xs = xs
                    self.phi = phi
                    self.e = e
                    self.psi = psi
                    self.exit_dist = exit_dist
                    self.ws = ws
                    self.angles = angles
                    self.exit_phi = exit_phi
                    # for it, t in enumerate(self.eval_array):
                    #     self.exit_phi[it, 0] = np.sum(np.multiply(self.ws, self.exit_dist[it, :, 0])) 
                    #     self.exit_phi[it, 1] = np.sum(np.multiply(self.ws, self.exit_dist[it, :, 1])) 
                   
                ##################################################################
                # if self.save_wave_loc == True:
                #     plt.figure(7)
                #     plt.plot(wave_tpnts[1:], wave_xpnts[1:], label = f'{N_space} spaces')
                #     plt.legend()
                #     plt.show()
                ##################################################################
                    
                if self.benchmarking == True:
                    if self.thermal_couple == 0:
                        if self.c_scaling == False:
                            benchmark_solution = benchmark(np.abs(xs))[0] #benchmark for RMS
                        elif self.c_scaling == True:
                            benchmark_solution = self.scattering_ratio*math.exp(-(1-self.scattering_ratio)*self.tfinal)*benchmark(np.abs(xs*self.scattering_ratio))[0] #benchmark for RMS

                        RMS = np.sqrt(np.mean((phi - benchmark_solution)**2))
                    
                    elif self.thermal_couple == 1:
                        if  (self.weights == "gauss_lobatto" and (self.sigma != 300 and self.x0[0] != 400)):
                            e_xs = benchmark(np.abs(xs))[2][:,0]
                            phi_bench = benchmark(np.abs(xs))[2][:,1]
                            e_bench = benchmark(np.abs(xs))[2][:,2]
                            
                      
                            # ##################################################################
                            # plt.figure(3)
                            # plt.plot(e_xs, phi_bench, "-k")
                            # plt.plot(-e_xs, phi_bench, "-k")
                            # plt.plot(e_xs, e_bench, "--k")
                            # plt.plot(-e_xs, e_bench, "--k")
                   

                            # plt.show()
                            ##################################################################
                            
                        elif self.weights == "gauss_legendre" or self.sigma == 300 or self.x0[0] == 400:
                            print("loading s2 bench")
                            phi_bench = benchmark(np.abs(xs))[0]
                            e_bench = benchmark_mat(np.abs(xs))[0]
                                # plot_p1_su_olson_mathematica()
                                            
                        RMS = np.sqrt(np.mean((phi - phi_bench)**2))
                        RMS_energy = np.sqrt(np.mean((e - e_bench)**2))
                        self.RMS_list_energy[count] = RMS_energy
                    
                    self.RMS_list[count] = RMS

                ##################################################################
                if self.find_wave_loc == True:
                    print('saving ')
                    saving.save_wave_loc(tpnts, left_edges, right_edges, T_front_location)

                    # plt.figure(2)
                    # plt.plot(tpnts, wavespeed_array, '-o', label = "calculated wavespeed")
                    # plt.plot(tpnts,  2/np.sqrt(tpnts + 1e-12), label = "2/sqrt(t)")
                    # plt.plot(tpnts,  8/np.sqrt(tpnts + 1e-12), label = "8/sqrt(t)")
                    # plt.plot(tpnts,  16/np.sqrt(tpnts + 1e-12), label = "16/sqrt(t)")
                    # plt.plot(tpnts,  32/np.sqrt(tpnts + 1e-12), label = "32/sqrt(t)")
                    # plt.plot(tpnts,  64/np.sqrt(tpnts + 1e-12), label = "64/sqrt(t)")
                    # plt.plot(tpnts,  128/np.sqrt(tpnts + 1e-12), label = "128/sqrt(t)")
                    # plt.plot(tpnts,  256/np.sqrt(tpnts + 1e-12), label = "256/sqrt(t)")
                    # # plt.plot(tpnts,  4/np.sqrt(tpnts)/math.sqrt(3))
                    # # plt.ylim(0,wavespeed_array[0])
                    # plt.legend()
                    # plt.show()


                    # plt.figure(5)
                    # plt.plot(tpnts, left_edges, '-o', label = 'left_edge')
                    # plt.plot(tpnts, right_edges, '-o', label = 'right_edge')
                    plt.plot(tpnts, T_front_location, '-o', label = 'wave temperature front')
                    plt.plot(tpnts, np.ones(tpnts.size)* xs[-1], 'k-')

                    plt.legend()
                    plt.show()
                ##################################################################
                
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
    
    
    
        # if self.benchmarking == True:
            # if self.bench_type == 'S2':
            #     saving.save_RMS_P1_su_olson(self.RMS_list, self.RMS_list_energy, self.N_angles, self.r_times, self.N_angles[0])
            # else:
            #     saving.save_RMS(self.RMS_list, self.RMS_list_energy, self.N_angles, self.r_times)
            
            # if ((self.benchmarking == True) and self.thermal_couple == 0):
                # plt.figure(3)
                # plt.plot(xsb, bench, "k-")#, label = "benchmark")
                # plt.plot(-xsb, bench, "k-")
                # plt.show()
                
            # elif (self.weights == "gauss_legendre" or self.sigma == 300 or self.x0[0] == 400) and self.thermal_couple == 1:
                # plt.figure(1)
                # phi_bench_plot = benchmark(np.abs(xsb))[0]
                # e_bench_plot = benchmark_mat(np.abs(xsb))[0]
                
                # plt.plot(xsb, phi_bench_plot, "-k")
                # plt.plot(-xsb, phi_bench_plot, "-k")
                # plt.plot(xsb, e_bench_plot, "--k")
                # plt.plot(-xsb, e_bench_plot, "--k")
                # if self.x0[0] == 400:
                #     plt.xlim(self.x0[0] -50, xs[-1])

                # plt.show()
                # if int(self.x0[0]) == 400:
                #     plt.xlim(self.x0[0]-10, self.x0[0] + 4 * math.sqrt(self.tfinal) /math.sqrt(3) + 10)
    
                    
                


        
        




        
    
    
    
    

            





