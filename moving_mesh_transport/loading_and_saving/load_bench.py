#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Feb  3 15:12:33 2022

@author: bennett
"""
import math
import h5py 
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

from pathlib import Path

###############################################################################

class load_bench:
    def __init__(self, source_type, tfinal, x0_or_sigma, c, c_scaling):
        data_folder = Path("moving_mesh_transport/benchmarks")
        benchmark_file_path = data_folder / "benchmarks.hdf5"

        self.c = c
        self.ask_for_bench = True
        self.source_type = source_type
        print(tfinal, 'tfinal')
        self.tfinal = float(tfinal)
        self.c_scaling = c_scaling
        if self.c_scaling == True:
            self.tfinal = self.tfinal * self.c
            print(self.tfinal, 'tfinal for benchmark')

        su_olson = 0
        if self.tfinal == 1.0:
            su_olson = np.loadtxt(data_folder / 'su_olson_1.txt')
        elif self.tfinal == 31.6228:
            su_olson = np.loadtxt(data_folder / 'su_olson_31.txt')
        
        f = h5py.File(benchmark_file_path, "r+")

        
        self.source_type_str = ["plane_IC", "square_IC", "square_source", 
                                "gaussian_IC", "MMS", "gaussian_source",
                                "gaussian_IC_2D", "line_source", "P1_su_olson_rad",
                                "P1_su_olson_mat", "P1_gaussian_rad", "P1_gaussian_mat", 
                                "P1_gaussian_rad_thick", "P1_gaussian_mat_thick",
                                "P1_su_olson_rad_thick", "P1_su_olson_mat_thick"]

        self.t_eval_str = ["t = 1", "t = 5", "t = 10", "t = 31.6228"]
        index_of_source_name = np.argmin(np.abs(np.array(self.source_type)-1))
        source_name = self.source_type_str[index_of_source_name]
        self.x0 = x0_or_sigma
        if source_name == "MMS":
            # self.ask_for_bench = False
            self.xs = np.linspace(0, tfinal + 1/10)
        
            
        if self.ask_for_bench == True:
            print("loading bench for ", source_name)
            # tstring = self.t_eval_str[self.t_string_index]
            if tfinal != 31.6228:
                tstring = "t = " + str((self.tfinal))
            elif tfinal == 31.6228:
                tstring = "t = " + str(float(tfinal))
                
            
            self.solution_dataset = f[source_name][tstring]
                
            self.xs = self.solution_dataset[0]
            
            self.phi = self.solution_dataset[1]
            self.phi_u = self.solution_dataset[2]
            
            self.scale_transport_benchmark()
            
            self.interpolated_solution = interp1d(self.xs, self.phi, kind = "cubic")
            self.interpolated_uncollided_solution = interp1d(self.xs, self.phi_u, kind = "cubic")
        
        self.e_sol = su_olson
            
        f.close()
        
    def stich_solution(self, xs):
        """ if an answer is requested outside of the solution interval, adds 
            zeros to the edge of the interpolated solution
        """
        stiched_solution = xs*0
        stiched_uncollided_solution = xs*0
        original_xs = self.xs
        

        if abs(xs[-1]-xs[0]) <= 1e-12:
            symmetric = True
        else:
            symmetric = False
        edge_of_interpolated_sol = original_xs[-1]
        if symmetric == False:
            edge_index = np.argmin(np.abs(xs-edge_of_interpolated_sol))
            self.xs_inside = xs[0:edge_index]
            stiched_solution[0:edge_index] = self.interpolated_solution(self.xs_inside)
            stiched_uncollided_solution[0:edge_index] = self.interpolated_uncollided_solution(self.xs_inside)
        elif symmetric == True:
            # assuming that the vector is even 
            middle_index = int(xs.size/2)
            edge_index = np.argmin(np.abs(xs[middle_index:-1]-edge_of_interpolated_sol))
            self.xs_inside = xs[middle_index-edge_index:edge_index+middle_index]
            stiched_solution[middle_index-edge_index:edge_index+middle_index] = self.interpolated_solution(self.xs_inside)
            stiched_uncollided_solution[middle_index-edge_index:edge_index+middle_index] = self.interpolated_uncollided_solution(self.xs_inside)
            
        return [stiched_solution, stiched_uncollided_solution]
    
    def scale_transport_benchmark(self):
        if self.c_scaling == True and 1 == 0:
            phi_new = self.phi*0
            x_new = self.xs/self.c
            phi_new = self.c*math.exp(-(1-self.c)*self.tfinal/self.c)*self.phi
            # self.xs = x_new
            self.phi = phi_new

    
    def __call__(self, xs):
        if self.ask_for_bench == True:
            original_xs = self.xs

            if xs[-1] > original_xs[-1]:
                beyond_solution_domain = True
            else:
                beyond_solution_domain = False
            if (beyond_solution_domain == False):
                return [self.interpolated_solution(xs), self.interpolated_uncollided_solution(xs), self.e_sol]
            elif (beyond_solution_domain == True):
                return [self.stich_solution(xs)[0], self.stich_solution(xs)[1], self.e_sol]

        elif self.ask_for_bench == False and self.source_type[4] == 1:
            return [np.exp(-xs*xs/2)/(1+self.tfinal) * np.heaviside(self.tfinal - np.abs(xs) + self.x0, 1),0,0]
        else:
            return xs * 0
        
        
        
        
    
    

