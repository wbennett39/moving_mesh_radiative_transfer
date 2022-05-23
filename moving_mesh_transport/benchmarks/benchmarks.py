#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 11:58:28 2022

@author: bennett
"""
import numpy as np
import matplotlib.pyplot as plt

from .benchmark_functions import make_benchmark_file_structure, write_to_file
from .benchmark_functions import check_gaussian_tail

from .uncollided import uncollided_class
from .collided import collided_class

from ..main_functions import plot_p1_su_olson_mathematica
from .test_benchmarks import test_P1_against_mathematica

###############################################################################

class make_benchmark:
    def __init__(self, source_type, x0, t0):
        self.x0 = x0
        self.t0 = t0
        self.source_type = source_type
        
        self.call_uncollided = uncollided_class(source_type, self.x0, self.t0)
        self.call_collided = collided_class(source_type, self.x0, self.t0)
    
    def integrate(self, t, npnts):
        self.t = t
        print("t = ", t)
        self.npnts = npnts
        self.xs = np.linspace(0, t + self.x0, npnts)
        if self.source_type == "gaussian_IC_2D":
            self.xs = np.linspace(0.0, t + 1/self.x0, npnts)
            
        elif self.source_type == "gaussian_IC":
            if t == 1:
                self.xs = np.linspace(0.0, 3.41, npnts)
            elif t == 5:
                self.xs = np.linspace(0.0, 7, npnts)
            elif t == 10:
                self.xs = np.linspace(0.0, 11.6, npnts)
        elif self.source_type == "gaussian_source":
            if t == 1:
                self.xs = np.linspace(0.0, 3.5, npnts)
            elif t == 5:
                self.xs = np.linspace(0.0, 6.8, npnts)
            elif t == 10:
                self.xs = np.linspace(0.0, 11.3, npnts)
            
        self.uncollided_sol = self.call_uncollided(self.xs, t)
        self.collided_sol = self.call_collided(self.xs, t)
        
        if self.source_type == "gaussian_IC" or self.source_type == "gaussian_source" or self.source_type == "gaussian_IC_2D":
            tol = 1e-12
            index_of_zero_phi = check_gaussian_tail(self.uncollided_sol + self.collided_sol, tol)
            print(f"solution goes to {tol} at", self.xs[index_of_zero_phi])
        
        if self.source_type == "P1_su_olson_rad":
            print(self.collided_sol, "collided solution")
        
    def save(self):
        phi = self.uncollided_sol + self.collided_sol
        write_to_file(self.xs, phi, self.uncollided_sol, self.t, self.source_type, self.npnts)
    
    def clear_file(self):
        make_benchmark_file_structure()
        
    def plot(self, fign):
        plt.ion()
        plt.figure(fign)
        plt.plot(self.xs, self.uncollided_sol, "--k")
        plt.plot(self.xs, self.uncollided_sol + self.collided_sol, "-k")
        if self.source_type == "P1_su_olson_rad" and self.t == 1:
                test_P1_against_mathematica(self.t, self.xs, self.collided_sol, "rad")
        elif self.source_type == "P1_su_olson_mat" and self.t == 1:
            test_P1_against_mathematica(self.t, self.xs, self.collided_sol, "mat")
            
        plt.show()
    
        
        
    
        
        
        
        
        
        
    
    