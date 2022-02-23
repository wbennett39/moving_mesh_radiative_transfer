#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Feb  3 15:12:33 2022

@author: bennett
"""

import h5py 
import numpy as np
from scipy.interpolate import interp1d

###############################################################################

class load_bench:
    def __init__(self, source_type, tfinal, x0):
        self.ask_for_bench = True
        self.source_type = source_type
        self.tfinal = tfinal
        f = h5py.File("benchmarks_plane.hdf5", "r")
        self.source_type_str = ["plane_IC", "square_IC", "square_source", "gaussian_IC", "MMS"]
        self.t_eval_str = ["t = 1", " t = 5", "t = 10"]
        index_of_source_name = np.argmin(np.abs(np.array(self.source_type)-1))
        source_name = self.source_type_str[index_of_source_name]
        self.x0 = x0
        if source_name == "MMS":
            self.ask_for_bench = False
        if tfinal == 1.0:
            self.t_string_index = 0
        elif tfinal == 5.0:
            self.t_string_index = 1
        elif tfinal == 10.0:
            self.t_string_index = 2
        else:
            self.ask_for_bench = False
            
        if self.ask_for_bench == True:
            tstring = self.t_eval_str[self.t_string_index]
            solution = f[source_name][tstring]
            self.interpolated_solution = interp1d(solution[0], solution[1], kind = "cubic")
            
        f.close()
    def __call__(self, xs):
        if self.ask_for_bench == True:
            return self.interpolated_solution(xs)
        elif self.ask_for_bench == False and self.source_type[4] == 1:
            return np.exp(-xs*xs/2)/(1+self.tfinal) * np.heaviside(self.tfinal - np.abs(xs) + self.x0, 1)
        else:
            return xs * 0
        
        
        
    
    

