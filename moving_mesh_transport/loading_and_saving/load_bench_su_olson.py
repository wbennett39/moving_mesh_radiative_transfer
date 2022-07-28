#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 12 13:31:28 2022

@author: bennett
"""

"""
Created on Thu Feb  3 15:12:33 2022

@author: bennett
"""

import h5py 
import numpy as np
from scipy.interpolate import interp1d

from pathlib import Path

###############################################################################

class load_bench_su_olson:
    def __init__(self, tfinal, x0):
        data_folder = Path("moving_mesh_transport/benchmarks")
        self.ask_for_bench = True
        self.tfinal = tfinal
        su_olson_1 = np.loadtxt(data_folder / 'su_olson_1.txt')
        self.x0 = x0
        if tfinal == 1.0:
            self.t_string_index = 0
            self.e_sol = su_olson_1
        else:
            self.ask_for_bench = False
            
    def __call__(self):
        if self.ask_for_bench == True:
                return [self.e_sol[:,0], self.e_sol[:,1], self.e_sol[:,2]]
        else:
            print("no solution file")

        
        
        
        
    
    

