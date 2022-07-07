#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 08:38:18 2022

@author: bennett
"""

import h5py 
from pathlib import Path


class load_sol:
    def __init__(self, source_name, rad_or_transfer, c, s2):

        data_folder = Path("moving_mesh_transport")
        self.data_file_path = data_folder / 'run_data.h5'
        self.source_name = source_name
        self.rad_or_transfer = rad_or_transfer
        self.c = c 
        self.s2 = s2
        
    
    def call_sol(self, tfinal, M, x0_or_sigma, N_space, mat_or_rad):
        full_str = self.rad_or_transfer
        if self.s2 == True:
            full_str += '_s2'
        full_str += '/N_space = ' + str(N_space) + '_t = ' + str(int(tfinal)) + '_c = ' + str(self.c) + '_x0_or_sigma = ' + str(x0_or_sigma)
        f = h5py.File(self.data_file_path, "r+")
        sol_data = f[full_str]
        self.xs = sol_data[0]
        self.phi = sol_data[1]
        self.e = sol_data[2]
        
        self.ws = f['weights/' + full_str][:]


        if self.rad_or_transfer == 'transfer':
            if mat_or_rad == 'rad':
                self.coeff_mat = f['coefficients/' + full_str][:-1,:,:]
                
            elif mat_or_rad =='mat':
                self.coeff_mat = f['coefficients/' + full_str][-1,:,:]
        else:
            self.coeff_mat = f['coefficients/' + full_str][:,:,:]

        f.close()
        
        
        
        