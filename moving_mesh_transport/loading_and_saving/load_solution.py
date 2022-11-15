#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 08:38:18 2022

@author: bennett
"""

import h5py 
from pathlib import Path
import numpy as np


class load_sol:
    def __init__(self, problem_name = 'transport', source_name = 'square_s', rad_or_transfer = 'rad', c = 1.0, s2 = False, cv0 = 0.0):

        data_folder = Path("moving_mesh_transport")

        # self.data_file_path = data_folder / 'run_data_crc.hdf5'
        self.data_file_path = data_folder / 'run_data.h5' 

        self.wavepoints_file_path = data_folder / 'wavepoints_crc.hdf5'
        print(self.wavepoints_file_path)
        self.source_name = source_name
        self.rad_or_transfer = rad_or_transfer
        self.c = c 
        self.s2 = s2
        self.problem_name = problem_name
        self.cv0 = cv0
        if self.problem_name == 'rad_transfer_const_cv_thick':
            self.problem_name = f'transfer_const_cv={self.cv0}_thick'
        elif self.problem_name == 'rad_transfer_const_cv':
            self.problem_name = f'transfer_const_cv={self.cv0}'
    
    def call_sol(self, tfinal, M, x0_or_sigma, N_space, mat_or_rad, uncollided, moving):
        # full_str = self.rad_or_transfer
        full_str = ''
        if self.s2 == True:
            full_str += '_s2'
        full_str += "/" + str(self.source_name) + '_uncollided_' * (uncollided) + 'moving_mesh_' * (moving) + 'N_space = ' + str(N_space) + '_t = ' + str(tfinal) + '_c = ' + str(self.c) + '_x0_or_sigma = ' + str(x0_or_sigma)
        f = h5py.File(self.data_file_path, "r+")

        if self.problem_name != 'su_olson_thick_s2': # FIX THIS LATER 
            # print(f.keys())

            sol_data = f[self.problem_name]['solution/'+full_str]

            self.xs = sol_data[0]
            self.phi = sol_data[1]
            self.e = sol_data[2]
      


        # print(f[self.problem_name]['coefficients/'].keys())

        self.ws = f[self.problem_name]['weights/' + full_str][:]
        if self.problem_name != 'rad_transfer_const_cv=0.03_thick':
            self.edges = f[self.problem_name]['edges/' + full_str][:]
        # print(f[self.problem_name]['weights'].keys())
        # print(f[self.problem_name]['coefficients'].keys())
        
        if self.rad_or_transfer == 'transfer':
            if mat_or_rad == 'rad':
                self.coeff_mat = f[self.problem_name]['coefficients/' + full_str][:-1,:,:]
                
            elif mat_or_rad =='mat':
                self.coeff_mat = f[self.problem_name]['coefficients/' + full_str][-1,:,:]
        else:
            self.coeff_mat = f[self.problem_name]['coefficients/' + full_str][:,:,:]

        f.close()
    
    def call_wavepoints(self, tfinal):

        f = h5py.File(self.wavepoints_file_path, "r+")

        full_str = str(self.source_name) + 't = ' + str((tfinal))    
        if self.problem_name == 'su_olson_thick':
                if not f.__contains__(self.problem_name):
                    f.create_group(self.problem_name)
                self.tpnts = f['su_olson_thick']['tpnts_' + full_str][:]
                self.left = f['su_olson_thick']['left_' + full_str][:]
                self.right = f['su_olson_thick']['right_' + full_str][:]
                self.T_wave = f[folder_name]['T_wave_' + full_str][:]

        elif self.problem_name == 'su_olson_thick_s2':
                if not f.__contains__(self.problem_name):
                    f.create_group(self.problem_name)
                if not f[self.problem_name].__contains__('tpnts_' + full_str):
                    self.tpnts = np.array([tfinal])
                    self.left = np.array([0])
                    self.right = np.array([0])
                    self.T_wave = np.array([0])
                else:
                    self.tpnts = f['su_olson_thick_s2']['tpnts_' + full_str][:]
                    self.left = f['su_olson_thick_s2']['left_' + full_str][:]
                    self.right = f['su_olson_thick_s2']['right_' + full_str][:]
                    self.T_wave = f[self.problem_name]['T_wave_' + full_str][:]




        elif self.problem_name == 'rad_transfer_constant_cv_thick':
            folder_name =  f"transfer_const_cv={self.cv0}_thick"
            if not f.__contains__(folder_name):
                f.create_group(folder_name)
            if not f[folder_name].__contains__('tpnts_' + full_str):
                self.tpnts = np.array([tfinal])
                self.left = np.array([0])
                self.right = np.array([0])
                self.T_wave = np.array([0])
            else:
                self.tpnts = f[folder_name]['tpnts_' + full_str][:]
                self.left = f[folder_name]['left_' + full_str][:]
                self.right = f[folder_name]['right_' + full_str][:]
                self.T_wave = f[folder_name]['T_wave_' + full_str][:]

        elif self.problem_name == 'rad_transfer_constant_cv_thick_s2':
            folder_name =  f"transfer_const_cv={self.cv0}_thick_s2"
            if not f.__contains__(folder_name):
                f.create_group(folder_name)

            if not f[folder_name].__contains__('tpnts_' + full_str):
                self.tpnts = np.array([tfinal])
                self.left = np.array([0])
                self.right = np.array([0])
                self.T_wave = np.array([0])
            else:
                self.tpnts = f[folder_name]['tpnts_' + full_str][:]
                self.left = f[folder_name]['left_' + full_str][:]
                self.right = f[folder_name]['right_' + full_str][:]
                self.T_wave = f[folder_name]['T_wave_' + full_str][:]


        else:
            # print(f.keys())
            print(self.problem_name, "pn")
            self.tpnts = f[self.problem_name]['tpnts_' + full_str][:]
            self.left = f[self.problem_name]['left_' + full_str][:]
            self.right = f[self.problem_name]['right_' + full_str][:]
            self.T_wave = f[folder_name]['T_wave_' + full_str][:]





        f.close()
    
        
        
        
        