#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 07:35:06 2022

@author: bennett
"""
import numpy as np
import h5py
from pathlib import Path

class save_output:
    def __init__(self, tfinal, N_spaces, Ms, source_type, moving, uncollided, major,
                 thermal_couple, temp_function, c, sigma, x0, cv_const, problem_type, N_angles, epsilon):
        data_folder = Path("moving_mesh_transport/local_run_data")
        self.solution_file_path = data_folder / 'run_data.hdf5'
        self.wavepoints_file_path = data_folder / 'wavepoints_crc.hdf5'
        self.problem_type = problem_type              
        self.Ms = Ms
        self.tfinal = tfinal
        self.moving = moving
        self.uncollided = uncollided
        self.major = major
        self.N_spaces = N_spaces
        self.x0 = x0
        self.thermal_couple = thermal_couple
        self.temp_function = temp_function
        self.c = c
        self.sigma = sigma
        self.cv_const = cv_const
        self.N_angles = N_angles
        self.epsilon = epsilon

        if self.problem_type == 'transport':
            self.config_file_path = data_folder / 'run_data_transport_RMS.h5'
        elif self.problem_type in ['su_olson', 'su_olson_s2', 'su_olson_thick', 'su_olson_thick_s2']:
            self.config_file_path = data_folder / 'run_data_su_olson_RMS.h5'
        else:
            print('no problem selected')

        source_name_list = ["plane_IC", "square_IC", "square_s", "gaussian_IC", "MMS", "gaussian_s"]
        gaussian_sources = ["gaussian_IC", "gaussian_s"]
        index_of_source_name = np.argmin(np.abs(np.array(source_type)-1))
        self.source_name = source_name_list[index_of_source_name]  
              
    def save_RMS(self, RMS_list, energy_RMS_list, N_angles, r_times):
        if (self.tfinal in [1,5,10,1.25, 0.8333333333333334]) or (self.thermal_couple == 1 and self.tfinal == 31.6228) :
            saving_condition = True
        else:
            saving_condition = False
            
        if self.major == 'cells':
            if saving_condition == True:
                print("saving")
                f = h5py.File(self.config_file_path, 'r+')
                if self.tfinal == 1.25:
                    tf = str(1.25)
                elif self.tfinal == 0.8333333333333334:
                    tf = str(0.8333333333333334)
                else:
                    tf = str(int(self.tfinal))
                dest_str = str(self.source_name + "/" + "t="  + tf + "/" + "RMS")
                print(dest_str)
                if not f.__contains__(dest_str):
                    f.create_group(dest_str)
                destination = f[dest_str]
                rms_str = self.uncollided * ("uncollided_")  + (not(self.uncollided))  * ("no_uncollided_")  + self.moving * ("moving_") + (not(self.moving)) * ("static_") + "M_" + str(self.Ms[0])
                if destination.__contains__(rms_str):
                    del destination[rms_str]
                dset = destination.create_dataset(rms_str, (6, len(self.N_spaces)) )
                dset[0] = self.N_spaces
                dset[1] = RMS_list
                dset[2] = N_angles
                dset[3] = r_times
                dset[4] = self.Ms
                dset[5] = energy_RMS_list
                f.close()
                
        elif self.major == 'Ms':
            if saving_condition == True :
                print("saving")
                f = h5py.File(self.config_file_path, 'r+')
                dest_str = str(self.source_name + "/" + "t="  + str(int(self.tfinal)) + "/" + "RMS")
                
                if not f.__contains__(dest_str):
                    f.create_group(dest_str)
                destination = f[dest_str]
                
                destination = f[dest_str]
                rms_str = 'Ms_' + self.uncollided * ("uncollided_")  + (not(self.uncollided))  * ("no_uncollided_")  + self.moving * ("moving") + (not(self.moving)) * ("static")
                if destination.__contains__(rms_str):
                    del destination[rms_str]
                dset = destination.create_dataset(rms_str, (6, len(self.N_spaces)) )
                dset[0] = self.Ms
                dset[1] = RMS_list
                dset[2] = N_angles
                dset[3] = r_times
                dset[4] = self.N_spaces
                dset[5] = energy_RMS_list
                f.close()
            
    def save_RMS_P1_su_olson(self, RMS_list, energy_RMS_list, N_angles, r_times, ang):
        saving_condition = True
        if int(self.sigma) == 300:
            self.source_name = 'gaussian_s_thick' 
        elif self.x0[0] == 400.0:
            self.source_name = 'su_olson_thick'
        print("saving source name", self.source_name)
        if saving_condition == True :
            if self.major == 'cells':
                print("saving RMSE")
                f = h5py.File(self.config_file_path, 'r+')
                if self.tfinal != 31.6228:
                    self.tfinal = int(self.tfinal)
                dest_str = str(self.source_name + "/" + "t="  + str(int(self.tfinal)) + "/" + "RMS" + "/" + f"S{ang}")
                print(dest_str)
                if not f.__contains__(dest_str):
                    f.create_group(dest_str)
                destination = f[dest_str]
                
                rms_str =  self.uncollided * ("uncollided_")  + (not(self.uncollided))  * ("no_uncollided_")  + self.moving * ("moving_") + (not(self.moving)) * ("static_") + "M_" + str(self.Ms[0])
                if destination.__contains__(rms_str):
                    del destination[rms_str]
                dset = destination.create_dataset(rms_str, (6, len(self.N_spaces)) )
                dset[0] = self.N_spaces
                dset[1] = RMS_list
                dset[2] = N_angles
                dset[3] = r_times
                dset[4] = self.Ms
                dset[5] = energy_RMS_list
                f.close()
            
            elif self.major == 'Ms':
                print("saving RMSE")
                f = h5py.File(self.config_file_path, 'r+')
                dest_str = str(self.source_name + "/" + "t="  + str(int(self.tfinal)) + "/" + "RMS" + "/" + f"S{ang}")
                print(dest_str)
                
                if not f.__contains__(dest_str):
                    f.create_group(dest_str)
                destination = f[dest_str]
                
                destination = f[dest_str]
                rms_str = 'Ms_' + self.uncollided * ("uncollided_")  + (not(self.uncollided))  * ("no_uncollided_")  + self.moving * ("moving") + (not(self.moving)) * ("static")
                if destination.__contains__(rms_str):
                    del destination[rms_str]
                dset = destination.create_dataset(rms_str, (6, len(self.N_spaces)))
                dset[0] = self.Ms
                dset[1] = RMS_list
                dset[2] = N_angles
                dset[3] = r_times
                dset[4] = self.N_spaces
                dset[5] = energy_RMS_list
                f.close()
                
    def save_solution(self, xs, phi, e, sol_matrix, edges,  x0_or_sigma, ws, N_space, s2, psi, epsilon, mus):
        print("saving solution")
        "transport or transfer/source_name/t = {tfinal}/c = {c}/ x0(or sigma) = {val}"
        epsilon = self.epsilon
        f = h5py.File(self.solution_file_path, 'r+')
        
        if self.problem_type == 'transport':
            folder_name = "transport"
        elif self.problem_type == 'rad_transfer_constant_cv':
            folder_name =  f"transfer_const_cv={self.cv_const}"
        elif self.problem_type == 'rad_transfer_constant_cv_s2':
            folder_name =  f"transfer_const_cv={self.cv_const}_s2"
        elif self.problem_type == 'rad_transfer_constant_cv_thick':
            folder_name =  f"transfer_const_cv={self.cv_const}_thick"
        elif self.problem_type == 'rad_transfer_constant_cv_thick_s2':
            folder_name =  f"transfer_const_cv={self.cv_const}_thick_s2"
        elif self.problem_type == "su_olson_s2":
            folder_name = "su_olson_s2"
        elif self.problem_type == "su_olson_thick_s2":
            folder_name = "su_olson_thick_s2"
        elif self.problem_type == "su_olson_thick":
            folder_name = "su_olson_thick"
        elif self.problem_type == "su_olson":
            folder_name = "su_olson"
        else:
            folder_name = 'none'
            print('not saving')
            print(self.problem_type)
        if not f.__contains__(folder_name):
            f.create_group(folder_name)
        full_str = ''

        full_str += "/" + str(self.source_name) + '_uncollided_' * (self.uncollided) + 'moving_mesh_' * (self.moving) + 'N_space = ' + str(N_space) + '_t = ' + str(self.tfinal) + '_c = ' + str(self.c) + '_x0_or_sigma = ' + str(x0_or_sigma)
        if epsilon != 1.0:
            full_str += '_epsilon=' + str(epsilon)

        if f[folder_name].__contains__('solution/' + full_str):
            del f[folder_name]['solution/'+full_str]
        print('###  ###  ###  ###  ###  ###  ###  ###  ###')
        print(folder_name + 'solution/' + full_str)
        print('###  ###  ###  ###  ###  ###  ###  ###  ###')
        dset = f[folder_name].create_dataset('solution/' + full_str, (4, len(xs)))
        
        if not f[folder_name].__contains__(full_str):
            f.create_group(full_str)



        if f[folder_name][full_str].__contains__('psi'):
            del f[folder_name][full_str]['psi']
        if f[folder_name][full_str].__contains__('mus'):
            del f[folder_name][full_str]['mus']
        print("saving solution")
        dset[0] = xs
        dset[1] = phi
        if folder_name != 'transport':
            dset[2] = e
        
        dsetpsi = f[folder_name][full_str].create_dataset('psi', data = psi)
        dsetmus = f[folder_name][full_str].create_dataset('mus', data = mus)
        # dsetpsi = psi
        # dset[3] = psi
        # dset[2] = self.ws
        # dset[3] = sol_matrix


        
        if f[folder_name].__contains__('coefficients/' + full_str):
            del f[folder_name]['coefficients/' + full_str]
        
        size = np.shape(sol_matrix)
        dset2  = f[folder_name].create_dataset('coefficients/' + full_str, data = sol_matrix)
        
        if f[folder_name].__contains__('weights/' + full_str):
            del f[folder_name]['weights/' + full_str]
            
        dset3  = f[folder_name].create_dataset('weights/' + full_str, data = ws) 

        if f[folder_name].__contains__('edges/' + full_str):
            del f[folder_name]['edges/' + full_str]
            
        dset4  = f[folder_name].create_dataset('edges/' + full_str, data = edges) 
     
        
        f.close()        
     

    def save_wave_loc(self, tpnts, leftpnts, rightpnts, T_wave_points):
    
        f = h5py.File(self.wavepoints_file_path, 'r+')


        if self.problem_type == 'transport':
            folder_name = "transport"
        elif self.problem_type == 'rad_transfer_constant_cv':
            folder_name =  f"transfer_const_cv={self.cv_const}"
        elif self.problem_type == 'rad_transfer_constant_cv_thick':
            folder_name =  f"transfer_const_cv={self.cv_const}_thick"
        elif self.problem_type == "su_olson_s2":
            folder_name = "su_olson_s2"
        elif self.problem_type == "su_olson_thick_s2":
            folder_name = "su_olson_thick_s2"
        elif self.problem_type == "su_olson_thick":
            folder_name = "su_olson_thick"
        elif self.problem_type == "su_olson":
            folder_name = "su_olson"
        else:
            print('not saving correctly')
            folder_name = 'none'


        if not f.__contains__(folder_name):
            f.create_group(folder_name)

        full_str = str(self.source_name) + 't = ' + str((self.tfinal))

        if f[folder_name].__contains__(full_str):
            del f[folder_name][full_str]


        if f[folder_name].__contains__('tpnts_' + full_str):
              del f[folder_name]['tpnts_' + full_str]
        
        if f[folder_name].__contains__('left_' + full_str):
              del f[folder_name]['left_' + full_str]
        if f[folder_name].__contains__('right_' + full_str):
              del f[folder_name]['right_' + full_str]
        if f[folder_name].__contains__('T_wave_' + full_str):
              del f[folder_name]['T_wave_' + full_str]

        dset1  = f[folder_name].create_dataset('tpnts_' + full_str, data = tpnts)
        dset2  = f[folder_name].create_dataset('left_' + full_str, data = leftpnts)
        dset3  = f[folder_name].create_dataset('right_' + full_str, data = rightpnts)
        dset3  = f[folder_name].create_dataset('T_wave_' + full_str, data = T_wave_points)



        f.close()
   
        
    
        
        
        
        
        
            
        
        
        
        
        
            
            
            
            
            
    

            
            
