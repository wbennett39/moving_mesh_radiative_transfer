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
    def __init__(self, tfinal, N_spaces, Ms, source_type, moving, uncollided, major, thermal_couple, temp_function):
        data_folder = Path("moving_mesh_transport")
        if thermal_couple == 0:
            self.config_file_path = data_folder / 'run_data_RMS.h5'
        elif thermal_couple == 1:
            self.config_file_path = data_folder / 'run_data_radiative_transfer_RMS.h5'
        self.Ms = Ms
        self.tfinal = int(tfinal)
        self.moving = moving
        self.uncollided = uncollided
        self.major = major
        self.N_spaces = N_spaces
        source_name_list = ["plane_IC", "square_IC", "square_s", "gaussian_IC", "MMS", "gaussian_s"]
        index_of_source_name = np.argmin(np.abs(np.array(source_type)-1))
        self.source_name = source_name_list[index_of_source_name]  
        if self.major == 'cells':
            if self.Ms[0] == 2:
                self.mkr = "o"
            elif self.Ms[0] == 4:
                self.mkr = "^"
            elif self.Ms[0] == 6:
                self.mkr = "s"
            else:
                self.mkr = "p"
        elif self.major == 'Ms':
            self.mkr = 'D'
        if self.moving == True:
            self.line_mkr = "-"
        elif self.moving == False:
            self.line_mkr = "--"
        if self.uncollided == True:
            self.clr = "c0"
        elif self.uncollided == False:
            self.clr = "c1"
        self.mkr_string = self.line_mkr + self.mkr
        self.tlist = [1,5,10]
        
        self.thermal_couple = thermal_couple
        self.temp_function = temp_function
        
        
        
    def save_RMS(self, RMS_list, energy_RMS_list, N_angles, r_times):
        print("###############")
        print("###############")
        if (self.tfinal == 1 or self.tfinal == 5 or self.tfinal == 10) or (self.thermal_couple == 1 and self.tfinal == 31.6228) :
            saving_condition = True
        else:
            saving_condition = False

        if self.major == 'cells':
            if saving_condition == True:
                print("saving")
                f = h5py.File(self.config_file_path, 'a')
                dest_str = str(self.source_name + "/" + "t="  + str(self.tfinal) + "/" + "RMS")
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
                f = h5py.File(self.config_file_path, 'a')
                dest_str = str(self.source_name + "/" + "t="  + str(self.tfinal) + "/" + "RMS")
                destination = f[dest_str]
                rms_str = 'Ms_' + self.uncollided * ("uncollided_")  + (not(self.uncollided))  * ("no_uncollided_")  + self.moving * ("moving") + (not(self.moving)) * ("static")
                if destination.__contains__(rms_str):
                    del destination[rms_str]
                dset = destination.create_dataset(rms_str, (5, len(self.N_spaces)) )
                dset[0] = self.Ms
                dset[1] = RMS_list
                dset[2] = N_angles
                dset[3] = r_times
                dset[4] = self.N_spaces
                dset[5] = energy_RMS_list
                f.close()
            
    def save_RMS_P1_su_olson(self, RMS_list, energy_RMS_list, N_angles, r_times):
        if ((self.thermal_couple == 1) and (self.tfinal == 31.6228) or (self.tfinal == 1)):
            saving_condition = True
       
        if saving_condition == True :
            if self.major == 'cells':
                print("saving")
                f = h5py.File(self.config_file_path, 'a')
                
                dest_str = str(self.source_name + "/" + "t="  + str(self.tfinal) + "/" + "RMS" + "/" + "S2")
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
                print("saving")
                f = h5py.File(self.config_file_path, 'a')
                dest_str = str(self.source_name + "/" + "t="  + str(self.tfinal) + "/" + "RMS" + "/" + "S2")
                destination = f[dest_str]
                rms_str = 'Ms_' + self.uncollided * ("uncollided_")  + (not(self.uncollided))  * ("no_uncollided_")  + self.moving * ("moving") + (not(self.moving)) * ("static")
                if destination.__contains__(rms_str):
                    del destination[rms_str]
                dset = destination.create_dataset(rms_str, (5, len(self.N_spaces)) )
                dset[0] = self.Ms
                dset[1] = RMS_list
                dset[2] = N_angles
                dset[3] = r_times
                dset[4] = self.N_spaces
                dset[5] = energy_RMS_list
                f.close()
            
            
            
            
            
    

            
            
