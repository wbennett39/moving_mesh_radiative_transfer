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
    def __init__(self, tfinal, M, source_type, moving, uncollided):
        data_folder = Path("moving_mesh_transport")
        self.config_file_path = data_folder / 'run_data_RMS.h5'
        self.M = M
        self.tfinal = int(tfinal)
        self.moving = moving
        self.uncollided = uncollided
        source_name_list = ["plane_IC", "square_IC", "square_s", "gaussian_IC", "MMS", "gaussian_s"]
        index_of_source_name = np.argmin(np.abs(np.array(source_type)-1))
        self.source_name = source_name_list[index_of_source_name]   
        if self.M == 2:
            self.mkr = "o"
        elif self.M == 4:
            self.mkr = "^"
        elif self.M == 6:
            self.mkr = "s"
        else:
            self.mkr = "p"
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
        
        
        
    def save_RMS(self, RMS_list, N_spaces, N_angles, r_times):
        print("###############")
        print("saving", self.source_name )
        print("###############")
        if self.tfinal == any(self.tlist):
            f = h5py.File(self.config_file_path, 'a')
            dest_str = str(self.source_name + "/" + "t="  + str(self.tfinal) + "/" + "RMS")
            print(dest_str)
            destination = f[dest_str]
            rms_str = self.uncollided * ("uncollided_")  + (not(self.uncollided))  * ("no_uncollided_")  + self.moving * ("moving_") + (not(self.moving)) * ("static_") + "M_" + str(self.M)
            if destination.__contains__(rms_str):
                del destination[rms_str]
            dset = destination.create_dataset(rms_str, (4, len(N_spaces)) )
            dset[0] = N_spaces
            dset[1] = RMS_list
            dset[2] = N_angles
            dset[3] = r_times
            f.close()
        
            
            
            
            
    

            
            
