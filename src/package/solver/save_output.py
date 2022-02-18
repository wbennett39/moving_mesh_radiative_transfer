#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 07:35:06 2022

@author: bennett
"""
import numpy as np
import h5py

class save_output:
    def __init__(self, tfinal, M, N_ang, source_type, moving, uncollided):
        self.M = M
        self.tfinal = int(tfinal)
        self.N_ang = N_ang
        self.moving = moving
        self.uncollided = uncollided
        source_name_list = ["plane_IC", "square_IC", "square_source", "truncated_gaussian_IC"]
        index_of_source_name = np.argmin(np.abs(np.array(source_type)-1))
        self.source_name = source_name_list[index_of_source_name]
        if self.M == 2:
            self.mkr = "o"
        elif self.M == 4:
            self.mkr = "^"
        elif self.M == 6:
            self.mkr = "s"
        if self.moving == True:
            self.line_mkr = "-"
        elif self.moving == False:
            self.line_mkr = "--"
        if self.uncollided == True:
            self.clr = "c0"
        elif self.uncollided == False:
            self.clr = "c1"
        self.mkr_string = self.line_mkr + self.mkr
        
        
        
    def save_RMS(self, RMS_list, N_spaces):
        f = h5py.File('run_data_RMS.h5', 'a')
        dest_str = str(self.source_name + "/" + "t="  + str(self.tfinal) + "/" + "RMS")
        # dest_str = "plane_IC"
        destination = f[dest_str]
        rms_str = self.uncollided * ("uncollided_")  + (not(self.uncollided))  * ("no_uncollided_")  + self.moving * ("moving_") + (not(self.moving)) * ("static_")  + str(self.N_ang) + "_angles" + "_M_" + str(self.M)
        if destination.__contains__(rms_str):
            del destination[rms_str]
        dset = destination.create_dataset(rms_str, (2, len(N_spaces)) )
        dset[0] = N_spaces
        dset[1] = RMS_list
        f.close()
    
            
            
            
            
    

            
            