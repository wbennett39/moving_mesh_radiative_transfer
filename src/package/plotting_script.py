#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 15:40:39 2022

@author: bennett
"""

import matplotlib.pyplot as plt
import h5py 


class rms_plotter:
    def __init__(self):
        self.name_of_run_data_file = 'run_data_RMS.h5'
        self.case_list = ["uncol_mov", "no_uncol_stat", "uncol_stat", "no_uncol_stat"]
    
    
    def load_RMS_data(self, Ms = [2,4,6], t = 1, fign = 1, case_ns = [1,1,1,1]):
        f = h5py.File('run_data_RMS.h5', 'r')
        for M in Ms:
            if M == 2:
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
        
        f.close()
    