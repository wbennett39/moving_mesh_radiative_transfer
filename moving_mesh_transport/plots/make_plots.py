#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 15:40:39 2022

@author: bennett
"""

import matplotlib.pyplot as plt
import h5py 
from pathlib import Path


class rms_plotter:
    
    def __init__(self, tfinal, M, source_name):
        data_folder = Path("moving_mesh_transport")
        self.data_file_path = data_folder / 'run_data_RMS.h5'
        self.plot_file_path = data_folder / "plots"
        # self.case_list = ["uncol_mov", "no_uncol_stat", "uncol_stat", "no_uncol_stat"]
        self.tfinal = tfinal
        self.M = M
        self.source_type_list  = ["plane_IC", "square_IC", "square_s", "gaussian_IC", "MMS", "gaussian_s"]
        self.source_name = source_name
        
    def load_RMS_data(self, uncollided = True, moving = True):

        self.uncollided = uncollided 
        self.moving = moving
        
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
            self.clr = "b"
            self.mfc = "b"
        elif self.uncollided == False:
            self.clr = "r"
            self.mfc = "none"
            
        f = h5py.File(self.data_file_path, 'r')
        # f = h5py.File("run_data_RMS.h5", 'r')
        
        self.dest_str = str(self.source_name + "/" + "t="  + str(self.tfinal) + "/" + "RMS")
        data_str = self.uncollided * ("uncollided_")  + (not(self.uncollided))  * ("no_uncollided_")  + self.moving * ("moving_") + (not(self.moving)) * ("static_") + "M_" + str(self.M)
        
        data = f[self.dest_str + '/' + data_str]
        self.cells = data[0]
        self.RMS = data[1]
        self.angles = data[2]
        self.times = data[3]
        f.close()
    
    def plot_RMS_vs_cells(self, fign = 1, clear = False):
        plt.ion()
        plt.figure(fign)
        if clear == True:
            plt.clf()
        plt.xlabel("cells")
        plt.ylabel("RMSE")
        plt.title(f"{self.source_name} t = {self.tfinal}")
        plt.loglog(self.cells, self.RMS, self.line_mkr + self.mkr, c = self.clr, mfc = self.mfc)
        plt.savefig(self.plot_file_path / f"{self.source_name}_t={self.tfinal}_RMSE_vs_cells.pdf")
        plt.show(block = False)
    


def plot_all_RMS(tfinal, M):
    source_type_list  = ["plane_IC", "square_IC", "square_s", "gaussian_IC", "gaussian_s"]
    case_list_1 = [True, True, False, False]
    case_list_2 = [True, False, True, False]
    
    for count1, source in enumerate(source_type_list):
        plotter = rms_plotter(tfinal, M, source)
        for count2, uncollided in enumerate(case_list_1):
            moving = case_list_2[count2]
            plotter.load_RMS_data(uncollided, moving)
            if count2 == 0:
                clear = False
            if source == "plane_IC" and uncollided == False and moving == True:
                print("skipping no-uncol moving case for plane IC")
            else:
                plotter.plot_RMS_vs_cells(count1+1, clear)
    plotter = rms_plotter(tfinal, M, "MMS")
    plotter.load_RMS_data(uncollided = False, moving = True)
    plotter.plot_RMS_vs_cells(6, clear)
        
# plot_all_RMS(1, 2)
    
    
    
    
    


