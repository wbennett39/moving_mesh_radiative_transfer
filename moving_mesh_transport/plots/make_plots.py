#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 15:40:39 2022

@author: bennett
"""

import matplotlib.pyplot as plt
import math
import h5py 
from pathlib import Path
from ..loading_and_saving.load_bench import load_bench
import numpy as np
from .plot_functions.show import show
from .plot_functions.show_loglog import show_loglog
from .plot_functions.sn_labels import logplot_sn_labels,logplot_sn_labels_2
from .plot_functions.show_loglog_timeplots import show_loglog_time
from scipy.stats import linregress 
from .plot_functions.order_triangle import order_triangle
from .plot_functions.coeff_con import coeff_con
from ..loading_and_saving.load_solution import load_sol
from ..solver_classes.functions import  convergence
class rms_plotter:
    
    def __init__(self, tfinal, M, source_name, major):
        self.data_folder = Path("moving_mesh_transport")
        self.data_file_path = self.data_folder / 'run_data_transport_RMS.h5'
        
        self.plot_file_path = self.data_folder / "plots" 
        # self.case_list = ["uncol_mov", "no_uncol_stat", "uncol_stat", "no_uncol_stat"]
        self.tfinal = tfinal
        self.source_type_list  = ["plane_IC", "square_IC", "square_s", "gaussian_IC", "MMS", "gaussian_s", "su_olson"]
        self.source_name = source_name
        self.major = major 
        self.M = M
        self.su_olson_data_file_path = self.data_folder / "run_data_su_olson_RMS.h5"
        if self.M == 2:
            self.mkr = "o"
        elif self.M == 4:        
            self.mkr = "^"
        elif self.M == 6:
            self.mkr = "s"
        else:
            self.mkr = "p"
        
    def find_c(self):
        """
        Finds the intercept from log log data. Skips the first data point
        """
        x = np.log(self.cells[1:])
        y = np.log(self.RMS[1:])
        y1 = (np.log(self.energy_RMS[1:]))
        slope, intercept, r, p, se = linregress(x, y)
        slope2, intercept2, r, p, se = linregress(x, y1)
        print('intercept', np.exp(intercept))
        print('energy intercept', np.exp(intercept2))
        return np.exp(intercept)
    
    def find_c_semilog(self):
        """
        Finds the intercept from log log data. Skips the first data point
        """
        x = self.Ms[1:]
        y = np.log(self.RMS[1:])
        y1 = np.log(self.energy_RMS[1:])
        slope, intercept, r, p, se = linregress(x, y)
        slope2, intercept, r, p, se = linregress(x, y1)
        print(slope, 'slope - c1')
        print(slope2, 'energy density slope - c1')

        if abs(r) < 0.9:
            print('bad correlation')
            print(r, 'r')
    
        return slope

        
    def load_RMS_data(self, uncollided = True, moving = True):

        self.uncollided = uncollided 
        self.moving = moving
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
        
        if self.source_name =='su_olson' or self.source_name == "su_olson_energy":
            self.dest_str = str('square_s' + "/" + "t="  + str(self.tfinal) + "/" + "RMS")
            
        if self.source_name =='su_olson_s2' or self.source_name == "su_olson_energy_s2":
            self.dest_str = str('square_s' + "/" + "t="  + str(self.tfinal) + "/" + "RMS" + "/" + "S2")
            
        elif self.source_name in ['gaussian_s2', "gaussian_energy_s2"]:
            self.dest_str = str('gaussian_s' + "/" + "t="  + str(self.tfinal) + "/" + "RMS" + "/" + "S2")
            
        elif self.source_name in ['gaussian_s_thick_s2', 'gaussian_s_thick_s2_energy']:
            self.dest_str = str("gaussian_s_thick" + "/" + "t="  + str(self.tfinal) + "/" + "RMS" + "/" + "S2")
            
        elif self.source_name in ['gaussian_s_thick_s8', 'gaussian_s_thick_s8_energy']:
            self.dest_str = str('gaussian_s_thick' + "/" + "t="  + str(self.tfinal) + "/" + "RMS" + "/" + "S8")
        
        elif self.source_name in ['su_olson_thick_s2', 'su_olson_thick_s2_energy']:
            self.dest_str = str("su_olson_thick" + "/" + "t="  + str(self.tfinal) + "/" + "RMS" + "/" + "S2")
            
        elif self.source_name in ['su_olson_thick_s8', 'su_olson_thick_s8_energy' ]:
            self.dest_str = str("su_olson_thick" + "/" + "t="  + str(self.tfinal) + "/" + "RMS" + "/" + "S8")
            
        if self.major == 'cells':
            data_str = self.uncollided * ("uncollided_")  + (not(self.uncollided))  * ("no_uncollided_")  + self.moving * ("moving_") + (not(self.moving)) * ("static_") + "M_" + str(self.M)
        elif self.major == 'Ms':
            data_str = 'Ms_' + self.uncollided * ("uncollided_")  + (not(self.uncollided))  * ("no_uncollided_")  + self.moving * ("moving") + (not(self.moving)) * ("static") 
        
        
        if self.source_name not in ["su_olson", "su_olson_energy", "su_olson_energy_s2",
                                    "su_olson_s2", "gaussian_s2", "gaussian_energy_s2", 
                                    'gaussian_s_thick_s2', 'gaussian_s_thick_s8', 'su_olson_thick_s2',
                                    'su_olson_thick_s8', 'gaussian_s_thick_s2_energy', 'gaussian_s_thick_s8_energy',
                                    'su_olson_thick_s2_energy','su_olson_thick_s8_energy']:
            print(f['gaussian_IC'].keys())
            data = f[self.dest_str + '/' + data_str]
    
            if self.major == 'cells':
                self.cells = data[0]
                # self.Ms = data[4]x
                # self.M = self.Ms[0]
            elif self.major == 'Ms':
                self.Ms = data[0]
                self.cells = data[4]
                self.N_spaces = self.cells[0]
            self.RMS = data[1]
            self.angles = data[2]
            self.times = data[3]
    
            f.close()
        
        if self.source_name in ["su_olson", "su_olson_energy", "su_olson_energy_s2",
                                    "su_olson_s2", "gaussian_s2", "gaussian_energy_s2", 
                                    'gaussian_s_thick_s2', 'gaussian_s_thick_s8', 'su_olson_thick_s2',
                                    'su_olson_thick_s8', 'gaussian_s_thick_s2_energy', 'gaussian_s_thick_s8_energy',
                                    'su_olson_thick_s2_energy','su_olson_thick_s8_energy']:
            # print("loading s2 RMS data")
            f_rad = h5py.File(self.su_olson_data_file_path, 'r')
            print(f_rad.keys())
            rad_data = f_rad[self.dest_str + '/' + data_str]
            if self.major == 'cells':
                self.cells = rad_data[0]
            elif self.major == 'Ms':
                self.Ms = rad_data[0]
                self.cells = rad_data[4]
                self.N_spaces = self.cells[0]
            self.RMS = rad_data[1]
            self.angles = rad_data[2]
            self.times = rad_data[3]
            self.energy_RMS = rad_data[5]
    
    def plot_RMS_vs_cells(self, fign = 1, clear = False, xlimright = 260):
        plt.ion()
        plt.figure(fign)

        print(self.cells, 'cells')
        print(self.RMS, 'RMSE')
        xlimleft = 2.5

        if clear == True:
            plt.clf()
    
        plt.xlabel("cells", fontsize = 20)
        plt.ylabel("RMSE", fontsize = 20)
        # plt.title(f"{self.source_name} t = {self.tfinal}")
        file_path_string = str(self.plot_file_path) + '/' "RMS_plots" '/' + f"{self.source_name}_t={self.tfinal}_M={self.M}_RMSE_vs_cells"

        #######################################################################
        if self.tfinal == 1:
            xlimleft = 1.5
            
            if self.source_name == "MMS":
                plt.ylim(5e-12,10e-3)
                if self.M == 6:
                    self.cells = self.cells[:3]
                    self.RMS = self.RMS[:3]
                    intercept = self.find_c()
                    order_triangle(3, 5, 7, intercept, 2, 1.5)
                if self.M == 2:
                    self.cells = self.cells[:5]
                    self.RMS = self.RMS[:5]
                    intercept = self.find_c()
                    order_triangle(9, 15, 3, intercept, 2, 1.05)
                    
            elif self.source_name == "square_s":
                if self.uncollided == True and self.moving == True:
                    intercept = self.find_c()
                    if self.M == 6:
                        order_triangle(9, 15, 2, intercept, 2, 1.7)
                    if self.M == 4:
                        order_triangle(9, 15, 2, intercept, 2, 1.05)
                self.cells = self.cells[1:5]
                self.RMS = self.RMS[1:5]
                self.angles = self.angles[1:5]
                plt.ylim(1e-8,1e-1)
                xlimleft = 3.5
                
            elif self.source_name =='square_IC':
                plt.ylim(1e-8,1e-1)
                xlimright = 70
                if self.uncollided == True and self.moving == True:
                    intercept = self.find_c()
                    if self.M == 6:
                        order_triangle(9, 15, 2, intercept, 2, 1.5)
                    elif self.M == 4:
                        order_triangle(9, 15, 2, intercept, 2, 1.5)
                        
            elif self.source_name == 'plane_IC':
                plt.ylim(1e-6, 1)
                if self.uncollided == True and self.moving == True:
                    
                    intercept = self.find_c()
                    ######## test line ##############
                    plt.plot(self.cells, intercept * self.cells**(-1.16))
                    
                    ################################
                 
                    if self.M == 6:
                        order_triangle(9, 15, 1, intercept, 2, 1.7)
                    elif self. M == 4:
                        order_triangle(9, 15, 1, intercept, 2, 1.7)
                        
            elif self.source_name == 'gaussian_IC':
              
                xlimright = 67
                if self.M == 6:
                    if self.M == 6:
                        xlimright = 18
                        self.cells = self.cells[:4]
                        self.RMS = self.RMS[:4]
                        self.angles = self.angles[:4]
                        
                elif self.M == 3:
                    xlimright = 67
                        
                        
                if self.uncollided == True and self.moving == True:
                    intercept = self.find_c()
                    
                    if self.M == 6:
                        order_triangle(3, 5, 6, intercept, 2, 1.7)
                        plt.ylim(1e-8,1e-2)
                        
                    elif self.M == 4:
                        order_triangle(5, 9, 4, intercept, 2, 1.7)
                        plt.ylim(1e-8,1e-1)
                        
                    elif self.M == 3:
                        order_triangle(5, 12, 4, intercept, 2, 1.7)
                        plt.ylim(1e-8,1e-1)
                        
            elif self.source_name == 'gaussian_s':
                self.cells = self.cells[:]
                self.RMS = self.RMS[:]
                self.angles = self.angles[:]
                
                
                if self.M == 6:
                   self.cells = self.cells[:4]
                   self.RMS = self.RMS[:4]
                   self.angles = self.angles[:4]
                   xlimright = 18
                
                elif self.M == 3:
                    xlimright = 67
                   
                   
                if self.uncollided == True and self.moving == True:
                    intercept = self.find_c()
                    if self.M == 6:
                        order_triangle(3, 5, 6, intercept, 2, 1.3)
                        order_triangle(7, 11, 7, intercept, 2, 2.2)
                        plt.ylim(1e-8,1e-1)
                        xlimright = 18
                    if self.M == 4:
                        order_triangle(5, 9, 4, intercept, 2, 1.7)
                        plt.ylim(1e-9,1e-1)
                    elif self.M == 3:
                        order_triangle(5, 9, 4, intercept, 2, 1.7)
                        plt.ylim(1e-8,1e-1)
                        
            elif self.source_name == 'su_olson' or self.source_name == "su_olson_energy":
                xlimright = 36
                plt.ylim(1e-5, 1e-2)
                if self.uncollided == False and self.moving == False:
                    intercept = self.find_c()
                    if self.M == 4:
                        order_triangle(5, 8, 1, intercept, 2, 1.2)


            elif self.source_name in ['su_olson_s2', 'su_olson_energy_s2']:
                xlimright = 36
                intercept = self.find_c()
                if self.uncollided == True:
                    # self.cells = self.cells[:-1]
                    # self.RMS = self.RMS[:-1]
                    intercept = self.find_c()

                    if self.M == 4 and self.moving == True:
                        order_triangle(6, 11, 2, intercept, 2, 1.2)
                plt.ylim(1e-7, 1e-2)

            elif self.source_name in ["gaussian_s2", "gaussian_energy_s2"]:
                xlimright = 24
                if self.uncollided == True and self.moving == True:
                    intercept = self.find_c()
                    if self.M == 6:
                        order_triangle(6, 11, 7, intercept, 1.5, 1.01)
            
            elif self.source_name in ['su_olson_thick_s2', 'su_olson_thick_s2_energy']:
                xlimright = 260
                xlimleft = 14
                if self.uncollided == True:
                    # self.cells = self.cells[:-1]
                    # self.RMS = self.RMS[:-1]
                    intercept = self.find_c()
                    if self.M == 4 and self.moving == True:
                        order_triangle(18, 24, 2, intercept, 2, 1.2)
                # plt.ylim(1e-7, 1e-2)

            elif self.source_name in ['gaussian_s_thick_s2', 'gaussian_s_thick_s2_energy']:
                xlimright = 24
                plt.ylim(5e-10,1e-3)
                
                self.cells = self.cells[:-1]
                self.RMS = self.RMS[:-1]
                self.energy_RMS = self.energy_RMS[:-1]
                if self.uncollided == True and self.moving == True:
                    intercept = self.find_c()
                    if self.M == 6:
                        order_triangle(4, 6, 7, intercept, 2.9, 2)



            energy_list = ["su_olson_energy", 'su_olson_energy_s2', "gaussian_energy_s2", 'gaussian_s_thick_s2_energy', 
                               'gaussian_s_thick_s8_energy','su_olson_thick_s2_energy','su_olson_thick_s8_energy' ]
            
            if self.source_name not in energy_list:

                plt.loglog(self.cells, self.RMS, self.line_mkr + self.mkr, c = self.clr, mfc = self.mfc)
            
            if self.source_name in energy_list:
                plt.loglog(self.cells, self.energy_RMS, self.line_mkr + self.mkr, c = self.clr, mfc = self.mfc)
            
            if  (self.source_name not in ['gaussian_IC', 'gauissian_s', 'gaussian_s2', 'gaussian_energy_s2', 'su_olson_s2', 'su_olson_energy_s2', 
            'su_olson', 'su_olson_energy', 'gaussian_s_thick_s2', 'gaussian_s_thick_s2_energy', 'su_olson_thick_s2', 'su_olson_thick_s2_energy']  and self.source_name != 'gaussian_s') and (self.uncollided == False and self.moving == False or self.source_name == "MMS" and self.M == 4):
                logplot_sn_labels(self.cells, self.RMS, self.angles, 0.3, fign )
            
            if self.source_name in ['su_olson', 'su_olson_energy']:
                if self.uncollided == False and self.moving == False:
                    logplot_sn_labels(self.cells, self.RMS, self.angles, 0.7, fign )


            elif (self.source_name == 'gaussian_IC' or self.source_name == 'gaussian_s') and (self.uncollided == False and self.moving == False):     
                logplot_sn_labels_2(self.cells, self.RMS, self.angles, 0.5, fign )

            # plt.savefig(self.plot_file_path / "RMS_plots" / f"{self.source_name}_t={self.tfinal}_RMSE_vs_cells.pdf")
            file_path_string = str(self.plot_file_path) + '/' "RMS_plots" '/' + f"{self.source_name}_t={self.tfinal}_M={self.M}_RMSE_vs_cells"
            show_loglog(file_path_string, xlimleft, xlimright)
            print(self.source_name)
            print('uncollided', self.uncollided)
            print('moving', self.moving)
            # print('intercept', self.find_c())
            print("--- --- --- --- --- ")

        elif self.tfinal in [1.25, 0.8333333333333334]:
            # file_path_string = str(self.plot_file_path) + '/' "RMS_plots" '/' + f"{self.source_name}_t={self.tfinal}_M={self.M}_RMSE_vs_cells"

            if self.source_name == 'square_IC':
                plt.ylim(1e-8,1e-1)
                xlimright = 36
                xlimleft = 2.5
                if self.uncollided == True and self.moving == True:
                    intercept = self.find_c()
                    if self.M == 6:
                        order_triangle(9, 15, 2, intercept, 2, 1.5)
                    elif self.M == 4:
                        order_triangle(9, 15, 2, intercept, 2, 1.5)
            elif self.source_name == 'gaussian_IC':
                plt.ylim(1e-8,1e-1)
                xlimright = 36
                xlimleft = 1.
                # if self.tfinal == 0.8333333333333334:
                #     xlimright = 18
                #     self.RMS = self.RMS[:-1]
                #     self.cells = self.cells[:-1]
                #     self.angles = self.angles[:-1]
                if self.uncollided == True and self.moving == True:
                    intercept = self.find_c()
                    if self.M == 6:
                        order_triangle(9, 15, 7, intercept, 2, 1.5)
                    elif self.M == 4:
                        order_triangle(5, 8, 5, intercept, 2, 1.9)
            
            plt.loglog(self.cells, self.RMS, self.line_mkr + self.mkr, c = self.clr, mfc = self.mfc)
            if self.uncollided == False and self.moving == False:
                logplot_sn_labels(self.cells, self.RMS, self.angles, 0.1, fign )
            show_loglog(file_path_string, xlimleft, xlimright)


        elif self.tfinal == 100:
            file_path_string = str(self.plot_file_path) + '/' "RMS_plots" '/' + f"{self.source_name}_t={self.tfinal}_M={self.M}_RMSE_vs_cells"
            xlimleft = 10
            xlimright = 160
            plt.ylim(1e-6,1e-2)
            if self.source_name in ['su_olson_thick_s2', 'su_olson_thick_s2_energy']:
                self.cells = self.cells[:-1]
                self.RMS = self.RMS[:-1]
                if self.uncollided == True:
                    if self.M == 4:
                        intercept = self.find_c()
                        order_triangle(12, 16, 5, intercept, 1.5, 1.5)

            plt.loglog(self.cells, self.RMS, self.line_mkr + self.mkr, c = self.clr, mfc = self.mfc)
            show_loglog(file_path_string, xlimleft, xlimright)
            

       
            
    
            plt.show(block = False)
        #######################################################################
        elif self.tfinal == 10:
            xlimleft = 1
            xlimright = 256
            plt.loglog(self.cells, self.RMS, self.line_mkr + self.mkr, c = self.clr, mfc = self.mfc)
            
            if self.source_name == "plane_IC":
                intercept = self.find_c()
                if self.M ==6 and self.uncollided == False and self.moving == False:
                    logplot_sn_labels(self.cells, self.RMS, self.angles, 0.3, fign )
                elif self.M == 6 and self.moving== True and self.uncollided == True:
                    order_triangle(3, 6, 6, intercept, 2, 1.7)
                
            file_path_string = str(self.plot_file_path) + '/' "RMS_plots" '/' + f"{self.source_name}_t={self.tfinal}_RMSE_vs_cells"
            show_loglog(file_path_string, xlimleft, xlimright)
            print(self.source_name)
            print('uncollided', self.uncollided)
            print('moving', self.moving)
            print('intercept', self.find_c())
            print("--- --- --- --- --- ")
        
    def plot_RMS_vs_Ms(self, source_type, fign = 1, clear = False):
        xlimleft = 1
        xlimright = self.Ms[-1] + 2
        
        plt.ion()
        plt.figure(fign)
        if clear == True:
            plt.clf()
    
        plt.xlabel("M", fontsize = 20)
        plt.ylabel("RMSE", fontsize = 20)
        # plt.title(f"{self.source_name} t = {self.tfinal}")
        #######################################################################
        
        if self.source_name == 'MMS':
            self.Ms = self.Ms[:4]
            self.RMS = self.RMS[:4]
            self.angles = self.angles[:4]
            xlimright = 12
            plt.ylim = 10e-14
            self.find_c_semilog()
        energy_list = ['su_olson_thick_s2_energy', 'gaussian_s_thick_s2_energy', 'gaussian_energy_s2']

        if self.source_name not in energy_list:    
            if self.source_name == "gaussian_IC" or self.source_name == "gaussian_s":
                self.Ms = self.Ms[:]
                self.RMS = self.RMS[:]
                self.angles = self.angles[:]
                plt.ylim = 10e-10
                if (self.uncollided == False) and (self.moving == False):
                    logplot_sn_labels_2(self.Ms, self.RMS, self.angles, 0.02, fign )
                self.find_c_semilog()
                
            elif self.source_name == 'square_IC' or self.source_name == 'square_s':
                if (self.uncollided == False) and (self.moving == False):
                    logplot_sn_labels(self.cells, self.RMS, self.angles, 0.1, fign )
            
            elif self.source_name in ['gaussian_s2','gaussian_energy_s2', 'gaussian_s_thick_s2']:
                self.find_c_semilog()
                if self.source_name in ['gaussian_s_thick_s2']:
                    if self.tfinal == 1:
                        plt.ylim(1e-8, 1e-1)
                    elif self.tfinal == 10:
                        plt.ylim(1e-6, 1e0)
                    elif self.tfinal == 100:
                        plt.ylim(1e-6,1e0)

            elif self.source_name == "su_olson_thick_s2":
                if (self.uncollided == False) and (self.moving == False):
                    if self.tfinal == 100:
                        plt.ylim(5e-6,1e0)
                
                    self.find_c_semilog()

               

            plt.semilogy(self.Ms, self.RMS, self.line_mkr + self.mkr, c = self.clr, mfc = self.mfc)


        elif self.source_name in energy_list:
            print('material energy density')
            if self.source_name in ['gaussian_s_thick_s2_energy']:
                    if self.tfinal == 1:
                        plt.ylim(1e-8, 1e-1)
                    elif self.tfinal == 10:
                        plt.ylim(1e-6, 1e0)
                    elif self.tfinal == 100:
                        plt.ylim(1e-6, 1e0)

            elif self.source_name == 'su_olson_thick_s2_energy':
                if self.tfinal == 100:
                    plt.ylim(5e-6,1e0)


            plt.semilogy(self.Ms, self.energy_RMS, self.line_mkr + self.mkr, c = self.clr, mfc = self.mfc)

            self.find_c_semilog()


        # plt.savefig(self.plot_file_path / "RMS_plots" / f"{self.source_name}_t={self.tfinal}_RMSE_vs_cells.pdf")
        file_path_string = str(self.plot_file_path) + '/' "RMS_plots" '/' + f"{self.source_name}_t={self.tfinal}_RMSE_vs_Ms"
        show_loglog(file_path_string, xlimleft, xlimright)
        print(self.source_name)
        print('uncollided', self.uncollided)
        print('moving', self.moving)
        print("--- --- --- --- --- ")
        

        plt.show(block = False)
        #######################################################################
            
        
    def plot_RMS_vs_times(self, fign = 1, clear = False):
        plt.ion()
        plt.figure(fign)
        if clear == True:
            plt.clf()
        plt.xlabel("average run time (s)", fontsize = 20)
        plt.ylabel("RMSE", fontsize = 20)
        self.times = self.times[1:]
        self.RMS = self.RMS[1:]
        plt.loglog(self.times, self.RMS, self.line_mkr + self.mkr, c = self.clr, mfc = self.mfc)
        file_path_string = str(self.plot_file_path) + '/' + f"{self.source_name}_t={self.tfinal}_times_vs_cells"
        xlimleft = 10e-1
        xlimright = 10e3
        show_loglog_time(file_path_string, xlimleft, xlimright)
        plt.show(block = False)
        
    def plot_RMS_vs_times_Ms(self, fign = 1, clear = False):
        plt.ion()
        plt.figure(fign)
        if clear == True:
            plt.clf()
        plt.xlabel("average run time (s)", fontsize = 20)
        plt.ylabel("RMSE", fontsize = 20)

        plt.loglog(self.times, self.RMS, self.line_mkr + self.mkr, c = self.clr, mfc = self.mfc)
        file_path_string = str(self.plot_file_path) + '/' + f"{self.source_name}_t={self.tfinal}_times_vs_cells"
        xlimleft = 10e-1
        xlimright = 10e3
        show_loglog_time(file_path_string, xlimleft, xlimright)
        plt.show(block = False)
        
    def plot_best(self, mkr, clr):
        plt.ion()
        plt.figure(1)
        
        plt.xlabel("cells", fontsize = 20)
        plt.ylabel("RMSE", fontsize = 20)
        
        if self.source_name == 'plane_IC':
            name = 'plane pulse'
        elif self.source_name == 'square_IC':
            name = 'square pulse'
        elif self.source_name == 'square_s':
            self.cells = self.cells[1:]
            self.RMS = self.RMS[1:]
            self.angles = self.angles[1:]
            
            name = 'square source'
        elif self.source_name =='gaussian_IC':
            name = 'Gaussian pulse'
        elif self.source_name == 'gaussian_s':
            name = 'Gaussian source'
            intercept = self.find_c()
            order_triangle(3, 6, 6, intercept, 2, 1.7)
        elif self.source_name == 'MMS':
            name = 'MMS'

            
        plt.loglog(self.cells, self.RMS, mkr, c = clr, label = name)
        
        file_path_string = str(self.plot_file_path) + '/' + "all_bestt={self.tfinal}_cells_vs_RMSE"
        xlimleft = 1
        xlimright = 33

        
        plt.ylim(10e-8, 10e-3)
        
        if self.source_name == 'gaussian_s':
            plt.legend()
            
        show_loglog(file_path_string, xlimleft, xlimright)
        plt.show(block = False)
     
    def plot_best_static(self, mkr, clr):
         plt.ion()
         plt.figure(1)
         
         plt.xlabel("cells", fontsize = 20)
         plt.ylabel("RMSE", fontsize = 20)
         
         if self.source_name == 'plane_IC':
             name = 'plane pulse'
         elif self.source_name == 'square_IC':
             name = 'square pulse'
         elif self.source_name == 'square_s':
             self.cells = self.cells[1:]
             self.RMS = self.RMS[1:]
             self.angles = self.angles[1:]
             
             name = 'square source'
         elif self.source_name =='gaussian_IC':
             name = 'Gaussian pulse'
         elif self.source_name == 'gaussian_s':
             name = 'Gaussian source'
             intercept = self.find_c()
             order_triangle(3, 6, 6, intercept, 2, 500)
         elif self.source_name == 'MMS':
             name = 'MMS'

             
         plt.loglog(self.cells, self.RMS, mkr, c = clr, label = name)
         
         file_path_string = str(self.plot_file_path) + '/' + "all_best_static_t={self.tfinal}_cells_vs_RMSE"
         xlimleft = 1
         xlimright = 33

         
         plt.ylim(10e-8, 10e-3)
         
         if self.source_name == 'gaussian_s':
             plt.legend(loc='lower right')
             
         show_loglog(file_path_string, xlimleft, xlimright)
         plt.show(block = False)
        
    def plot_coefficients(self, tfinal,  M, source_name,  N_spaces, problem_name, rad_or_transport,
     x0_or_sigma, c, cv0, uncollided, s2, mat_or_rad, moving, line, legend = True, pc = 0, fign = 1, ifshow = False):

        print('here')
        # if tfinal == 100.0 and problem_name == 'transfer_const_cv=0.03':
        #     file_name = 'run_data_crc_nov15.hdf5'
        #     data = load_sol(problem_name, source_name, rad_or_transport, c, s2, cv0)
        # else:
        if (source_name == 'square_s' or x0_or_sigma == 0.375) and problem_name != 'su_olson_thick_s2':
            data = load_sol(problem_name, source_name, rad_or_transport, c, s2, cv0)
        else:
            if source_name == 'gaussian_s' and x0_or_sigma == 0.5:
                if problem_name in ['transfer_const_cv=0.03', 'transfer_const_cv=0.03_s2']:
                    if tfinal >=10.0:
                        file_name = 'run_data_crc_dec7-4.hdf5'
                    else:
                        file_name = 'run_data_crc_nov23.hdf5'
                else:
                    file_name = 'run_data_crc_nov23.hdf5'
                    
            elif problem_name == 'su_olson_thick_s2':
                file_name = 'run_data_crc_dec7-4.hdf5'
            # elif self.tfinal >= 10.0 and self.source_name == 'gaussian_s' and x0_or_sigma == 0.5 and problem_name in ['transfer_const_cv=0.03', 'transfer_const_cv=0.03_s2'] :
                
            print('loading', file_name)
            data = load_sol(problem_name, source_name, rad_or_transport, c, s2, cv0, file_name)

        

        self.M_coeff = M
        self.j_matrix = np.zeros((len(N_spaces), (M+1)))
        self.label_list = N_spaces
        self.line = line
        self.legend = legend
        self.N_angles = []
        self.mat_or_rad = mat_or_rad
        self.fign = fign
        self.tfinal = tfinal
        self.plot_counter = pc
        self.N_spaces = N_spaces
        self.ifshow = ifshow
        print(self.tfinal, ifshow)

        for count, N_space in enumerate(N_spaces):
            print(uncollided,'uncollided',moving,'moving')
            data.call_sol(tfinal, M, x0_or_sigma, N_space, mat_or_rad, uncollided, moving)
            N_ang = np.shape(data.coeff_mat)[0]
            self.N_angles.append(N_ang)

            
            for k in range(N_space):
                coeff_data = coeff_con(data.ws, data.coeff_mat, N_ang, M, k)
                self.j_matrix[count] += np.abs(coeff_data)/N_space
            
        
        for j in range(M+1):
            plt.ion()
            plt.figure(fign)  
            xdata = np.array(N_spaces)
            ydata = np.array(self.j_matrix[:,j])

            # plt.loglog(xdata, np.abs(ydata), "-o", mfc = 'none', label =f"j={j}")
        self.file_path_string = str(self.plot_file_path) + '/'  + "final_coefficient_convergence" + "/" + str(self.source_name) + '/' + problem_name + "_"  +'_' + self.mat_or_rad + "_" + source_name + "_M=" + str(M) + "x0_or_sigma_" + str(x0_or_sigma) + "_S2" * s2
        # show_loglog(self.file_path_string, N_spaces[0]-1, N_spaces[-1] + 2)
        # plt.show()
        
        
        # plt.legend()
        # plt.show()
    
    def plot_coeff_boyd(self):

        plt.figure(self.fign)
        xdata = np.linspace(0,self.M_coeff, self.M_coeff+1)

        label_list = self.label_list
        mkr_list_old = ['o', '^', 's', 'p', '*', 'x', 'd'] 

            
        mkr_list = [x + self.line for x in mkr_list_old]

        
        self.Ms = xdata
        
        if self.tfinal == 0.1:
            self.plot_counter = 0
        elif self.tfinal == 0.31623:
            self.plot_counter = 1
        elif self.tfinal == 1.0:
            self.plot_counter = 2
        elif self.tfinal == 3.16228:
            self.plot_counter = 3
        elif self.tfinal == 10.0:
            self.plot_counter = 4
        elif self.tfinal == 31.6228:
            self.plot_counter = 5
        elif self.tfinal == 100.0:
            self.plot_counter = 6
        
        
        
        for ij in range(len(self.j_matrix[:,0])):
            _ = plt.semilogy(xdata, self.j_matrix[ij], mkr_list[self.plot_counter], label = str(self.N_spaces[0]) + ' cells,' + ' t = ' + str(self.tfinal), mfc = 'none', c = 'b')
            self.RMS = self.j_matrix[ij]        # this is a hack
            self.energy_RMS = self.j_matrix[ij] # this is a hack
            self.find_c_semilog()
            # if self.mat_or_rad == 'rad':
            #     _ = logplot_sn_labels([xdata[-1]], [self.j_matrix[ij][-1]], np.ones(1)*self.N_angles[ij], 0.3, 3)

        plt.xlabel("M", fontsize = 20)
        plt.ylabel(r"avg. $|c_n|$", fontsize = 20)
        
        if self.legend == True:
            plt.legend()

        plt.ylim(1e-14, 1e-0)
        if self.ifshow == True: 
            _ = show_loglog(self.file_path_string + "_boyd", 0.001, self.M_coeff + 4)
        plt.show()


        # plt.figure(self.fign + 2)
        # plt.ylim(1e-14, 1e-0)
        # for ij in range(len(self.j_matrix[:,0])):
        #     _ = plt.loglog(xdata, self.j_matrix[ij], mkr_list[self.plot_counter], label = 't = ' + str(self.tfinal), mfc = 'none', c = 'b')
        #     self.RMS = self.j_matrix[ij]
        #     if self.mat_or_rad == 'rad':
        #         logplot_sn_labels([xdata[-1]], [self.j_matrix[ij][-1]], np.ones(1)*self.N_angles[ij], 0.3, 3)
        # plt.xlabel("M", fontsize = 20)
        # plt.ylabel("RMSE", fontsize = 20)
        # if self.legend == True:
        #     plt.legend()
        # # show_loglog(self.file_path_string + "_boyd_loglog" ,1, self.M_coeff + 4)
        # # print("convergence order", "%.2f" % convergence(self.RMS[2], xdata[2], self.RMS[3], xdata[3]))

        # plt.show()

    
    def plot_bench(self, tfinal, source_name, fign):
        plt.figure(fign)
        fntsize = 30
        plt.ion()
        plt.xlabel("x", fontsize = fntsize)
        plt.ylabel("scalar flux", fontsize = fntsize)
        file_path_string = str(self.plot_file_path) + "/" + "benchmark_plots"
        npnts = 10000
        if source_name == "plane_IC":
            source_type = np.array([1,0,0,0,0,0,0,0])
            x0 = 1e-11
            bench = load_bench(source_type, tfinal, x0, 1.0, False)
            xs = np.linspace(0, tfinal, npnts)
            interp_bench = bench(xs)[0]
            uncol = bench(xs)[1]
            interp_bench2 = np.append(interp_bench, np.array([0.0]))
            uncol2 = np.append(uncol, np.array([0.0]))
            xs2 = np.append(xs, np.array([tfinal + .0000001]))
            plt.plot(xs2, interp_bench2, "-k")
            plt.plot(-xs2, interp_bench2, "-k")
            if max(uncol2)>= 1e-5:
                plt.plot(xs2, uncol2, "-.k")
                plt.plot(-xs2, uncol2, "-.k")
            
            show(file_path_string + f"/plane_IC_t_{tfinal}_benchmark")
            plt.show(block = False)
        elif source_name == "square_IC":
            source_type = np.array([0,1,0,0,0,0,0,0])
            x0 = 0.5
            bench = load_bench(source_type, tfinal, x0, 1.0, False)
            xs = np.linspace(0, tfinal + x0, npnts)
            interp_bench = bench(xs)[0]
            uncol = bench(xs)[1]
            plt.plot(xs, interp_bench, "-k")
            plt.plot(-xs, interp_bench, "-k")
            if max(uncol)>= 1e-5:
                plt.plot(xs, uncol, "-.k")
                plt.plot(-xs, uncol, "-.k")
            show(file_path_string + f"/square_IC_t_{tfinal}_benchmark")
            plt.show()
        elif source_name == "square_source":
            source_type = np.array([0,0,1,0,0,0,0,0])
            x0 = 0.5
            bench = load_bench(source_type, tfinal, x0, 1.0, False)
            xs = np.linspace(0, tfinal + x0, npnts)
            interp_bench = bench(xs)[0]
            uncol = bench(xs)[1]
            plt.plot(xs, interp_bench, "-k")
            plt.plot(-xs, interp_bench, "-k")
            if tfinal == 1 or tfinal == 5:
                plt.plot(xs, uncol, "-.k")
                plt.plot(-xs, uncol, "-.k")
            show(file_path_string + f"/square_source_t_{tfinal}_benchmark")
            plt.show()
        elif source_name == "square_source_s2":
            source_type = np.array([0,0,0,0,0,0,0,0,1,0])
            x0 = 0.5
            bench = load_bench(source_type, tfinal, x0, 1.0, False)
            xs = np.linspace(0, tfinal + x0, npnts)
            interp_bench = bench(xs)[0]
            uncol = bench(xs)[1]
            plt.plot(xs, interp_bench, "-k")
            plt.plot(-xs, interp_bench, "-k")
            # if tfinal == 1 or tfinal == 5:
            #     plt.plot(xs, uncol, "-.k")
            #     plt.plot(-xs, uncol, "-.k")
            source_type = np.array([0,0,0,0,0,0,0,0,0,1])
            x0 = 0.5
            bench = load_bench(source_type, tfinal, x0, 1.0, False)
            xs = np.linspace(0, tfinal + x0, npnts)
            interp_bench = bench(xs)[0]
            uncol = bench(xs)[1]
            plt.plot(xs, interp_bench, ":k")
            plt.plot(-xs, interp_bench, ":k")
            show(file_path_string + f"/s2_square_source_t_{tfinal}_benchmark")
            plt.show()
        elif source_name == "square_source_s2_thick":
            source_type = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0])
            x0 = 400
            bench = load_bench(source_type, tfinal, x0, 1.0, False)
            xs = np.linspace(0, 500, npnts)
            interp_bench = bench(xs)[0]
            uncol = bench(xs)[1]
            plt.plot(xs, interp_bench, "-k")
            plt.plot(-xs, interp_bench, "-k")
            # if tfinal == 1 or tfinal == 5:
            #     plt.plot(xs, uncol, "-.k")
            #     plt.plot(-xs, uncol, "-.k")
            source_type = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1])
            x0 = 400
            bench = load_bench(source_type, tfinal, x0, 1.0, False)
            xs = np.linspace(0, 500, npnts)
            interp_bench = bench(xs)[0]
            uncol = bench(xs)[1]
            plt.plot(xs, interp_bench, ":k")
            plt.plot(-xs, interp_bench, ":k")
            show(file_path_string + f"/s2_square_source_thick_t_{tfinal}_benchmark_view1")
            plt.show()

            plt.figure(fign+100)
            source_type = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0])
            x0 = 400
            bench = load_bench(source_type, tfinal, x0, 1.0, False)
            interp_bench = bench(xs)[0]
            plt.plot(xs, interp_bench, "-k")
            plt.plot(-xs, interp_bench, "-k")
            source_type = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1])
            x0 = 400
            bench = load_bench(source_type, tfinal, x0, 1.0, False)
            interp_bench = bench(xs)[0]
            plt.plot(xs, interp_bench, ":k")
            plt.plot(-xs, interp_bench, ":k")
            plt.xlim(380, 420)
            show(file_path_string + f"/s2_square_source_thick_t_{tfinal}_benchmark_view2")
            plt.show()
        elif source_name == "gaussian_source_s2":
            source_type = np.array([0,0,0,0,0,0,0,0,0,0,1,0])
            x0 = 0.5
            bench = load_bench(source_type, tfinal, x0, 1.0, False)
            xs = np.linspace(0, tfinal + x0, npnts)
            interp_bench = bench(xs)[0]
            uncol = bench(xs)[1]
            plt.plot(xs, interp_bench, "-k")
            plt.plot(-xs, interp_bench, "-k")
            # if tfinal == 1 or tfinal == 5:
            #     plt.plot(xs, uncol, "-.k")
            #     plt.plot(-xs, uncol, "-.k")
            source_type = np.array([0,0,0,0,0,0,0,0,0,0,0,1])
            x0 = 0.5
            bench = load_bench(source_type, tfinal, x0, 1.0, False)
            xs = np.linspace(0, tfinal + x0, npnts)
            interp_bench = bench(xs)[0]
            uncol = bench(xs)[1]
            plt.plot(xs, interp_bench, ":k")
            plt.plot(-xs, interp_bench, ":k")
            show(file_path_string + f"/s2_gaussian_source_t_{tfinal}_benchmark")
            plt.show()
        elif source_name == "gaussian_source_s2_thick":
            source_type = np.array([0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0])
            
            x0 = 1500
            bench = load_bench(source_type, tfinal, x0, 1.0, False)
            xs = np.linspace(0, tfinal + x0, npnts)
            interp_bench = bench(xs)[0]
            uncol = bench(xs)[1]
            plt.plot(xs, interp_bench, "-k")
            plt.plot(-xs, interp_bench, "-k")
            # if tfinal == 1 or tfinal == 5:
            #     plt.plot(xs, uncol, "-.k")
            #     plt.plot(-xs, uncol, "-.k")
            source_type = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0])
            x0 = 1500
            bench = load_bench(source_type, tfinal, x0, 1.0, False)
            xs = np.linspace(0, tfinal + x0, npnts)
            interp_bench = bench(xs)[0]
            uncol = bench(xs)[1]
            plt.plot(xs, interp_bench, ":k")
            plt.plot(-xs, interp_bench, ":k")

            show(file_path_string + f"/s2_gaussian_source_thick_t_{tfinal}_benchmark")
            plt.show()


        elif source_name == "gaussian_IC":
            source_type = np.array([0,0,0,1,0,0,0,0])
            x0 = 4
            bench = load_bench(source_type, tfinal, x0, 1.0, False)
            if tfinal == 1:
                xs = np.linspace(0, 3.0, npnts)
            elif tfinal == 5:
                xs= np.linspace(0, 5.0, npnts)
            elif tfinal == 10:
                xs = np.linspace(0, 10.0, npnts)
            interp_bench = bench(xs)[0]
            uncol = bench(xs)[1]
            if tfinal == 1:
                plt.plot(xs, uncol, "-.k")
                plt.plot(-xs, uncol, "-.k")
            plt.plot(xs, interp_bench, "-k")
            plt.plot(-xs, interp_bench, "-k")
            show(file_path_string + f"/gaussian_IC_t_{tfinal}_benchmark")
        elif source_name == "gaussian_source":
            source_type = np.array([0,0,0,0,0,1,0,0])
            x0 = 4
            bench = load_bench(source_type, tfinal, x0, 1.0, False)
            if tfinal == 1:
                xs = np.linspace(0,3.0,npnts)
            elif tfinal == 5:
                xs = np.linspace(0, 5.0,npnts)
            elif tfinal == 10:
                xs = np.linspace(0, 10.0,npnts)
            interp_bench = bench(xs)[0]
            uncol = bench(xs)[1]
            if tfinal == 1 or tfinal == 5:
                plt.plot(xs, uncol, "-.k")
                plt.plot(-xs, uncol, "-.k")
            plt.plot(xs, interp_bench, "-k")
            plt.plot(-xs, interp_bench, "-k")
            show(file_path_string + f"/gaussian_source_t_{tfinal}_benchmark")
        elif source_name == "gaussian_IC_2D":
            plt.xlabel("r", fontsize = fntsize)
            source_type = np.array([0,0,0,0,0,0,1,0])
            x0 = 0.5
            xs = np.linspace(0, tfinal + 1/x0, npnts)
            bench = load_bench(source_type, tfinal, x0, 1.0, False)
            interp_bench = bench(xs)[0]
            uncol = bench(xs)[1]
            if tfinal == 1:
                plt.plot(xs, uncol, "--k")
            plt.plot(xs, interp_bench, "-k")
            show(file_path_string + f"/gaussian_IC_2D_t_{tfinal}_benchmark")
        elif source_name == "line_source":
            plt.xlabel("r", fontsize = fntsize)
            source_type = np.array([0,0,0,0,0,0,0,1])
            x0 = 0.5
            xs = np.linspace(0, tfinal, npnts)
            bench = load_bench(source_type, tfinal, x0, 1.0, False)
            interp_bench = bench(xs)[0]
            uncol = bench(xs)[1]
            if tfinal == 1 or tfinal == 5:
                plt.plot(xs, uncol, "--k")
            plt.plot(xs, interp_bench, "-k")
            show(file_path_string + f"/line_source_t_{tfinal}_benchmark")
        # elif source_name == "MMS":
        #     plt.xlabel("x", fontsize = fntsize)
        #     source_type = np.array([0,0,0,0,1,0,0,0])
        #     x0 = 0.1
        #     xs = np.linspace(0, tfinal + x0, npnts)
        #     bench = load_bench(source_type, tfinal, x0, 1.0, False)
        #     interp_bench = bench(xs)[0]
        #     plt.plot(xs, interp_bench, "-k")
        #     plt.plot(-xs, interp_bench, "-k")
        #     show(file_path_string + f"/MMS_t_{tfinal}_benchmark")
        elif source_name == "square_IC_c_not_one":
            plt.figure(fign)
            plt.xlabel("x", fontsize = fntsize)
            source_type = np.array([0,1,0,0,0,0,0,0])
            x0 = 0.5
            xs = np.linspace(0, 1.25 + 0.625, npnts)
            bench = load_bench(source_type, 1.25, x0, 0.8, True)
            interp_bench = 0.8*math.exp(-(1-0.8)*1.25)*bench(xs*0.8)[0]
            uncol = 0.8*math.exp(-(1-0.8)*1.25)*bench(0.8*xs)[1]
            plt.plot(xs, interp_bench, "-k")
            plt.plot(-xs, interp_bench, "-k")
            plt.plot(xs, uncol, "--k")
            plt.plot(-xs, uncol, "--k")
            show(file_path_string + "/sq_IC_c=0.8_t_1.25_benchmark")
            plt.show()

            plt.figure(fign + 10)
            plt.xlabel("x", fontsize = fntsize)
            source_type = np.array([0,1,0,0,0,0,0,0])
            x0 = 0.5
            tactual = 1/1.2
            xs = np.linspace(0, tactual + 0.5/1.2, npnts)
            bench = load_bench(source_type, tactual, x0, 1.2, True)
            interp_bench = 1.2*math.exp(-(1-1.2)*tactual)*bench(xs*1.2)[0]
            uncol = 1.2*math.exp(-(1-1.2)*tactual)*bench(1.2*xs)[1]
            plt.plot(xs, interp_bench, "-k")
            plt.plot(-xs, interp_bench, "-k")
            plt.plot(xs, uncol, "--k")
            plt.plot(-xs, uncol, "--k")
            show(file_path_string + f"/sq_IC_c=1.2_t_{tactual}_benchmark")
            plt.show()

        elif source_name == "gaussian_IC_c_not_one":
            plt.figure(fign)
            plt.xlabel("x", fontsize = fntsize)
            source_type = np.array([0,0,0,1,0,0,0,0])
            x0 = 4
            xs = np.linspace(0, 1.25 + x0, npnts)
            bench = load_bench(source_type, 1.25, x0, 0.8, True)
            interp_bench = 0.8*math.exp(-(1-0.8)*1.25)*bench(xs*0.8)[0]
            uncol = 0.8*math.exp(-(1-0.8)*1.25)*bench(0.8*xs)[1]
            plt.plot(xs, interp_bench, "-k")
            plt.plot(-xs, interp_bench, "-k")
            plt.plot(xs, uncol, "--k")
            plt.plot(-xs, uncol, "--k")
            show(file_path_string + "/gauss_IC_c=0.8_t_1.25_benchmark")
            plt.show()

            plt.figure(fign + 10)
            plt.xlabel("x", fontsize = fntsize)
            source_type = np.array([0,0,0,1,0,0,0,0])
            x0 = 4
            tactual = 1/1.2
            xs = np.linspace(0, tactual + x0, npnts)
            bench = load_bench(source_type, tactual, x0, 1.2, True)
            interp_bench = 1.2*math.exp(-(1-1.2)*tactual)*bench(xs*1.2)[0]
            uncol = 1.2*math.exp(-(1-1.2)*tactual)*bench(1.2*xs)[1]
            plt.plot(xs, interp_bench, "-k")
            plt.plot(-xs, interp_bench, "-k")
            plt.plot(xs, uncol, "--k")
            plt.plot(-xs, uncol, "--k")
            show(file_path_string + f"/gauss_IC_c=1.2_t_{tactual}_benchmark")
            plt.show()





        

        
    
        
        
        
        
        
        
            

    
    
    
    
    
    
    


