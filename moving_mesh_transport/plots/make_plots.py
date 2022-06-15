#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 15:40:39 2022

@author: bennett
"""

import matplotlib.pyplot as plt
import h5py 
from pathlib import Path
from ..load_bench import load_bench
import numpy as np
from .plot_functions.show import show
from .plot_functions.show_loglog import show_loglog
from .plot_functions.sn_labels import logplot_sn_labels,logplot_sn_labels_2
from .plot_functions.show_loglog_timeplots import show_loglog_time
from scipy.stats import linregress 
from .plot_functions.order_triangle import order_triangle

class rms_plotter:
    
    def __init__(self, tfinal, M, source_name, major):
        data_folder = Path("moving_mesh_transport")
        self.data_file_path = data_folder / 'run_data_RMS.h5'
        
        self.plot_file_path = data_folder / "plots" 
        # self.case_list = ["uncol_mov", "no_uncol_stat", "uncol_stat", "no_uncol_stat"]
        self.tfinal = tfinal
        self.source_type_list  = ["plane_IC", "square_IC", "square_s", "gaussian_IC", "MMS", "gaussian_s", "su_olson"]
        self.source_name = source_name
        self.major = major 
        self.M = M
        self.su_olson_data_file_path = data_folder / "run_data_radiative_transfer_RMS.h5"
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
        slope, intercept, r, p, se = linregress(x, y)
        print(slope, 'slope')
        print(np.exp(slope), 'exp slope')
        return np.exp(intercept)
    
    def find_c_semilog(self):
        """
        Finds the intercept from log log data. Skips the first data point
        """
        x = self.Ms[1:]
        y = np.log(self.RMS[1:])
        slope, intercept, r, p, se = linregress(x, y)
        print(slope, 'slope - c1')
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
        elif self.source_name =='gaussian_s2' or self.source_name == "gaussian_energy_s2":
            self.dest_str = str('gaussian_s' + "/" + "t="  + str(self.tfinal) + "/" + "RMS" + "/" + "S2")
            
        if self.major == 'cells':
            data_str = self.uncollided * ("uncollided_")  + (not(self.uncollided))  * ("no_uncollided_")  + self.moving * ("moving_") + (not(self.moving)) * ("static_") + "M_" + str(self.M)
        elif self.major == 'Ms':
            data_str = 'Ms_' + self.uncollided * ("uncollided_")  + (not(self.uncollided))  * ("no_uncollided_")  + self.moving * ("moving") + (not(self.moving)) * ("static") 
        
        if self.source_name != "su_olson" and self.source_name != "su_olson_energy" and self.source_name != "su_olson_energy_s2" and self.source_name != "su_olson_s2" and self.source_name != "gaussian_s2" and self.source_name != "gaussian_energy_s2":
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
        
        if self.source_name == "su_olson" or self.source_name == "su_olson_energy" or self.source_name == "su_olson_energy_s2" or self.source_name == "su_olson_s2" or self.source_name == "gaussian_s2" or self.source_name == "gaussian_energy_s2":
            print("loading s2 RMS data")
            f_rad = h5py.File(self.su_olson_data_file_path, 'r')
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
    
    def plot_RMS_vs_cells(self, fign = 1, clear = False):
        plt.ion()
        plt.figure(fign)
        if clear == True:
            plt.clf()
    
        plt.xlabel("cells", fontsize = 20)
        plt.ylabel("RMSE", fontsize = 20)
        # plt.title(f"{self.source_name} t = {self.tfinal}")
        #######################################################################
        if self.tfinal == 1:
            xlimleft = 1.5
            xlimright = 40
            
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
                        
            elif self.source_name == 'su_olson' or self.source_name == "su_olson_energy" or self.source_name == "gaussian_s2" or self.source_name == "gaussian_energy_s2":
                xlimright = 65

            if self.source_name !="su_olson_energy" and self.source_name !="gaussian_energy_s2":
                plt.loglog(self.cells, self.RMS, self.line_mkr + self.mkr, c = self.clr, mfc = self.mfc)
            
            if self.source_name == "su_olson_energy" or self.source_name == "gaussian_energy_s2":
                plt.loglog(self.cells, self.energy_RMS, self.line_mkr + self.mkr, c = self.clr, mfc = self.mfc)
            
            if  (self.source_name != 'gaussian_IC' and self.source_name != 'gaussian_s') and (self.uncollided == False and self.moving == False or self.source_name == "MMS" and self.M == 4):
                logplot_sn_labels(self.cells, self.RMS, self.angles, 0.3, fign )
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
            
    
            plt.show(block = False)
        #######################################################################
        elif self.tfinal == 10:
            xlimleft = 1
            xlimright = 35
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
        xlimleft = 1.5
        xlimright = 18
        
        plt.ion()
        plt.figure(fign)
        if clear == True:
            plt.clf()
    
        plt.xlabel("M", fontsize = 20)
        plt.ylabel("RMSE", fontsize = 20)
        # plt.title(f"{self.source_name} t = {self.tfinal}")
        #######################################################################
        
        if self.source_name == 'MMS':
            self.Ms = self.Ms[:]
            self.RMS = self.RMS[:]
            self.angles = self.angles[:]
            xlimright = 10
            plt.ylim = 10e-14
            self.find_c_semilog()
            
        elif self.source_name == "gaussian_IC" or self.source_name == "gaussian_s":
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
        plt.semilogy(self.Ms, self.RMS, self.line_mkr + self.mkr, c = self.clr, mfc = self.mfc)
        
        
        
       
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
            bench = load_bench(source_type, tfinal, x0)
            xs = np.linspace(0, tfinal, npnts)
            interp_bench = bench(xs)[0]
            uncol = bench(xs)[1]
            interp_bench2 = np.append(interp_bench, np.array([0.0]))
            uncol2 = np.append(uncol, np.array([0.0]))
            xs2 = np.append(xs, np.array([tfinal + .0000001]))
            plt.plot(xs2, interp_bench2, "-k")
            plt.plot(-xs2, interp_bench2, "-k")
            if tfinal == 1:
                plt.plot(xs2, uncol2, "-.k")
                plt.plot(-xs2, uncol2, "-.k")
            
            show(file_path_string + f"/plane_IC_t_{tfinal}_benchmark")
            plt.show(block = False)
        elif source_name == "square_IC":
            source_type = np.array([0,1,0,0,0,0,0,0])
            x0 = 0.5
            bench = load_bench(source_type, tfinal, x0)
            xs = np.linspace(0, tfinal + x0, npnts)
            interp_bench = bench(xs)[0]
            uncol = bench(xs)[1]
            plt.plot(xs, interp_bench, "-k")
            plt.plot(-xs, interp_bench, "-k")
            if tfinal == 1:
                plt.plot(xs, uncol, "-.k")
                plt.plot(-xs, uncol, "-.k")
            show(file_path_string + f"/square_IC_t_{tfinal}_benchmark")
            plt.show()
        elif source_name == "square_source":
            source_type = np.array([0,0,1,0,0,0,0,0])
            x0 = 0.5
            bench = load_bench(source_type, tfinal, x0)
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
        elif source_name == "gaussian_IC":
            source_type = np.array([0,0,0,1,0,0,0,0])
            x0 = 4
            bench = load_bench(source_type, tfinal, x0)
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
            bench = load_bench(source_type, tfinal, x0)
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
            bench = load_bench(source_type, tfinal, x0)
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
            bench = load_bench(source_type, tfinal, x0)
            interp_bench = bench(xs)[0]
            uncol = bench(xs)[1]
            if tfinal == 1 or tfinal == 5:
                plt.plot(xs, uncol, "--k")
            plt.plot(xs, interp_bench, "-k")
            show(file_path_string + f"/line_source_t_{tfinal}_benchmark")
        elif source_name == "MMS":
            plt.xlabel("x", fontsize = fntsize)
            source_type = np.array([0,0,0,0,1,0,0,0])
            x0 = 0.1
            xs = np.linspace(0, tfinal + x0, npnts)
            bench = load_bench(source_type, tfinal, x0)
            interp_bench = bench(xs)[0]
            plt.plot(xs, interp_bench, "-k")
            plt.plot(-xs, interp_bench, "-k")
            show(file_path_string + f"/MMS_t_{tfinal}_benchmark")



        

        
    
        
        
        
        
        
        
            

    
    
    
    
    
    
    


