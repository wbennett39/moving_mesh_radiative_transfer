#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 24 13:42:55 2022

@author: bennett
"""
import matplotlib.pyplot as plt
from ..solver import main_class
from pathlib import Path
import yaml


class run:
    def __init__(self):
        self.data_folder = Path("moving_mesh_transport/input_scripts")
    
    def load(self, problem_type = 'transport'):
        config_file_path = self.data_folder / f"{problem_type}.yaml"
        with open(config_file_path, 'r') as file:
            self.parameters = yaml.safe_load(file)
            file.close()
   
    def h(self):
        print(" choose problem type : 'transport','rad_transfer','su_olson','s2_rad_transfer','s2_rad_transfer_thick','rad_transfer_thick','config'")
        
    def plane_IC(self, uncollided = True, moving = True, All = False):
        plt.ion()
        plt.figure(1)
        source_name = "plane_IC"
        print("---  ---  ---  ---  ---  ---  ---")
        print("running plane IC")
        print("---  ---  ---  ---  ---  ---  ---")
        
        solver = main_class(source_name, self.parameters) 
        if All == True:
            solver.main(True, True)
            solver.main(False, True)
            solver.main(True, False)
            solver.main(False, False)
        else:
            solver.main(uncollided, moving)
            self.get_results(solver)
        plt.title("plane IC")
        plt.legend()
        plt.show(block = False)
        
    def square_IC(self, uncollided = True, moving = True, All = False):
        plt.ion()
        plt.figure(2)
        source_name = "square_IC"
        print("---  ---  ---  ---  ---  ---  ---")
        print("running square IC")
        print("---  ---  ---  ---  ---  ---  ---")
        solver = main_class(source_name, self.parameters) 
        if All == True:
            solver.main(True, True)
            solver.main(False, True)
            solver.main(True, False)
            solver.main(False, False)
        else:
            solver.main(uncollided, moving)
            self.get_results(solver)
        plt.title("square IC")
        plt.legend()
        plt.show(block = False)
        
    def square_source(self, uncollided = True, moving = True, All = False):
        plt.ion()
        plt.figure(3)
        source_name = "square_source"
        print("---  ---  ---  ---  ---  ---  ---")
        print("running square source")
        print("---  ---  ---  ---  ---  ---  ---")
        solver = main_class(source_name, self.parameters) 
        if All == True:
            solver.main(True, True)
            solver.main(False, True)
            solver.main(True, False)
            solver.main(False, False)
        else:
            solver.main(uncollided, moving)
            self.get_results(solver)
        plt.title("square source")
        plt.legend()
        plt.show(block = False)
        
    def gaussian_IC(self, uncollided = True, moving = True, All = False):
        plt.ion()
        plt.figure(4)
        source_name = "gaussian_IC"
        print("---  ---  ---  ---  ---  ---  ---")
        print("running Gaussian IC")
        print("---  ---  ---  ---  ---  ---  ---")
        solver = main_class(source_name, self.parameters) 
        if All == True:
            solver.main(True, True)
            solver.main(False, True)
            solver.main(True, False)
            solver.main(False, False)
        else:
            solver.main(uncollided, moving)
            self.get_results(solver)
        plt.title("Gaussian IC")
        plt.legend()
        plt.show(block = False)
        
    def gaussian_source(self, uncollided = True, moving = True, All = False):
        plt.ion()
        plt.figure(5)
        source_name = "gaussian_source"
        print("---  ---  ---  ---  ---  ---  ---")
        print("running Gaussian source")
        print("---  ---  ---  ---  ---  ---  ---")
        solver = main_class(source_name, self.parameters) 
        if All == True:
            solver.main(True, True)
            solver.main(False, True)
            solver.main(True, False)
            solver.main(False, False)
        else:
            solver.main(uncollided, moving)
            self.get_results(solver)
        plt.title("Gaussian source")
        plt.legend()
        plt.show(block = False)
        
    def MMS(self, uncollided = False, moving = True, All = False):
        plt.ion()
        plt.figure(6)
        source_name = "MMS"
        print("---  ---  ---  ---  ---  ---  ---")
        print("running MMS problem")
        print("---  ---  ---  ---  ---  ---  ---")
        solver = main_class(source_name, self.parameters) 
        if All == True:
            solver.main(True, True)
            solver.main(False, True)
            solver.main(True, False)
            solver.main(False, False)
        else:
            solver.main(uncollided, moving)
            self.get_results(solver)
        plt.title("MMS")
        plt.legend()
        plt.show(block = False)

    def get_results(self, solver):
        self.xs = solver.xs
        self.phi = solver.phi
        self.e = solver.e

        
    def run_all(self):
        # self.plane_IC(True, True)
        # self.plane_IC(True, False)
        # # # self.plane_IC(False, True)        # this doesn't converge
        # self.plane_IC(False, False)
        
        self.square_IC(True, True)
        self.square_IC(True, False)
        self.square_IC(False, True)
        self.square_IC(False, False)
        
        # self.square_source(True, True)
        # self.square_source(True, False)
        # self.square_source(False, True)
        # self.square_source(False, False)
        
        # self.gaussian_IC(True, True)
        # self.gaussian_IC(True, False)
        # self.gaussian_IC(False, True)
        # self.gaussian_IC(False, False)
        
        # self.gaussian_source(True, True)
        # self.gaussian_source(True, False)
        # self.gaussian_source(False, True)
        # self.gaussian_source(False, False)
        
        # self.MMS(False, True)            # only one case is possible for the MMS


    
###############################################################################

###############################################################################


    
