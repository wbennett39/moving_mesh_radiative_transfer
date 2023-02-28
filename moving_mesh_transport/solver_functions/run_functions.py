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
from scipy.special import erf
import math
import numpy as np

# I could take all of the plotting out of the solver, group problems here by 
# infinite medium, finite etc. to simplify things

class run:
    def __init__(self):
        self.data_folder = Path("moving_mesh_transport/input_scripts")
    
    def load(self, problem_type = 'transport'):
        config_file_path = self.data_folder / f"{problem_type}.yaml"
        mesh_config_file_path = self.data_folder / "mesh_parameters.yaml"
        with open(config_file_path, 'r') as file:
            self.parameters = yaml.safe_load(file)
            file.close()

        with open(mesh_config_file_path, 'r') as file:
            self.mesh_parameters = yaml.safe_load(file)
            file.close()
   
    def h(self):
        print(" choose problem type : 'transport','rad_transfer','su_olson','s2_rad_transfer','s2_rad_transfer_thick','rad_transfer_thick','config'")
        
    def plane_IC(self, uncollided = True, moving = True, All = False):
        plt.ion()
        # plt.figure(1)
        source_name = "plane_IC"
        print("---  ---  ---  ---  ---  ---  ---")
        print("running plane IC")
        print("---  ---  ---  ---  ---  ---  ---")
        
        solver = main_class(source_name, self.parameters, self.mesh_parameters) 
        if All == True:
            solver.main(True, True)
            solver.main(False, True)
            solver.main(True, False)
            solver.main(False, False)
        else:
            solver.main(uncollided, moving)
            self.get_results(solver)
        # plt.title("plane IC")
        # plt.legend()
        # plt.show(block = False)
        
    def square_IC(self, uncollided = True, moving = True, All = False):
        plt.ion()
        # plt.figure(2)
        source_name = "square_IC"
        print("---  ---  ---  ---  ---  ---  ---")
        print("running square IC")
        print("---  ---  ---  ---  ---  ---  ---")
        solver = main_class(source_name, self.parameters, self.mesh_parameters) 
        if All == True:
            solver.main(True, True)
            solver.main(False, True)
            solver.main(True, False)
            solver.main(False, False)
        else:
            solver.main(uncollided, moving)
            self.get_results(solver)
        # plt.title("square IC")
        # plt.legend()
        # plt.show(block = False)
        
    def square_source(self, uncollided = True, moving = True, All = False):
        plt.ion()
        # plt.figure(3)
        source_name = "square_source"
        print("---  ---  ---  ---  ---  ---  ---")
        print("running square source")
        print("---  ---  ---  ---  ---  ---  ---")
        solver = main_class(source_name, self.parameters, self.mesh_parameters) 
        if All == True:
            solver.main(True, True)
            solver.main(False, True)
            solver.main(True, False)
            solver.main(False, False)
        else:
            solver.main(uncollided, moving)
            self.get_results(solver)
        # plt.title("square source")
        # plt.legend()
        # plt.show(block = False)
        
    def gaussian_IC(self, uncollided = True, moving = True, All = False):
        plt.ion()
        # plt.figure(4)
        source_name = "gaussian_IC"
        print("---  ---  ---  ---  ---  ---  ---")
        print("running Gaussian IC")
        print("---  ---  ---  ---  ---  ---  ---")
        solver = main_class(source_name, self.parameters, self.mesh_parameters) 
        if All == True:
            solver.main(True, True)
            solver.main(False, True)
            solver.main(True, False)
            solver.main(False, False)
        else:
            solver.main(uncollided, moving)
            self.get_results(solver)
        # plt.title("Gaussian IC")
        # plt.legend()
        # plt.show(block = False)
        
    def gaussian_source(self, uncollided = True, moving = True, All = False):
        plt.ion()
        plt.figure(5)
        source_name = "gaussian_source"
        print("---  ---  ---  ---  ---  ---  ---")
        print("running Gaussian source")
        print("---  ---  ---  ---  ---  ---  ---")
        solver = main_class(source_name, self.parameters, self.mesh_parameters) 
        if All == True:
            solver.main(True, True)
            solver.main(False, True)
            solver.main(True, False)
            solver.main(False, False)
        else:
            solver.main(uncollided, moving)
            self.get_results(solver)
        # plt.title("Gaussian source")
        # plt.legend()
        # plt.show(block = False)
        
    def MMS(self, uncollided = False, moving = True, All = False):
        plt.ion()
        # plt.figure(6)
        source_name = "MMS"
        print("---  ---  ---  ---  ---  ---  ---")
        print("running MMS problem")
        print("---  ---  ---  ---  ---  ---  ---")
        solver = main_class(source_name, self.parameters, self.mesh_parameters) 
        if All == True:
            solver.main(True, True)
            solver.main(False, True)
            solver.main(True, False)
            solver.main(False, False)
        else:
            solver.main(uncollided, moving)
            self.get_results(solver)
        # plt.title("MMS")
        # plt.legend()
        # plt.show(block = False)


    def boundary_source(self, uncollided = False, moving = True, All = False):
        plt.ion()
        # plt.figure(6)
        source_name = "boundary_source"
        print("---  ---  ---  ---  ---  ---  ---")
        print("running boundary source problem")
        print("---  ---  ---  ---  ---  ---  ---")
        solver = main_class(source_name, self.parameters, self.mesh_parameters) 
        if All == True:
            solver.main(True, True)
            solver.main(False, True)
            solver.main(True, False)
            solver.main(False, False)
        else:
            solver.main(uncollided, moving)
            self.get_results(solver)
            import numpy as np
            f = lambda x: np.exp(-x**2 /(2 * 0.5**2))
            fsol = lambda x, mu: np.exp(x * 1/mu) 
            plt.figure(3)
            if solver.sigma_func[0] == 1:
                plt.plot(self.xs, self.psi[-1,:], '-^')
                plt.plot(self.xs, fsol(self.xs+self.x0, -1), 'rx')
                plt.show()
            elif solver.sigma_func[1] == 1:
                self.steady_state_gaussian_benchmark()

                
      


    def get_results(self, solver):
        self.xs = solver.xs
        self.phi = solver.phi
        self.e = solver.e
        self.psi = solver.psi
        self.ws = solver.ws
        self.mus = solver.angles
        self.x0 = solver.x0

        
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


    
    def steady_state_gaussian_benchmark(self):
        f = lambda x1, x2, sigma: np.sqrt(np.pi/2)*sigma*(erf(x2/math.sqrt(2)/sigma) - erf(x1/math.sqrt(2)/sigma))
        fsol = lambda x, x0, mu: np.exp((f(-x0, x, 0.5)) / mu) 
        phi_sol = self.xs * 0 
        for ix in range(self.xs.size):
            for l in range(self.ws.size):
                if self.mus[l] < 0:
                    phi_sol[ix] += self.ws[l] * fsol(self.xs[ix], self.x0, self.mus[l])
            


        # plt.plot(self.xs, self.psi[-1,:], '-^')
        plt.plot(self.xs, phi_sol, 'ko', mfc = 'none', label = 'scalar flux')
        plt.legend()
