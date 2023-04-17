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
        if self.x0 == 2.5:
            self.olson_henderson_bench(self.tfinal)
            plt.figure(9)
            plt.plot(self.xs, self.phi, '-.', label = 'scalar flux', mfc = 'none')
            plt.legend()
            plt.show()
     
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
                # plt.plot(self.xs, self.psi[-1,:], '-^')
                plt.plot(self.xs, fsol(self.xs+self.x0, -1), 'rx')
                plt.show()
            elif solver.sigma_func[1] == 1:
                self.steady_state_gaussian_benchmark()
            
            elif solver.sigma_func[2] == 1 or solver.sigma_func[3] == 1:
                self.siewert_bench(solver.sigma_func)

                
      


    def get_results(self, solver):
        self.xs = solver.xs
        self.phi = solver.phi
        self.e = solver.e
        self.psi = solver.psi
        self.ws = solver.ws
        self.mus = solver.angles
        self.x0 = solver.x0
        self.exit_dist = solver.exit_dist
        self.tfinal = solver.tfinal

        
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
        plt.plot(self.xs, phi_sol, 'kx', mfc = 'none', label = 'scalar flux benchmark')
        plt.legend()

    def siewert_bench(self, sigma):
        self.psibenchpsis = np.array([0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        self.psi0bench = np.array([0.58966, 0.53112, 0.44328, 0.38031, 0.33296, 
                                0.29609, 0.26656, 0.24239, 0.22223, 0.20517, 0.19055])
        self.psi1bench = np.array([0.6075e-5, 0.62952e-5, 0.96423e-5, 0.16234e-4, 
        0.43858e-4, 0.16937e-3, 0.57347e-3, 0.15128e-2, 0.32437e-2, 0.59604e-2, 0.97712e-2 ])

        self.psi0benchinf = np.array([0.89780, 0.88784, 0.86958, 0.85230, 0.83550, 0.81900, 0.80278, 
                                    0.78649, 0.77043, 0.75450, 0.73872])

        self.psi1benchinf = np.array([0.10220, 0.11216, 0.13042, 0.14770, 0.16450, 0.18100, 0.19732, 
                                    0.21351, 0.22957, 0.24550, 0.26128])



        if sigma[2] == 1: 
            resultsfigns = [9,10]
            plt.figure(9)
            plt.plot(self.psibenchpsis, self.psi0bench, 'kx', label = 'benchmark s = 1')
            plt.legend()
            plt.show()
            plt.figure(10)
            plt.plot(self.psibenchpsis, self.psi1bench, 'kx', label = 'benchmark s = 1')
            plt.legend()
            plt.show()
        elif sigma[3] == 1: 
            resultsfigns = [11,12]
            plt.figure(11)
            plt.plot(self.psibenchpsis, self.psi0benchinf, 'kx', label = 'benchmark s = inf')
            plt.legend()
            plt.show()
            plt.figure(12)
            plt.plot(self.psibenchpsis, self.psi1benchinf, 'kx', label = 'benchmark s = inf')
            plt.legend()
            plt.show()





        plt.figure(resultsfigns[0])
        plt.plot(-self.mus, self.exit_dist[:,0], '--b', mfc = 'none', label = 'left exit distribution')
        # plt.plot(self.mus, self.exit_dist[:,-1], '-or', mfc = 'none', label = 'right exit distribution')
        plt.xlabel(r'$\mu$')
        plt.ylabel(r'$\phi$')
        plt.xlim(0.0, 1.1)
        plt.legend()
        plt.show()


        plt.figure(resultsfigns[1])
        # plt.plot(self.mus, self.exit_dist[:,0], '-ob', mfc = 'none', label = 'left exit distribution')
        plt.plot(self.mus, self.exit_dist[:,-1], '--r', mfc = 'none', label = 'right exit distribution')
        plt.xlabel(r'$\mu$')
        plt.ylabel(r'$\phi$')
        plt.legend()
        plt.xlim(0.1, 1.1)
        plt.show()


    def olson_henderson_bench(self, tfinal):
        self.xs_bench = np.array([0, 0.5, 1.0, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5,
                                    9, 9.5, 10 ]) - 5.0 
        self.phi_bench = self.xs_bench*0

        if tfinal == 1.0:
            self.phi_bench = np.array([0.0, 0.0, 0, 0, 0.052, 0.476, 0.899, 0.952, 0.952, 0.952, 0.952, 0.952, 0.952, 0.952, 0.899, 0.476, 0.052, 0,
                                    0, 0, 0]) 
        elif tfinal == 5.0:
            self.phi_bench=np.array([0.051, 0.138, 0.290, 0.562, 1.035, 1.968, 2.900, 3.371, 3.636, 3.771, 
            3.812, 3.771, 3.636, 3.371, 2.900, 1.968, 1.035, 0.562, 0.290, 0.138, 0.051])
        plt.figure(9)
        plt.plot(self.xs_bench, self.phi_bench, 'kx', label = 'benchmark scalar flux')
        plt.legend()
        plt.xlabel('x')
        plt.ylabel(r'$\phi$')
        plt.show()