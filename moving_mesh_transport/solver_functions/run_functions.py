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
from scipy import integrate
import scipy.interpolate as interp

# I could take all of the plotting out of the solver, group problems here by 
# infinite medium, finite etc. to simplify things
#-----------------------------------------------------------------------------
# Next steps for this code:

# Special mesh function for Fake Sedov  [x] doesn't seem to work
# Get Sedov loaded  []
# Update readme for info on every problem   []
# Clean up the output so people can use it  []
# Put these benchmarks in benchmark module  []
# New saving routine?   []
# Unit tests!

#-----------------------------------------------------------------------------




def RMSE(list1, list2):
    return np.sqrt(np.mean((list1-list2)**2))

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
            
            elif solver.sigma_func[4] == 1:
                plt.figure(4)
                mu2 = self.mus[np.argmin(np.abs(self.mus - 0.6))]
                # plt.plot(self.xs[-1], self.psi[-1,-1], '--', label = 'mu = 1 from solver')
                # plt.plot(self.xs[-1], self.psi[-1,np.argmin(np.abs(self.mus - 0.6))], '--', label = f'mu = {round(mu2,2)} from solver')
                plt.plot(self.xs[-1], self.phi[-1])
                plt.legend()
                plt.show()
                self.fake_sedov_benchmark()
                print(self.mus[-1], 'last mu')

                plt.figure(5)
                plt.plot(self.eval_array, self.exit_phi[:,0], '--b', mfc = 'none', label = 'left exit distribution')
                plt.xlabel('t')
                plt.legend()
                plt.show()
                # plt.plot(self.exit_dist[-1], "right exit")
                # plt.show()

                plt.figure(6)
                plt.plot(self.eval_array, self.exit_phi[:,-1], '--b', mfc = 'none', label = f'right exit distribution')
                plt.legend()
                plt.xlabel('t')
                plt.show()
                # plt.plot(self.exit_dist[0], "left exit")



    def dipole(self, uncollided = True, moving = True, All = False):
        plt.ion()
        # plt.figure(1)
        source_name = "dipole"
        print("---  ---  ---  ---  ---  ---  ---")
        print("running dipole")
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
    def self_sim_plane(self, uncollided = True, moving = True, All = False):
        plt.ion()
        # plt.figure(1)
        source_name = "self_sim_plane"
        print("---  ---  ---  ---  ---  ---  ---")
        print("running self_sim_plane")
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
        self.v0 = solver.fake_sedov_v0
        self.t0_source = solver.t0
        self.eval_array = solver.eval_array
        self.exit_phi = solver.exit_phi

        
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
    
    def moving_gaussian_shock(self, tfinal, x, mu):
        sigma_v = -0.25
        
        left_bound = - 5.0
        x0 = -5 
        t0 =  (x0-x)/mu + tfinal # time the particle is emitted
        # right_bound = x - sigma_v * (tfinal - t0)
        right_bound = x

        eta_func = lambda s, mu: (s - x0) / mu + t0
        sigma_func = lambda s, t: np.exp(-(s-sigma_v * t)**2/(2*2**2))
        integrand = lambda s, t, mu: sigma_func(s, eta_func(s, mu))

        mfp = integrate.quad(integrand, left_bound, right_bound, args = (tfinal, mu))[0]

        return np.exp(-mfp / abs(mu)) * np.heaviside(abs(mu) - abs((x - x0)/ tfinal), 0)

    def fake_sedov(self, mu, tfinal, x):
        c1 = 1.0
        v0 = self.v0
        x0 = -5 
        t0 =  (x0-x)/mu + tfinal # time the particle is emitted
        x02 = 0.0
        sqrt_pi = math.sqrt(math.pi)
        kappa = 2
        rho0 = 0.2
        # beta = c1 * (v0-1) - v0 * (x0/mu + t0)
        
        # b2 =  v0 * (-x0/mu - t0 + c1) / (1+v0/mu)
        b2 = ((v0*x0) - t0*v0*mu)/(v0 + mu)
        b1 = max(x, b2)
        # b2 = 0

        b4 = x0
        # b3 = min(x,0)

        b3 =  min(x, b2)

        # print(b1, b2, b3, b4, 'bs', x, 'x', t0, 't0')

        # t1 = lambda s: -0.5*(mu*sqrt_pi*kappa*erf((beta - (s*(mu + v0))/mu)/kappa))/(mu + v0)
        t1 = lambda s: (sqrt_pi*kappa*mu*erf((v0*(s - x0) + (c1 + s + t0*v0)*mu)/(kappa*mu)))/(2.*(v0 + mu))
        t2 = lambda s: rho0 * s

        mfp = t1(b1) - t1(b2) + t2(b3) - t2(b4) 
        # mfp = rho0 * x - rho0 * (-x0)
        # print(mfp, x, 'mfp')

        return np.exp(-mfp / mu) * np.heaviside(mu - abs(x - x0)/ (tfinal), 0) * np.heaviside(abs(x0-x) - (tfinal-self.t0_source)*mu,0)

 

    def moving_gaussian_shock_benchmark(self):
        psi1 = np.zeros(25)
        psi2 = np.zeros(25)
        mu1 = 1
        mu2 = self.mus[np.argmin(np.abs(self.mus - 0.6))]
        sparse_xs = np.linspace(self.xs[0], self.xs[-1], 25)
        for ix, xx in enumerate(sparse_xs):
            psi1[ix] = self.moving_gaussian_shock(self.tfinal, xx, mu1)
            psi2[ix] = self.moving_gaussian_shock(self.tfinal, xx, mu2)

        plt.figure(4)
        plt.plot(sparse_xs, psi1,'o', mfc = 'none', label = 'benchmark mu = 1')
        plt.plot(sparse_xs, psi2,'o', mfc = 'none', label = f'benchmark mu = {round(mu2,2)}')
        plt.legend()
        plt.xlabel('x', fontsize = 16)
        plt.ylabel(r'$\psi$', fontsize = 16)
        plt.show()

    def fake_sedov_benchmark(self):
        psi1 = np.zeros(50)
        psi2 = np.zeros(50)
        mu1 = 1
        mu2 = self.mus[np.argmin(np.abs(self.mus - 0.6))]
        sparse_xs = np.linspace(self.xs[-1,0], self.xs[-1,-1], 50)

        for ix, xx in enumerate(sparse_xs):
            psi1[ix] = self.fake_sedov(mu1, self.tfinal, xx)
            psi2[ix] = self.fake_sedov(mu2, self.tfinal, xx)


        psi_all = np.zeros((sparse_xs.size, self.mus.size))
        for imu, mu in enumerate(self.mus):
            for ix, xx in enumerate(sparse_xs):
                psi_all[ix, imu] = self.fake_sedov(mu, self.tfinal, xx)
        phi = np.zeros(sparse_xs.size)
        for ix in range(sparse_xs.size):
            # phi[ix] = integrate.quad(self.fake_sedov, -1, 1, args = (self.tfinal, xx))[0]
            # print(phi[ix])
            phi[ix] = np.sum(self.ws * psi_all[ix, :])
        
        plt.figure(4)
        # plt.plot(sparse_xs, psi1,'o', mfc = 'none', label = 'benchmark mu = 1')
        # plt.plot(sparse_xs, psi2,'o', mfc = 'none', label = f'benchmark mu = {round(mu2,2)}')
        plt.plot(sparse_xs, phi, 's', mfc = 'none', label = r'$\phi$ benchmark ' + f't = {self.tfinal}')
        print(phi)
        print('--- --- --- --- --- --- --- ')
        print(psi_all[:,:])
        plt.legend()
        plt.xlabel('x', fontsize = 16)
        plt.ylabel(r'$\phi$', fontsize = 16)
        plt.show()

        print("#--- --- --- --- --- --- --- --- ---#")
        phi_interp = interp.interp1d(self.xs[-1], self.phi)
        phi_eval = phi_interp(sparse_xs)
        print('RMSE', RMSE(phi_eval, phi))
        print()



