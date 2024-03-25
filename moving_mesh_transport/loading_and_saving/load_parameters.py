#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 10:49:08 2022

@author: bennett
"""

import numpy as np
import math
import numba as nb



class parameter_load_class:
    def __init__(self, source_name, parameters, mesh_parameters):
        # with open(config_file_path, 'r') as file:
        #    parameters = yaml.safe_load(file)
        # file.close()
        #'C:/users/bennett/Documents/GitHub/MovingMesh/moving_mesh_radiative_transfer/src/package/
        self.tfinal = float(parameters['all']['tfinal'])
        self.N_spaces = np.array(parameters['all']['N_spaces'])
        self.Ms = np.array(parameters['all']['Ms'])
        self.N_runs = int(parameters['all']['N_runs'])
        self.t_nodes = int(parameters['all']['tnodes'])
        self.c_scaling = int(parameters['all']['c_scaling'])
        self.rt = float(parameters['all']['rt'])
        self.at = float(parameters['all']['at'])
        self.t0 = float(parameters['all']['t0'])
        # self.scattering_ratio = float(parameters['all']['c'])
        self.major = str(parameters['all']['major'])
        # self.thermal_couple = int(parameters['all']['radiative_transfer'])
        self.thermal_couple = nb.typed.Dict.empty(key_type=nb.typeof('par_1'), value_type=nb.typeof(1))
        dictionary_loader(parameters['all']['radiative_transfer'], self.thermal_couple)   
        self.temp_function = np.array(parameters['all']['temperature_dependence'])
        self.e_initial = float(parameters['all']['e_initial'])
        self.weights = str(parameters['all']['weights'])
        self.particle_v = str(parameters['all']['particle_v'])
        self.edge_v = str(parameters['all']['edge_v'])
        self.cv0 = float(parameters['all']['cv_const'])
        self.problem_type = str(parameters['all']['problem_name'])
        self.sigma_t = float(parameters['all']['sigma_t'])
        self.sigma_s = float(parameters['all']['sigma_s'])
        self.scattering_ratio = self.sigma_s/self.sigma_t
        self.integrator = str(parameters['all']['integrator'])
        self.epsilon = float(parameters['all']['epsilon'])
        self.geometry = nb.typed.Dict.empty(key_type=nb.typeof('par_1'), value_type=nb.typeof(1))
        dictionary_loader(parameters['all']['geometry'], self.geometry)   
        

        self.thick = int(parameters['all']['thick'])
        if self.thick == True:
            self.l = float(parameters[source_name]['l'])
        else:
            self.l = 1.0

        self.mxstp = float(parameters['all']['mxstp'])


        if self.edge_v == 'one':
            self.edge_v = 1.0
        elif self.edge_v == 'sqrt_3':
            self.edge_v = 1.0/math.sqrt(3)
        if self.particle_v == 'one':
            self.particle_v = 1.0 
        elif self.particle_v == 'sqrt_3':
            self.particle_v = 1.0/math.sqrt(3)

        

        self.saving = int(parameters['all']['save_solution'])
    
        self.N_angles = np.array(parameters[source_name]['N_angles'])
        
        self.source_type = (np.array(parameters[source_name]['source_type']))
        
        for iele, ele, in enumerate(self.source_type):
            self.source_type[iele] = np.int64(self.source_type[iele])
        self.source_strength = float(parameters[source_name]['source_strength'])
        
        if self.source_type[0] != 0:
            self.x0 = np.zeros(np.array(parameters[source_name]['x0']).size)
            for ix, xx in enumerate(self.x0):
                self.x0[ix] = float(np.array(parameters[source_name]['x0'])[ix])

        else:
            self.x0 = np.array(parameters[source_name]['x0'])
        self.move_type = np.array(parameters[source_name]['move_type'])
        self.benchmarking = int(parameters[source_name]['benchmarking'])
        self.r_times = np.zeros(len(self.N_angles))
        self.RMS_list = np.zeros(len(self.N_angles))
        self.RMS_list_energy = np.zeros(len(self.N_angles))
        
        
        if source_name in ["gaussian_IC", "gaussian_source"]:
            self.sigma = float(parameters[source_name]['sigma'])
            self.x0_or_sigma = self.sigma
        else:
            self.x0_or_sigma = self.x0
            self.sigma = 0
            
        if source_name in ["square_source", "gaussian_source"]:
            self.bench_type = str(parameters[source_name]['bench_type'])
        else:
            self.bench_type = 'full'


        """ 
        mesh
        """

        self.find_edges_tol = float(mesh_parameters['find_edges_tol'])
        self.move_factor = float(mesh_parameters['sqrt_t_move_factor'])
        self.save_wave_loc = int(mesh_parameters['save_wave_loc'])
        self.find_wave_loc = int(mesh_parameters['find_wave_loc'])
        self.estimate_wavespeed = int(mesh_parameters['estimate_wavespeed'])
        self.choose_xs = int(mesh_parameters['choose_xs'])
        self.specified_xs = np.array(mesh_parameters['xs_list'])
        self.pad = float(mesh_parameters['pad'])
        self.leader_pad = float(mesh_parameters['leader_pad'])
        self.xs_quad = int(mesh_parameters['xs_quad'])
        self.eval_times = int(mesh_parameters['eval_times'])
        self.eval_array = np.array(mesh_parameters['eval_array'])
        self.boundary_on = np.array(mesh_parameters['boundary_on'])
        self.boundary_source = int(mesh_parameters['boundary_source'])
        self.boundary_source_strength = float(mesh_parameters['boundary_source_strength'])
        self.sigma_func = nb.typed.Dict.empty(key_type=nb.typeof('par_1'), value_type=nb.typeof(1))
        # for key in mesh_parameters['sigma_func'].keys():
        #     self.sigma_func[key] = mesh_parameters['sigma_func'][key]
        # 
        dictionary_loader(mesh_parameters['sigma_func'], self.sigma_func)
        print(self.sigma_func['constant'])            
        self.Msigma = int(mesh_parameters['Msigma'])
        self.finite_domain = int(mesh_parameters['finite'])
        self.domain_width = -1
        self.fake_sedov_v0 = float(mesh_parameters['fake_sedov_v0'])
        if self.finite_domain == True:
            self.domain_width = float(mesh_parameters['domain_width'])

        if not (len(self.N_spaces) == len(self.N_angles) == len(self.Ms)):
            print('Spaces, Ms, and N_angles should be the same length')
            assert(0)

        # if (0.82 <self.tfinal <0.84):
        #     self.tfinal = 1./1.2
        # if 0.415 <self.x0_or_sigma <0.417:
        #     self.x0_or_sigma = 0.5/1.2
        #     self.sigma = 0.5/1.2
        if int(parameters['all']['epsilon_scaling']) == True:
            #  self.particle_v = self.particle_v * 29.998
            #  self.sigma_s = self.sigma_s / self.epsilon
            # self.edge_v = self.edge_v / self.epsilon**2
            self.test_dimensional_rhs = True
        else:
            self.test_dimensional_rhs = False



            
def dictionary_loader(inputdict, outputdict):
    counter = 0
    for key in inputdict.keys():
        outputdict[key] = inputdict[key]
        #check if only one element of the list is true
        if inputdict[key] == True: 
            counter += 1
    if counter >= 2:
        print('Two noncomplementary parameters have been selected simultaneously.')
        assert(0)

      
