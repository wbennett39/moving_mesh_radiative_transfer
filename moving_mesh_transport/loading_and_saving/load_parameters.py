#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 10:49:08 2022

@author: bennett
"""

import numpy as np
import math


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
        self.scattering_ratio = float(parameters['all']['c'])
        self.major = str(parameters['all']['major'])
        self.thermal_couple = int(parameters['all']['radiative_transfer'])
        self.temp_function = np.array(parameters['all']['temperature_dependence'])
        self.e_initial = float(parameters['all']['e_initial'])
        self.weights = str(parameters['all']['weights'])
        self.particle_v = str(parameters['all']['particle_v'])
        self.edge_v = str(parameters['all']['edge_v'])
        self.cv0 = float(parameters['all']['cv_const'])
        self.problem_type = str(parameters['all']['problem_name'])

        self.thick = int(parameters['all']['thick'])
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
        self.x0 = np.array(parameters[source_name]['x0'])
        self.source_type = np.array(parameters[source_name]['source_type'])
        if self.source_type[0] == 1:
            self.x0 = np.array([1e-14])
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
        self.estimate_wavespeed = int(mesh_parameters['estimate_wavespeed'])
        self.find_wave_loc = int(mesh_parameters['find_wave_loc'])
        self.choose_xs = int(mesh_parameters['choose_xs'])
        self.specified_xs = np.array(mesh_parameters['xs_list'])


        if not (len(self.N_spaces) == len(self.N_angles) == len(self.Ms)):
            print('Spaces, Ms, and N_angles should be the same length')
            assert(0)




            
