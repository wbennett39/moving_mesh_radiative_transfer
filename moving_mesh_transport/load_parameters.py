#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 10:49:08 2022

@author: bennett
"""
import yaml
from pathlib import Path
import numpy as np

###############################################################################
data_folder = Path("moving_mesh_transport")
config_file_path = data_folder / "config.yaml"
###############################################################################



class parameter_load_class:
    def __init__(self, source_name):
        with open(config_file_path, 'r') as file:
           parameters = yaml.safe_load(file)
        file.close()
        #'C:/users/bennett/Documents/GitHub/MovingMesh/moving_mesh_radiative_transfer/src/package/
        self.tfinal = float(parameters['all']['tfinal'])
        self.N_spaces = np.array(parameters['all']['N_spaces'])
        self.Ms = np.array(parameters['all']['Ms'])
        self.N_runs = int(parameters['all']['N_runs'])
        self.t_nodes = int(parameters['all']['tnodes'])
        self.rt = float(parameters['all']['rt'])
        self.at = float(parameters['all']['at'])
        self.t0 = float(parameters['all']['t0'])
        self.scattering_ratio = float(parameters['all']['c'])
        self.major = str(parameters['all']['major'])
        self.thermal_couple = int(parameters['all']['radiative_transfer'])
        self.temp_function = np.array(parameters['all']['temperature_dependence'])
        self.e_initial = float(parameters['all']['e_initial'])
        self.weights = str(parameters['all']['weights'])
        self.benchmarking = int(parameters['all']['benchmarking'])
        self.saving = int(parameters['all']['save_solution'])
        
        
        self.N_angles = np.array(parameters[source_name]['N_angles'])
        self.x0 = np.array(parameters[source_name]['x0'])
        self.source_type = np.array(parameters[source_name]['source_type'])
        self.move_type = np.array(parameters[source_name]['move_type'])
        
        self.r_times = np.zeros(len(self.N_angles))
        self.RMS_list = np.zeros(len(self.N_angles))
        self.RMS_list_energy = np.zeros(len(self.N_angles))
        
        
        if source_name in ["gaussian_IC", "gaussian_source"]:
            self.sigma = float(parameters[source_name]['sigma'])
            self.x0_or_sigma = self.sigma
        else:
            self.x0_or_sigma = self.x0
            self.sigma = 0