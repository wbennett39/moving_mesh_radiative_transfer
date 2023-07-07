#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 09:14:26 2022

@author: bennett
"""
import matplotlib.pyplot as plt
from .make_plots import rms_plotter
from ..loading_and_saving.load_solution import load_sol
from ..solver_classes.functions import  convergence
from .plot_functions.coeff_con import coeff_con
import numpy as np

def plot_all_rms_cells(tfinal, M):
    major = 'cells'
    # source_type_list  = ["plane_IC", "square_IC", "square_s", "gaussian_IC", "gaussian_s",
    #                       "su_olson", "su_olson_energy", "su_olson_s2", "su_olson_energy_s2",
    #                       "gaussian_s2", "gaussian_energy_s2"]
    
    source_type_list  = ["plane_IC", "square_IC", "square_s", "gaussian_IC", "gaussian_s"
          , "su_olson_s2", "su_olson_energy_s2", "gaussian_s2", "gaussian_energy_s2",
          "gaussian_s_thick_s2", "gaussian_s_thick_s2", "su_olson_thick_s2", "su_olson_thick_s8" ]
    
    case_list_1 = [True, True, False, False]
    case_list_2 = [True, False, True, False]
    
    for count1, source in enumerate(source_type_list):
        print(source)
        plotter = rms_plotter(tfinal, M, source, major)
        print("loaded")
        for count2, uncollided in enumerate(case_list_1):
            moving = case_list_2[count2]
            if source == "plane_IC" and uncollided == False and moving == True:
                print("skipping no-uncol moving case for plane IC")
            else:
                
                plotter.load_RMS_data(uncollided, moving)
                if count2 == 0:
                    clear = False
                plotter.plot_RMS_vs_cells(count1+1, clear)
                
    plotter = rms_plotter(tfinal, 2, "MMS", major)
    plotter.load_RMS_data(uncollided = False, moving = True)
    plotter.plot_RMS_vs_cells(25, clear)
    
    plotter = rms_plotter(tfinal, 4, "MMS", major)
    plotter.load_RMS_data(uncollided = False, moving = True)
    plotter.plot_RMS_vs_cells(25, clear)
    
    plotter = rms_plotter(tfinal, 6, "MMS", major)
    plotter.load_RMS_data(uncollided = False, moving = True)
    plotter.plot_RMS_vs_cells(25, clear)

def plot_rms_Ms(tfinal, source_name, fign = 1 ):
    major = 'Ms'
    case_list_1 = [True, True, False, False]
    case_list_2 = [True, False, True, False]
    
    plotter = rms_plotter(tfinal, 1,  source_name, major)

    if source_name == "MMS":
        plotter.load_RMS_data(False, True)
        plotter.plot_RMS_vs_Ms(source_name, fign, False)

    elif source_name  in [ 'su_olson_thick_s2','su_olson_thick_s8','su_olson_thick_s2_energy','su_olson_thick_s8_energy', 
    'gaussian_s_thick_s2', 'gaussian_s_thick_s2_energy']:
        plotter.load_RMS_data(True, False)
        plotter.plot_RMS_vs_Ms(source_name, fign, False)
        plotter.load_RMS_data(False, False)
        plotter.plot_RMS_vs_Ms(source_name, fign, False)
    else:

        for count2, uncollided in enumerate(case_list_1):
            moving = case_list_2[count2]
            plotter.load_RMS_data(uncollided, moving)
            plotter.plot_RMS_vs_Ms(source_name, fign, False)

        
def plot_rms_cells(tfinal, source_name, M, fign = 1):
    major = 'cells'
    case_list_1 = [True, True, False, False]
    case_list_2 = [True, False, True, False]
    
    plotter = rms_plotter(tfinal, M,  source_name, major)
    if source_name == "MMS":
        plotter.load_RMS_data(False, True)
        plotter.plot_RMS_vs_cells()

    elif source_name in ['su_olson_thick_s2', 'su_olson_thick', 'su_olson_thick_s2_energy']:
        plotter.load_RMS_data(True, False)
        plotter.plot_RMS_vs_cells()
        plotter.load_RMS_data(False, False)
        plotter.plot_RMS_vs_cells()
    else:
        for count2, uncollided in enumerate(case_list_1):
            moving = case_list_2[count2]
            plotter.load_RMS_data(uncollided, moving)
            plotter.plot_RMS_vs_cells()
        
        
def plot_all_rms_Ms(tfinal):
    # plot_rms_Ms(tfinal, "plane_IC")
    # plot_rms_Ms(tfinal, "square_IC", 1)
    # plot_rms_Ms(tfinal, "square_s", 2)
    plot_rms_Ms(tfinal, "gaussian_IC", 3)
    plot_rms_Ms(tfinal, "gaussian_s", 4)
    plot_rms_Ms(tfinal, "MMS", 5)
    
    
def plot_best():
    tfinal = 1
    M = 6
    major = 'cells'
    plotter = rms_plotter(tfinal, M, 'plane_IC', major)
    plotter.load_RMS_data(True, True)
    plotter.plot_best('--o', 'darkorange')
    
    plotter = rms_plotter(tfinal, M, 'square_IC', major)
    plotter.load_RMS_data(True, True)
    plotter.plot_best('--s', 'r')
    
    plotter = rms_plotter(tfinal, M, 'square_s', major)
    plotter.load_RMS_data(True, True)
    plotter.plot_best('-s', 'r')
    
    plotter = rms_plotter(tfinal, M, 'gaussian_IC', major)
    plotter.load_RMS_data(True, True)
    plotter.plot_best('--p', 'b')
    
    plotter = rms_plotter(tfinal, M, 'gaussian_s', major)
    plotter.load_RMS_data(True, True)
    plotter.plot_best('-p', 'b')
    
    # plotter = rms_plotter(tfinal, M, 'MMS', major)
    # plotter.load_RMS_data(False, True)
    # plotter.plot_best('-^', 'm')
    
def plot_best_static():
    tfinal = 1
    M = 6
    major = 'cells'
    plotter = rms_plotter(tfinal, M, 'plane_IC', major)
    plotter.load_RMS_data(True, False)
    plotter.plot_best_static('--o', 'darkorange')
    
    plotter = rms_plotter(tfinal, M, 'square_IC', major)
    plotter.load_RMS_data(True, False)
    plotter.plot_best_static('--s', 'r')
    
    plotter = rms_plotter(tfinal, M, 'square_s', major)
    plotter.load_RMS_data(True, False)
    plotter.plot_best_static('-s', 'r')
    
    plotter = rms_plotter(tfinal, M, 'gaussian_IC', major)
    plotter.load_RMS_data(True, False)
    plotter.plot_best_static('--p', 'b')
    
    plotter = rms_plotter(tfinal, M, 'gaussian_s', major)
    plotter.load_RMS_data(True, False)
    plotter.plot_best_static('-p', 'b')
    
    
    
    
    
    
def plot_all_rms_times(tfinal, M, major):
    source_type_list  = ["plane_IC", "square_IC", "square_s", "gaussian_IC", "gaussian_s"]
    case_list_1 = [True, True, False, False]
    case_list_2 = [True, False, True, False]
    
    for count1, source in enumerate(source_type_list):
        plotter = rms_plotter(tfinal, M, source, major)
        for count2, uncollided in enumerate(case_list_1):
            moving = case_list_2[count2]
            plotter.load_RMS_data(uncollided, moving)
            if count2 == 0:
                clear = False
            if source == "plane_IC" and uncollided == False and moving == True:
                print("skipping no-uncol moving case for plane IC")
            else:
                plotter.plot_RMS_vs_times(count1+1, clear)
                
    plotter = rms_plotter(tfinal, 2, "MMS")
    plotter.load_RMS_data(uncollided = False, moving = True)
    plotter.plot_RMS_vs_times(6, clear)
    
    plotter = rms_plotter(tfinal, 4, "MMS")
    plotter.load_RMS_data(uncollided = False, moving = True)
    plotter.plot_RMS_vs_times(6, clear)
    
    plotter = rms_plotter(tfinal, 6, "MMS")
    plotter.load_RMS_data(uncollided = False, moving = True)
    plotter.plot_RMS_vs_times(6, clear)
    
def plot_rms_Ms_times(tfinal, source_name, fign):
    major = 'Ms'
    case_list_1 = [True, True, False, False]
    case_list_2 = [True, False, True, False]
    
    plotter = rms_plotter(tfinal, 1,  source_name, major)
    if source_name != "MMS":
        for count2, uncollided in enumerate(case_list_1):
            moving = case_list_2[count2]
            plotter.load_RMS_data(uncollided, moving)
            plotter.plot_RMS_vs_times()
    else:
        plotter.load_RMS_data(False, True)
        plotter.plot_RMS_vs_Ms(source_name, 1, False)
    
def compare_rms(tfinal, M):
    source_type_list  = ["plane_IC", "square_IC", "square_s", "gaussian_IC", "gaussian_s"]
    case_list_1 = [True, True, False, False]
    case_list_2 = [True, False, True, False]
    
    for count1, source in enumerate(source_type_list):
        plotter = rms_plotter(tfinal, M, source)
        plotter.load_RMS_data(uncollided = True, moving = True)
        RMS_list_case_1 = plotter.RMS
        if source != 'plane_IC':
            plotter.load_RMS_data(uncollided = False, moving = True)
            RMS_list_case_2 = plotter.RMS
        plotter.load_RMS_data(uncollided = True, moving = False)
        RMS_list_case_3 = plotter.RMS
        plotter.load_RMS_data(uncollided = False, moving = False)
        RMS_list_case_4 = plotter.RMS
        diff2 = 0
        diff3 = 0
        diff4 = 0
        for i in range(RMS_list_case_1.size): ## calculate percent error
            if source != 'plane_IC':
                diff2 += (RMS_list_case_2[i] - RMS_list_case_1[i])/RMS_list_case_1[i]
            diff3 += (RMS_list_case_3[i] - RMS_list_case_1[i])/RMS_list_case_1[i]
            diff4 += (RMS_list_case_4[i] - RMS_list_case_1[i])/RMS_list_case_1[i]
        diff2 = diff2/RMS_list_case_1.size
        diff3 = diff3/RMS_list_case_1.size
        diff4 = diff4/RMS_list_case_1.size
        
        print("--------------")
        print(source)
        print("case 1 is moving with uncollided")
        print(diff2, "average percent difference between case 1 and moving w/o uncollided ")
        print(diff3, "average percent difference between case 1 and static w/ uncollided ")
        print(diff4, "average percent difference between case 1 and static w/o uncollided ")
        
def plot_all_benchmarks(tfinal):
    M = 3
    if tfinal <100:
        # source_list = ["plane_IC", "square_IC", "square_source", "gaussian_IC",
        # "gaussian_source", "MMS", "gaussian_IC_2D", "line_source", "square_IC_c_not_one", 
        # "gaussian_IC_c_not_one", 'square_source_s2','gaussian_source_s2', 'square_source_s2_thick', 'gaussian_source_s2_thick' ]
        # source_list = ['square_source_s2','gaussian_source_s2']
        source_list = ['square_IC', 'plane_IC']
    elif tfinal >=100:
        source_list = ['square_source_s2_thick', 'gaussian_source_s2_thick']

    for count, source in enumerate(source_list):
        print(source)
        plotter = rms_plotter(tfinal, M, source, "cells")
        plotter.plot_bench(tfinal, source, count)
        
def plot_coefficients(tfinals = [3.0],  Ms=[8], source_name = 'gaussian_s',  N_spaces = [16], problem_name = 'transfer_const_cv=0.03_thick',
rad_or_transport ='transfer', x0_or_sigma = 0.375 , c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = False, s2 = False, moving = False, line = '-',legend = True, fign = 1):
    
    ifshow = False

    for count, tfinal in enumerate(tfinals):

        if tfinal == 100.0 or tfinal == 30.0:
            ifshow = True
        # ifshow = True
        plotter = rms_plotter(tfinal, Ms[count], source_name, 'cells')
        plotter.plot_coefficients(tfinal,  Ms[count], source_name,  [N_spaces[count]], problem_name, rad_or_transport,
        x0_or_sigma, c, cv0, uncollided, s2, mat_or_rad, moving, line, legend, count, fign, ifshow)

        plotter.plot_coeff_boyd()

        plotter = rms_plotter(tfinal, Ms[count], source_name, 'cells')
        plotter.plot_coefficients(tfinal,  Ms[count], source_name,  [N_spaces[count]], problem_name, rad_or_transport,
        x0_or_sigma, c, cv0, uncollided, s2, 'mat', moving, line, legend, count, 3, ifshow)

        plotter.plot_coeff_boyd()

    
    
    
    
    
    
        
