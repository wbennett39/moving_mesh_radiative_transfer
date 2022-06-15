#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 24 13:42:55 2022

@author: bennett
"""
import matplotlib.pyplot as plt
from .solver import main_class

def run_plane_IC(uncollided = True, moving = True):
    plt.ion()
    plt.figure(1)
    source_name = "plane_IC"
    print("---  ---  ---  ---  ---  ---  ---")
    print("running plane IC")
    print("---  ---  ---  ---  ---  ---  ---")
    solver = main_class(source_name) 
    solver.main(uncollided, moving)
    plt.title("plane IC")
    plt.legend()
    plt.show(block = False)
    
def run_square_IC(uncollided = True, moving = True):
    plt.ion()
    plt.figure(2)
    source_name = "square_IC"
    print("---  ---  ---  ---  ---  ---  ---")
    print("running square IC")
    print("---  ---  ---  ---  ---  ---  ---")
    solver = main_class(source_name) 
    solver.main(uncollided, moving)
    plt.title("square IC")
    plt.legend()
    plt.show(block = False)
    
def run_square_source(uncollided = True, moving = True):
    plt.ion()
    plt.figure(3)
    source_name = "square_source"
    print("---  ---  ---  ---  ---  ---  ---")
    print("running square source")
    print("---  ---  ---  ---  ---  ---  ---")
    solver = main_class(source_name) 
    solver.main(uncollided, moving)
    plt.title("square source")
    plt.legend()
    plt.show(block = False)
    
def run_gaussian_IC(uncollided = True, moving = True):
    plt.ion()
    plt.figure(4)
    source_name = "gaussian_IC"
    print("---  ---  ---  ---  ---  ---  ---")
    print("running Gaussian IC")
    print("---  ---  ---  ---  ---  ---  ---")
    solver = main_class(source_name) 
    solver.main(uncollided, moving)
    plt.title("Gaussian IC")
    plt.legend()
    plt.show(block = False)
    
def run_gaussian_source(uncollided = True, moving = True):
    plt.ion()
    plt.figure(5)
    source_name = "gaussian_source"
    print("---  ---  ---  ---  ---  ---  ---")
    print("running Gaussian source")
    print("---  ---  ---  ---  ---  ---  ---")
    solver = main_class(source_name) 
    solver.main(uncollided, moving)
    plt.title("Gaussian source")
    plt.legend()
    plt.show(block = False)
    
def run_MMS(uncollided = False, moving = True):
    plt.ion()
    plt.figure(6)
    source_name = "MMS"
    print("---  ---  ---  ---  ---  ---  ---")
    print("running MMS problem")
    print("---  ---  ---  ---  ---  ---  ---")
    solver = main_class(source_name) 
    solver.main(uncollided, moving)
    plt.title("MMS")
    plt.legend()
    plt.show(block = False)
    
def run_all():
    # run_plane_IC(True, True)
    # run_plane_IC(True, False)
    # # # run_plane_IC(False, True)        # this doesn't converge
    # run_plane_IC(False, False)
    
    run_square_IC(True, True)
    run_square_IC(True, False)
    run_square_IC(False, True)
    run_square_IC(False, False)
    
    # run_square_source(True, True)
    # run_square_source(True, False)
    # run_square_source(False, True)
    # run_square_source(False, False)
    
    # run_gaussian_IC(True, True)
    # run_gaussian_IC(True, False)
    # run_gaussian_IC(False, True)
    # run_gaussian_IC(False, False)
    
    # run_gaussian_source(True, True)
    # run_gaussian_source(True, False)
    # run_gaussian_source(False, True)
    # run_gaussian_source(False, False)
    
    # run_MMS(False, True)            # only one case is possible for the MMS