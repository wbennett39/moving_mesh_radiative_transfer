#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 24 13:42:55 2022

@author: bennett
"""
import matplotlib.pyplot as plt
from .solver import main_class
from pathlib import Path
import yaml

import argparse
parser = argparse.ArgumentParser(prog='PROG', description='description')
parser.add_argument('cmd', choices=['transport','rad_transfer','s2_rad_transfer',
                                    'rad_transfer_thick','config', 'help','quit'])

print('choose problem type: [transport, rad_transfer, s2_rad_transfer, rad_transfer_thick, config]')
while True:
    astr = input('$: ')
    # print astr
    try:
        args = parser.parse_args(astr.split())
    except SystemExit:
        # trap argparse error message
        print('error')
        continue
    if args.cmd in ['transport','rad_transfer','s2_rad_transfer',
                                        'rad_transfer_thick','config']:
        print('loading', args.cmd)
        data_folder = Path("moving_mesh_transport")
        config_file_path = data_folder / f"{args.cmd}.yaml"
        with open(config_file_path, 'r') as file:
            parameters = yaml.safe_load(file)
            file.close()
        break
    elif args.cmd == 'help':
        parser.print_help()
    else:
        print('error')
        break
    

###############################################################################

###############################################################################

def run_plane_IC(uncollided = True, moving = True, All = False):
    plt.ion()
    plt.figure(1)
    source_name = "plane_IC"
    print("---  ---  ---  ---  ---  ---  ---")
    print("running plane IC")
    print("---  ---  ---  ---  ---  ---  ---")
    
    solver = main_class(source_name, parameters) 
    if All == True:
        solver.main(True, True)
        solver.main(False, True)
        solver.main(True, False)
        solver.main(False, False)
    else:
        solver.main(uncollided, moving)
    plt.title("plane IC")
    plt.legend()
    plt.show(block = False)
    
def run_square_IC(uncollided = True, moving = True, All = False):
    plt.ion()
    plt.figure(2)
    source_name = "square_IC"
    print("---  ---  ---  ---  ---  ---  ---")
    print("running square IC")
    print("---  ---  ---  ---  ---  ---  ---")
    solver = main_class(source_name, parameters) 
    if All == True:
        solver.main(True, True)
        solver.main(False, True)
        solver.main(True, False)
        solver.main(False, False)
    else:
        solver.main(uncollided, moving)
    plt.title("square IC")
    plt.legend()
    plt.show(block = False)
    
def run_square_source(uncollided = True, moving = True, All = False):
    plt.ion()
    plt.figure(3)
    source_name = "square_source"
    print("---  ---  ---  ---  ---  ---  ---")
    print("running square source")
    print("---  ---  ---  ---  ---  ---  ---")
    solver = main_class(source_name, parameters) 
    if All == True:
        solver.main(True, True)
        solver.main(False, True)
        solver.main(True, False)
        solver.main(False, False)
    else:
        solver.main(uncollided, moving)
    plt.title("square source")
    plt.legend()
    plt.show(block = False)
    
def run_gaussian_IC(uncollided = True, moving = True, All = False):
    plt.ion()
    plt.figure(4)
    source_name = "gaussian_IC"
    print("---  ---  ---  ---  ---  ---  ---")
    print("running Gaussian IC")
    print("---  ---  ---  ---  ---  ---  ---")
    solver = main_class(source_name, parameters) 
    if All == True:
        solver.main(True, True)
        solver.main(False, True)
        solver.main(True, False)
        solver.main(False, False)
    else:
        solver.main(uncollided, moving)
    plt.title("Gaussian IC")
    plt.legend()
    plt.show(block = False)
    
def run_gaussian_source(uncollided = True, moving = True, All = False):
    plt.ion()
    plt.figure(5)
    source_name = "gaussian_source"
    print("---  ---  ---  ---  ---  ---  ---")
    print("running Gaussian source")
    print("---  ---  ---  ---  ---  ---  ---")
    solver = main_class(source_name, parameters) 
    if All == True:
        solver.main(True, True)
        solver.main(False, True)
        solver.main(True, False)
        solver.main(False, False)
    else:
        solver.main(uncollided, moving)
    plt.title("Gaussian source")
    plt.legend()
    plt.show(block = False)
    
def run_MMS(uncollided = False, moving = True, All = False):
    plt.ion()
    plt.figure(6)
    source_name = "MMS"
    print("---  ---  ---  ---  ---  ---  ---")
    print("running MMS problem")
    print("---  ---  ---  ---  ---  ---  ---")
    solver = main_class(source_name, parameters) 
    if All == True:
        solver.main(True, True)
        solver.main(False, True)
        solver.main(True, False)
        solver.main(False, False)
    else:
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