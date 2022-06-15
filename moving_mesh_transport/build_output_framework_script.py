#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 09:55:54 2022

@author: bennett
"""

import h5py


def make_RMS_file():
    
    f = h5py.File('run_data_RMS.h5', 'a')
    plane = f.create_group("plane_IC")
    square_ic = f.create_group("square_IC")
    square_s = f.create_group("square_s")
    trunc_gauss_ic = f.create_group("gaussian_IC")
    trunc_gauss_s = f.create_group("gaussian_s")
    MMS = f.create_group("MMS")
    
    plane.create_group("t=1/RMS")
    plane.create_group("t=5/RMS")
    plane.create_group("t=10/RMS")
    
    square_ic.create_group("t=1/RMS")
    square_ic.create_group("t=5/RMS")
    square_ic.create_group("t=10/RMS")
    
    square_s.create_group("t=1/RMS")
    square_s.create_group("t=5/RMS")
    square_s.create_group("t=10/RMS")
    
    trunc_gauss_ic.create_group("t=1/RMS")
    trunc_gauss_ic.create_group("t=5/RMS")
    trunc_gauss_ic.create_group("t=10/RMS")
    
    trunc_gauss_s.create_group("t=1/RMS")
    trunc_gauss_s.create_group("t=5/RMS")
    trunc_gauss_s.create_group("t=10/RMS")
    
    MMS.create_group("t=1/RMS")
    MMS.create_group("t=5/RMS")
    MMS.create_group("t=10/RMS")
    
    f.close()
    
    f = h5py.File('run_data_radiative_transfer_RMS.h5', 'a')
    square_s = f.create_group("gaussian_s")
    
    square_s.create_group("t=1/RMS")
    square_s.create_group("t=5/RMS")
    square_s.create_group("t=10/RMS")
    square_s.create_group("t=31.6228/RMS")
    
    square_s.create_group("t=1/RMS/S2")
    square_s.create_group("t=5/RMS/S2")
    square_s.create_group("t=10/RMS/S2")
    square_s.create_group("t=31.6228/RMS/S2")
    
    
    
    f.close()
    
def make_solution_file():
    
    f = h5py.File('run_data.h5', 'a')
    
    transport = f.create_group("transport")
    rad_transfer = f.create_group("transfer")
    
    
    f.close()
    

make_solution_file()


