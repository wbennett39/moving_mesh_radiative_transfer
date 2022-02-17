#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 09:55:54 2022

@author: bennett
"""

import h5py

f = h5py.File('run_data_RMS.h5', 'w')
plane = f.create_group("plane_IC")
square_ic = f.create_group("square_IC")
square_s = f.create_group("square_s")
trunc_gauss_ic = f.create_group("truncated_gaussian_IC")
trunc_gauss_s = f.create_group("truncated_gaussian_s")

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



f.close()




