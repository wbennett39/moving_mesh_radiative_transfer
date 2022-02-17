#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 09:55:54 2022

@author: bennett
"""

import h5py

f = h5py.File('run_data.h5', 'w')
plane = f.create_group("plane_IC")
square_ic = f.create_group("square_IC")
square_s = f.create_group("square_s")
trunc_gauss_ic = f.create_group("truncated_gaussian_IC")
trunc_gauss_s = f.create_group("truncated_gaussian_s")

f.plane.create_group("t=1")
f.plane.create_group("t=5")
f.plane.create_group("t=10")

f.square_ic.create_group("t=1")
f.square_ic.create_group("t=5")
f.square_ic.create_group("t=10")

f.square_s.create_group("t=1")
f.square_s.create_group("t=5")
f.square_s.create_group("t=10")

f.trunc_gauss_ic.create_group("t=1")
f.trunc_gauss_ic.create_group("t=5")
f.trunc_gauss_ic.create_group("t=10")

f.trunc_gauss_s.create_group("t=1")
f.trunc_gauss_s.create_group("t=5")
f.trunc_gauss_s.create_group("t=10")




f.close()




