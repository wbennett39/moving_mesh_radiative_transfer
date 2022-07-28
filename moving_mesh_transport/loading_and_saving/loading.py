#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 21:05:57 2022

@author: bennett
"""

import h5py 

f = h5py.File("run_data_RMS.h5", "r")
print(f["square_s"]['t=10']['RMS'].keys())
f.close()