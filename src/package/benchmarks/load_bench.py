#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Feb  3 15:12:33 2022

@author: bennett
"""

import h5py 
import numpy as np
from scipy.interpolate import interp1d

###############################################################################

class load_bench:
    def __init__(self, source_type, tfinal):
        self.source_type = source_type
        self.tfinal = tfinal
        f = h5py.File("benchmarks.hdf5", "r")
        
    
    

