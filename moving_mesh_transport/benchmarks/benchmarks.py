#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 11:58:28 2022

@author: bennett
"""
import numpy as np

from .benchmark_functions import make_benchmark_file_structure, write_to_file

from .uncollided import uncollided_class
from .collided import collided_class

###############################################################################


class make_benchmark:
    def __init__(self, source_type, x0, t0):
        self.x0 = x0
        self.t0 = t0
        self.source_type = source_type
        
        self.call_uncollided = uncollided_class(source_type, self.x0, self.t0)
        self.call_collided = collided_class(source_type, self.x0, self.t0)
    
    def integrate(self, t, npnts):
        self.t = t
        self.npnts = npnts
        self.xs = np.linspace(-t, t, npnts)
        self.uncollided_sol = self.call_uncollided(self.xs, t)
        self.collided_sol = self.call_collided(self.xs, t)
        
    def save(self):
        phi = self.uncollided_sol + self.collided_sol
        write_to_file(self.xs, phi, self.uncollided_sol, self.t, self.source_name, self.npnts)
    
    def clear_file(self):
        make_benchmark_file_structure()
    
        
        
    
        
        
        
        
        
        
    
    