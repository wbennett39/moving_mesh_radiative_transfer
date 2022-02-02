#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 09:26:00 2022

@author: bennett
"""
import numpy as np
from build_problem import build
###############################################################################

class scalar_flux:
    def __init__(self, build):
        self.P = np.zeros(build.M+1).transpose()
        self.M = build.M
        self.ws = build.ws
    def __call__(self, u):
        for i in range(0, self.M+1):
            self.P[i]  = np.sum(np.multiply(u[:,i],self.ws))
        return self.P
        
