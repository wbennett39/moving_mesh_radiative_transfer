#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 15:34:07 2022

@author: bennett
"""
from numba import njit, jit


class IC_func:
    def __init__(self, source, uncollided):
        self.source = source
        self.uncollided = uncollided
    def function(self, x):
        # if self.uncollided == True:
        return 0