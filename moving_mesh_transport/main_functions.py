#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 16 07:02:38 2022

@author: bennett
"""
import numpy as np
import scipy.integrate as integrate
import quadpy
import matplotlib.pyplot as plt





def parameter_function(major, N_spaces, Ms, count):
    if major == 'cells':
        M = Ms[0]
        N_space = N_spaces[count]
    elif major == 'Ms':
        N_space = N_spaces[1]
        M = Ms[count]
    return N_space, M


    