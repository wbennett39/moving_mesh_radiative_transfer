#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 07:47:25 2022

@author: bennett
"""
import numpy as np

def coeff_con(ws, sol_matrix, N_ang, M, k):
    
    weight_avg = np.zeros(M+1)
    
    for j in range(M+1):
    
        weight_avg[j] = np.sum(ws*sol_matrix[:,k,j])/np.sum(ws)
    
    return weight_avg
    
    
    
    