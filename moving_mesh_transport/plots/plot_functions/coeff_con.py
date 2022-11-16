#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 07:47:25 2022

@author: bennett
"""
import numpy as np

def coeff_con(ws, sol_matrix, N_ang, M, k):
    
    weight_avg = np.zeros(M+1)


    if len(np.shape(sol_matrix)) == 3:

        for j in range(M+1):

            weight_avg[j] = np.sum(ws*sol_matrix[:,k,j])/np.sum(ws)
        
        return weight_avg
    
    elif len(np.shape(sol_matrix)) == 2:

        for j in range(M+1):

            weight_avg[j] = sol_matrix[k,j]
            max_cell = max(sol_matrix[:,j])
            index_of_max = np.argmin(np.abs(sol_matrix[:,j]-max_cell))
            # print('j=', j, 'cell with highest coefficient', index_of_max, "value", max_cell)
        return weight_avg


    