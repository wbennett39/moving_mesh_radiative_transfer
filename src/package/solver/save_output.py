#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 07:35:06 2022

@author: bennett
"""
import numpy as np

class save_output:
    def __init__(self, build):
        self.M = build.M
        self.tfinal = int(build.tfinal)
        source_name_list = ["plane_IC", "square_IC", "square_source", "truncated_gaussian_IC"]
        index_of_source_name = np.argmin(np.abs(np.array(build.source_type)-1))
        self.source_type = source_name_list[index_of_source_name]
        if self.M == 2:
            self.mkr = "o"
        elif self.M == 4:
            self.mkr = "^"
        elif self.M == 6:
            self.mkr = "s"
        if build.moving == True:
            self.line_mkr = "-"
        elif build.moving == False:
            self.line_mkr = "--"
        if build.uncollided == True:
            self.clr = "c0"
        elif build.uncollided == False:
            self.clr = "c1"
            
    

            
            