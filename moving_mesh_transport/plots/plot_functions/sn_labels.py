#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 13:34:59 2022

@author: bennett
"""
import matplotlib.pyplot as plt
import math

def logplot_sn_labels(xdata, ydata, angles_list, dy, fign):
    for count, angle in enumerate(angles_list):
        sn = "$S_{" + str(int(angle)) + "}$";
        xcoord = (xdata[count]);
        ycoord = (ydata[count]) - math.exp(math.log(ydata[count]) - dy);
        plt.text(xcoord, ycoord, sn);        
        
def logplot_sn_labels_2(xdata, ydata, angles_list, dy, fign):
    for count, angle in enumerate(angles_list):
        sn = "$S_{" + str(int(angle)) + "}$"
        xcoord = (xdata[count])
        ycoord = (ydata[count]) + math.exp(math.log(ydata[count]) - dy)
        plt.text(xcoord, ycoord, sn)        