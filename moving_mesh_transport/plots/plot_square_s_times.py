#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 07:25:33 2022

@author: bennett
"""

from .make_plots import rms_plotter

def main():
    M = 6
    tfinal = 1
    
    plotter = rms_plotter(tfinal, M, "square_s")
    plotter.load_RMS_data(uncollided = True, moving = True)
    plotter.plot_RMS_vs_times(10, False)
    
    plotter = rms_plotter(tfinal, M, "square_s")
    plotter.load_RMS_data(uncollided = True, moving = False)
    plotter.plot_RMS_vs_times(10, False)
    
    plotter = rms_plotter(tfinal, M, "square_s")
    plotter.load_RMS_data(uncollided = False, moving = True)
    plotter.plot_RMS_vs_times(10, False)
    
    plotter = rms_plotter(tfinal, M, "square_s")
    plotter.load_RMS_data(uncollided = False, moving = False)
    plotter.plot_RMS_vs_times(10, False)