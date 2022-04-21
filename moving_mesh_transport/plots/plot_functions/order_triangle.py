#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 07:03:53 2022

@author: bennett
"""
import numpy as np
import matplotlib.pyplot as plt


def order_triangle(xstart, xend, slope, intercept, dy, dx):
    p1 = [xstart, (intercept-intercept/dx)*xstart**(-slope)]
    p2 = [xend, (intercept-intercept/dx)*(xend)**(-slope)]
    
    xs = np.array([p1[0], p2[0]])
    ys = np.array([p1[1], p2[1]])
    
    xs2 = np.array([p1[0],p1[0]])
    ys2 = ys
    
    xs3 = xs
    ys3 = np.array([p2[1], p2[1]])
    
    
    lwid = 2.5
    plt.loglog(xs,ys, 'k', linewidth = lwid)
    plt.loglog(xs2,ys2, 'k', linewidth = lwid)
    plt.loglog(xs3,ys3, 'k', linewidth = lwid)
    
    
    x1 = p1[0]
    x2 = p2[0]
    xmid = np.exp(np.log(x1) + np.log(x2/x1)/5) 
    
    plt.text(xmid, p2[1]-p2[1]/dy,f"ORDER {int(slope)}", fontsize = 8, fontweight = 'bold')
    
    plt.show()
    
    
# order_triangle(2,4,2,1e-3,1e-5,1e-5)
    
    
    