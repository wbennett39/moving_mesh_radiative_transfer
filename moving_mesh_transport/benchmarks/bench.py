#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 20:25:02 2022

@author: bennett
"""

from .benchmarks import make_benchmark

def plane_IC(t, npnts):
    bench_class = make_benchmark('plane_IC', 1e-16, 1e-16)
    bench_class.integrate(t, npnts)
    bench_class.save()
    bench_class.plot()
    
def square_IC(t, npnts):
    bench_class = make_benchmark('square_IC', 0.5, 1e-16)
    bench_class.integrate(t, npnts)
    bench_class.save()
    bench_class.plot()
    
def square_source(t, npnts):
    bench_class = make_benchmark('square_source', 0.5, 5.0)
    bench_class.integrate(t, npnts)
    bench_class.save()
    bench_class.plot()
    
def gaussian_IC(t, npnts):
    bench_class = make_benchmark('gaussian_IC', 4.0, 1e-16)
    bench_class.integrate(t, npnts)
    bench_class.save()
    bench_class.plot()
    
def gaussian_source(t, npnts):
    bench_class = make_benchmark('gaussian_source', 4.0, 5.0)
    bench_class.integrate(t, npnts)
    bench_class.save()
    bench_class.plot()
    
def do_all(npnts = [5000, 5000, 5000, 500, 500]):
    plane_IC(1, npnts[0])
    plane_IC(5, npnts[0])
    plane_IC(10, npnts[0])
    print("plane finished")
    square_IC(1, npnts[1])
    square_IC(5, npnts[1])
    square_IC(10, npnts[1])
    print("square IC finished")
    gaussian_IC(1, npnts[2])
    gaussian_IC(5, npnts[2])
    gaussian_IC(10, npnts[2])
    print("gaussian IC finished")
    square_source(1, npnts[3])
    square_source(5, npnts[3])
    square_source(10, npnts[3])
    print("square source finished")
    gaussian_source(1, npnts[4])
    gaussian_source(5, npnts[4])
    gaussian_source(10, npnts[4])
    print("gaussian source finished")
    