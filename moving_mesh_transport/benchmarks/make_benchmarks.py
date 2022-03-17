#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 12 11:17:06 2022

@author: bennett
"""

import numpy as np
import math as math
import scipy.integrate as integrate 
# from numba import njit, cfunc, jit
import matplotlib.pyplot as plt
from timeit import default_timer as timer
import h5py
from .benchmark_functions import F, F1, F_gaussian_source, F1_integrand, uncollided_square_s2, pyF
from pathlib import Path

# @cfunc("complex128(float64, float64)")
# @jit

    
def opts0(*args, **kwargs):
       return {'limit':100000}
   
def opts1(*args, **kwargs):
       return {'limit':100000}
   
def opts2(*args, **kwargs):
       return {'limit':100000}

def do_ganapol(x, tfinal, x0):
    integral_1 = pyF(0.0, 0.0, tfinal, x)
    integral_2 = integrate.nquad(F1, [[0, math.pi]], args =  (0.0, 0.0, x, tfinal, 0), opts = [opts0])[0]
    return integral_1 + integral_2

def do_square_ic(x, tfinal, x0):
    integral_1 = integrate.nquad(F, [[-x0, x0]], args =  ([0.0, tfinal, x, 0]), opts = [opts0])[0]
    integral_2 = integrate.nquad(F1, [[0, math.pi], [-x0, x0]], args =  (0.0, x, tfinal,0), opts = [opts0, opts0, opts0])[0]
    return integral_1 + integral_2

def do_square_source(x, tfinal, x0):
    """ clean up and comment this"""
    # integral_1 = 0*integrate.nquad(F, [[-x0, x0], [0, tfinal]], args =  (tfinal, x), opts = [opts0, opts1, opts2])[0]
#    integral_1 = integrate.nquad(F1_integrand, [[0.0, tfinal]], args =  (tfinal, x, 0.5))[0]
    integral_1 = uncollided_square_s2(x, tfinal, x0, tfinal)
    integral_2 = integrate.nquad(F1, [[0, math.pi], [-x0, x0], [0, tfinal]], args =  (x, tfinal, 0), opts = [opts0, opts1, opts2])[0]
    # integral_1 = 0.0
    # integral_2 = 0
    return integral_1 + integral_2

def do_gaussian_ic(x, tfinal):
    integral_1 = integrate.nquad(F, [[-np.inf, np.inf]], args =  ([0.0, tfinal, x, 1]), opts = [opts0])[0]
    integral_2 = integrate.nquad(F1, [[0, math.pi], [-np.inf, np.inf]], args =  (0.0, x, tfinal, 1), opts = [opts0, opts1])[0]
    return integral_1 + integral_2
    

def do_gaussian_source(x, tfinal):
    sqrtpi = math.sqrt(math.pi)
    # integral_1 = integrate.nquad(F, [[-np.inf, np.inf], [0, tfinal]], args =  (tfinal, x, 1), opts = [opts0, opts1, opts2])[0]
    integral_1 = integrate.nquad(F_gaussian_source, [[0, tfinal]], args =  (tfinal, x), opts = [opts0])[0]

    integral_2 = integrate.nquad(F1, [[0, math.pi], [-np.inf, np.inf], [0, tfinal]], args =  (x, tfinal, 1), opts = [opts0, opts1, opts2])[0]
    return sqrtpi/8 *integral_1 + integral_2



    
def make_benchmark_file_structure():
    data_folder = Path("moving_mesh_transport/benchmarks")
    bench_file_path = data_folder / 'benchmarks.hdf5'
    source_name_list = ['plane_IC', 'square_IC', 'square_source', 'gaussian_IC', 'gaussian_source']
    
    f = h5py.File(bench_file_path, "a")
    
    for source_name in source_name_list:
        if f.__contains__(source_name):
            del f[source_name]
        f.create_group(source_name)
    
    f.close()

def write_to_file(xs, phi, tfinal, source_name, npnts):
    data_folder = Path("moving_mesh_transport/benchmarks")
    bench_file_path = data_folder / 'benchmarks.hdf5'
    
    with h5py.File(bench_file_path,'r+') as f:
        if f.__contains__(source_name + f'/t = {tfinal}'):
            del f[source_name + f'/t = {tfinal}'] 
        f.create_dataset(source_name + f'/t = {tfinal}', (2, npnts), dtype = "f", data=(xs, phi))
    f.close()
    

def make_benchmarks(tfinal, x0, npnts = [10000, 1000, 500, 500, 500]):
    print("t = ", tfinal)
    xs1 = np.linspace(0, tfinal, npnts[0])
    xs2 = np.linspace(0, tfinal + x0, npnts[1])
    xs3 = np.linspace(0, tfinal + x0, npnts[2])
    xs4 = np.linspace(0, tfinal + 5, npnts[3])
    xs5 = np.linspace(0, tfinal + 5, npnts[4])
    # xs3 = np.array([0.176,1.5])
    times = np.zeros(5)
    phi_pl = xs1*0
    phi_sq = xs2*0
    phi_sqs = xs3*0
    phi_gss = xs4*0
    phi_gss_s = xs5*0
    
    start = timer()
    for i in range(npnts[0]):
        phi_pl[i] = do_ganapol(xs1[i], tfinal, 0.0)
    times[0] = timer() - start
    start = timer()
    print("plane finished")
    for j in range(npnts[1]):
        phi_sq[j] = do_square_ic(xs2[j], tfinal, x0)
    times[1] = timer()  - start
    start = timer()
    print("square IC finished")
    for k in range(npnts[2]):
        phi_sqs[k] = do_square_source(xs3[k], tfinal, x0)
    times[2] = timer()  - start
    start = timer()
    print("square source finished")
    for h in range(npnts[3]):
        phi_gss[h] = do_gaussian_ic(xs4[h], tfinal)
    times[3] = timer() - start
    start = timer()
    print("Gauss IC finished")
    for l in range(npnts[4]):
        phi_gss_s[l] = do_gaussian_source(xs5[l], tfinal)
    times[4] = timer() - start
    print("Gauss source finished")
    
    
    # plt.plot(xs1, phi_pl, "-.",label = "plane")
    # plt.plot(xs2, phi_sq, "--", label = "square IC")
    # plt.plot(xs3, phi_sqs, ":", label = "square source")
    # plt.plot(xs4, phi_gss, "--*", label = "Gaussian IC")
    # plt.plot(xs5, phi_gss_s, "--x", label = "Gaussian source")
    # plt.xlabel("x")
    # plt.ylabel("scalar flux")
    # plt.xlim(0,tfinal + x0 + 1)
    # plt.legend()
    
    
    print("-   -   -   -   -   -   -   -   -")
    print("time elapsed")
    print(times)
    print("-   -   -   -   -   -   -   -   -")
    print("time per evaluation point")
    print(times/np.array(npnts))
    print("-   -   -   -   -   -   -   -   -")

    

    write_to_file(xs1, phi_pl, tfinal, 'plane_IC', npnts[0])
    write_to_file(xs2, phi_sq, tfinal, 'square_IC', npnts[1])
    write_to_file(xs3, phi_sqs, tfinal, 'square_source', npnts[2])
    write_to_file(xs4, phi_gss, tfinal, 'gaussian_IC', npnts[3])
    write_to_file(xs5, phi_gss_s, tfinal, 'gaussian_source', npnts[4])

def make_square_source(tfinal, x0, npnts = [100]):
    print("t = ", tfinal)
    xs1 = np.linspace(0, tfinal + x0, npnts[0])
    # xs3 = np.array([0.176,1.5])
    times = np.zeros(1)
    phi_sqs = xs1*0

    start = timer()
    for k in range(npnts[0]):
        phi_sqs[k] = do_square_source(xs1[k], tfinal, x0)
    times[0] = timer()  - start
    print("square source finished")
    
    plt.figure(-1)
    # plt.plot(xs1, phi_pl, "-.",label = "plane")
    # plt.plot(xs2, phi_sq, "--", label = "square IC")
    plt.plot(xs1, phi_sqs, "-", label = "square source")
    plt.show()
    # plt.plot(xs4, phi_gss, "--*", label = "Gaussian IC")
    # plt.plot(xs5, phi_gss_s, "--x", label = "Gaussian source")
    # plt.xlabel("x")
    # plt.ylabel("scalar flux")
    # plt.xlim(0,tfinal + x0 + 1)
    # plt.legend()
    
    
    print("-   -   -   -   -   -   -   -   -")
    print("time elapsed")
    print(times)
    print("-   -   -   -   -   -   -   -   -")
    print("time per evaluation point")
    print(times/np.array(npnts))
    print("-   -   -   -   -   -   -   -   -")
    write_to_file(xs1, phi_sqs, tfinal, 'square_source', npnts[0])
    

def make_all():
    x0 = 0.5
    make_benchmarks(1, x0)
    make_benchmarks(5, x0)
    make_benchmarks(10, x0)
