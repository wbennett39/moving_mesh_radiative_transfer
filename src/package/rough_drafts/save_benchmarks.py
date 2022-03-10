#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 11:56:33 2022

@author: bennett
"""
import h5py
import numpy as np

from create_benchmark import integrate_ganapol

###############################################################################

f = h5py.File("benchmarks_2_23.hdf5", "a")
npnts = 10000
x0 = 1/2

bench_maker = integrate_ganapol(npnts, x0)

bench_maker.plane_source(1.0)
p1 = bench_maker.sol_list
p1xs = bench_maker.x_list
# bench_maker.plane_source(5.0)
# p5 = bench_maker.sol_list
# p5xs = bench_maker.x_list
bench_maker.plane_source(10.0)
p10 = bench_maker.sol_list
p10xs = bench_maker.x_list

plane = f.create_group("plane_IC")
plane.create_dataset("t = 1", (2, npnts), dtype = "f", data = (p1xs, p1))
# plane.create_dataset("t = 5", (2, npnts), dtype = "f", data = (p5xs, p5))
plane.create_dataset("t = 10", (2, npnts), dtype = "f", data = (p10xs, p10))
print("plane finished")

npnts = 175
x0 = 1/2

bench_maker = integrate_ganapol(npnts, x0)

bench_maker.square_IC(1.0)
si1 = bench_maker.sol_list
si1xs = bench_maker.x_list
# bench_maker.square_IC(5.0)
# si5 = bench_maker.sol_list
# si5xs = bench_maker.x_list
bench_maker.square_IC(10.0)
si10 = bench_maker.sol_list
si10xs = bench_maker.x_list

print("square IC finished")

x0 = 5
bench_maker = integrate_ganapol(npnts, x0)

bench_maker.gaussian_IC(1.0)
tg1 = bench_maker.sol_list
tg1xs = bench_maker.x_list
# bench_maker.gaussian_IC(5.0)
# tg5 = bench_maker.sol_list
# tg5xs = bench_maker.x_list
bench_maker.gaussian_IC(10.0)
tg10 = bench_maker.sol_list
tg10xs = bench_maker.x_list

print("gauss finished")


square_IC = f.create_group("square_IC")
# square_source = f.create_group("square_source")
gaussian_IC = f.create_group("gaussian_IC")


square_IC.create_dataset("t = 1", (2, npnts), dtype = "f", data = (si1xs, si1))
# square_IC.create_dataset("t = 5", (2, npnts), dtype = "f", data = (si5xs, si5))
square_IC.create_dataset("t = 10", (2, npnts), dtype = "f", data = (si10xs, si10))

# square_source.create_dataset("t = 1", (2, npnts), dtype = "f")
# # square_source.create_dataset("t = 5", (2, npnts), dtype = "f")
# square_source.create_dataset("t = 10", (2, npnts), dtype = "f")

gaussian_IC.create_dataset("t = 1", (2, npnts), dtype = "f", data = (tg1xs, tg1))
# gaussian_IC.create_dataset("t = 5", (2, npnts), dtype = "f", data = (tg5xs, tg5))
gaussian_IC.create_dataset("t = 10", (2, npnts), dtype = "f", data = (tg10xs, tg10))

f.close()