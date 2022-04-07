#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 07:07:19 2022

@author: bennett
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d as interp

from ..load_bench import load_bench
from pathlib import Path

data_folder = Path("moving_mesh_transport/benchmarks")

bench_object = load_bench(np.array([0,0,0,0,0,0,1]), 1, 0.5)

rhos = np.linspace(0, 2.9)

phi_min = np.loadtxt(data_folder / 'cut_rho200_p99.dat', unpack = True)
rhos_min = np.loadtxt(data_folder / 'cut_x200.dat', unpack = True)

my_solution = bench_object(rhos)[0]
plt.plot(rhos_min, phi_min, "-k")

minwoos = interp(rhos_min, phi_min, kind = "cubic")

RMS = np.sqrt(np.mean(((my_solution) - minwoos(rhos))**2))
print("RMS = ", RMS)

plt.plot(rhos_min, phi_min, "-k")
plt.plot(rhos, my_solution, "-s", mfc = 'none')

plt.show(block=False)