#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 21:12:54 2022

@author: bennett
"""

import matplotlib.pyplot as plt
import h5py 

"""
- : moving mesh
-- : static
o : M = 2
^ : M = 4
s : M = 6
b : uncollided 
r : no-uncollided
"""

f = h5py.File('run_data_RMS.h5', 'r')

plt.figure(1)

n_ang = 512
plane_source_t_1_dataset = f[f"plane_IC/t=1/RMS/uncollided_moving_{n_ang}_angles_M_4"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "-^b")

plane_source_t_1_dataset = f[f"plane_IC/t=1/RMS/uncollided_static_{n_ang}_angles_M_4"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "--^b")

plane_source_t_1_dataset = f[f"plane_IC/t=1/RMS/no_uncollided_moving_{n_ang}_angles_M_4"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "-^r")

plane_source_t_1_dataset = f[f"plane_IC/t=1/RMS/no_uncollided_static_{n_ang}_angles_M_4"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "--^r")

plt.xlabel("cells")
plt.ylabel("root mean square error")
# plt.savefig("plane_IC_t=1.pdf")
plt.show()

plt.figure(2)

n_ang = 512
plane_source_t_1_dataset = f[f"square_IC/t=1/RMS/uncollided_moving_{n_ang}_angles_M_2"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "-ob")

plane_source_t_1_dataset = f[f"square_IC/t=1/RMS/uncollided_static_{n_ang}_angles_M_2"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "--ob")

plane_source_t_1_dataset = f[f"square_IC/t=1/RMS/no_uncollided_moving_{n_ang}_angles_M_2"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "-or")

plane_source_t_1_dataset = f[f"square_IC/t=1/RMS/no_uncollided_static_{n_ang}_angles_M_2"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "--or")

plt.xlabel("cells")
plt.ylabel("root mean square error")
# plt.savefig("square_IC_t=1_M4.pdf")

plt.show()


plt.figure(3)
plane_source_t_1_dataset = f[f"square_IC/t=1/RMS/uncollided_moving_{n_ang}_angles_M_4"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "-^b")

plane_source_t_1_dataset = f[f"square_IC/t=1/RMS/uncollided_static_{n_ang}_angles_M_4"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "--^b")

plane_source_t_1_dataset = f[f"square_IC/t=1/RMS/no_uncollided_moving_{n_ang}_angles_M_4"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "-^r")

plane_source_t_1_dataset = f[f"square_IC/t=1/RMS/no_uncollided_static_{n_ang}_angles_M_4"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "--^r")

plt.xlabel("cells")
plt.ylabel("root mean square error")
# plt.savefig("square_IC_t=1_M4.pdf")
plt.show()

plt.figure(4)
plane_source_t_1_dataset = f[f"square_IC/t=1/RMS/uncollided_moving_{n_ang}_angles_M_6"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "-sb")

plane_source_t_1_dataset = f[f"square_IC/t=1/RMS/uncollided_static_{n_ang}_angles_M_6"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "--sb")

plane_source_t_1_dataset = f[f"square_IC/t=1/RMS/no_uncollided_moving_{n_ang}_angles_M_6"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "-sr")

plane_source_t_1_dataset = f[f"square_IC/t=1/RMS/no_uncollided_static_{n_ang}_angles_M_6"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "--sr")

plt.xlabel("cells")
plt.ylabel("root mean square error")
# plt.savefig("square_IC_t=1_M6.pdf")
plt.show()



f.close()
