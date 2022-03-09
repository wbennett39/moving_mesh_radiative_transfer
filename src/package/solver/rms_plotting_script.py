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

# n_ang = 512
plane_source_t_1_dataset = f[f"plane_IC/t=1/RMS/uncollided_moving_M_4"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "-^b")

plane_source_t_1_dataset = f[f"plane_IC/t=1/RMS/uncollided_static_M_4"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "--^b")

plane_source_t_1_dataset = f[f"plane_IC/t=1/RMS/no_uncollided_moving_M_4"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "-^r", mfc = "none")

plane_source_t_1_dataset = f[f"plane_IC/t=1/RMS/no_uncollided_static_M_4"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "--^r", mfc = "none")

plt.xlabel("cells")
plt.ylabel("root mean square error")
# plt.savefig("plane_IC_t=1.pdf")
plt.show()

plt.figure(2)

n_ang = 512
plane_source_t_1_dataset = f[f"square_IC/t=1/RMS/uncollided_moving_M_2"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "-ob")

plane_source_t_1_dataset = f[f"square_IC/t=1/RMS/uncollided_static_M_2"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "--ob")

plane_source_t_1_dataset = f[f"square_IC/t=1/RMS/no_uncollided_moving_M_2"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "-or", mfc = "none")

plane_source_t_1_dataset = f[f"square_IC/t=1/RMS/no_uncollided_static_M_2"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "--or", mfc = "none")

plt.xlabel("cells")
plt.ylabel("root mean square error")
# plt.savefig("square_IC_t=1_M2.pdf")

plt.show()


plt.figure(3)
plane_source_t_1_dataset = f[f"square_IC/t=1/RMS/uncollided_moving_M_4"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "-^b", mfc = "none")

plane_source_t_1_dataset = f[f"square_IC/t=1/RMS/uncollided_static_M_4"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "--^b", mfc = "none")

plane_source_t_1_dataset = f[f"square_IC/t=1/RMS/no_uncollided_moving_M_4"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "-^r", mfc = "none")

plane_source_t_1_dataset = f[f"square_IC/t=1/RMS/no_uncollided_static_M_4"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "--^r", mfc = "none")

plt.xlabel("cells")
plt.ylabel("root mean square error")
# plt.savefig("square_IC_t=1_M4.pdf")
plt.show()

plt.figure(4)
plane_source_t_1_dataset = f[f"square_IC/t=1/RMS/uncollided_moving_M_6"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "-sb")

plane_source_t_1_dataset = f[f"square_IC/t=1/RMS/uncollided_static_M_6"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "--sb")

plane_source_t_1_dataset = f[f"square_IC/t=1/RMS/no_uncollided_moving_M_6"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "-sr", mfc = "none")

plane_source_t_1_dataset = f[f"square_IC/t=1/RMS/no_uncollided_static_M_6"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "--sr", mfc = "none")

plt.xlabel("cells")
plt.ylabel("root mean square error")
# plt.savefig("square_IC_t=1_M6.pdf")
plt.show()

plt.figure(5)
t_1_dataset = f[f"gaussian_IC/t=1/RMS/uncollided_moving_M_2"]
spaces = t_1_dataset[0]
RMS = t_1_dataset[1]
plt.loglog(spaces, RMS, "-ob")

# t_1_dataset = f[f"gaussian_IC/t=1/RMS/uncollided_static_M_2"]
# spaces = t_1_dataset[0]
# RMS = t_1_dataset[1]
# plt.loglog(spaces, RMS, "--ob")

# t_1_dataset = f[f"gaussian_IC/t=1/RMS/no_uncollided_moving_M_2"]
# spaces = t_1_dataset[0]
# RMS = t_1_dataset[1]
# plt.loglog(spaces, RMS, "-or")

# t_1_dataset = f[f"gaussian_IC/t=1/RMS/no_uncollided_static_M_2"]
# spaces = t_1_dataset[0]
# RMS = t_1_dataset[1]
# plt.loglog(spaces, RMS, "--or")




t_1_dataset = f[f"gaussian_IC/t=1/RMS/uncollided_moving_M_4"]
spaces = t_1_dataset[0]
RMS = t_1_dataset[1]
plt.loglog(spaces, RMS, "-^b")

# t_1_dataset = f[f"gaussian_IC/t=1/RMS/uncollided_static_M_4"]
# spaces = t_1_dataset[0]
# RMS = t_1_dataset[1]
# plt.loglog(spaces, RMS, "--^b")

# t_1_dataset = f[f"gaussian_IC/t=1/RMS/no_uncollided_moving_M_4"]
# spaces = t_1_dataset[0]
# RMS = t_1_dataset[1]
# plt.loglog(spaces, RMS, "-^r")

# t_1_dataset = f[f"gaussian_IC/t=1/RMS/no_uncollided_static_M_4"]
# spaces = t_1_dataset[0]
# RMS = t_1_dataset[1]
# plt.loglog(spaces, RMS, "--^r")




t_1_dataset = f[f"gaussian_IC/t=1/RMS/uncollided_moving_M_6"]
spaces = t_1_dataset[0]
RMS = t_1_dataset[1]
plt.loglog(spaces, RMS, "-sb")

# t_1_dataset = f[f"gaussian_IC/t=1/RMS/uncollided_static_M_6"]
# spaces = t_1_dataset[0]
# RMS = t_1_dataset[1]
# plt.loglog(spaces, RMS, "--sb")

# t_1_dataset = f[f"gaussian_IC/t=1/RMS/no_uncollided_moving_M_6"]
# spaces = t_1_dataset[0]
# RMS = t_1_dataset[1]
# plt.loglog(spaces, RMS, "-sr")

# t_1_dataset = f[f"gaussian_IC/t=1/RMS/no_uncollided_static_M_6"]
# spaces = t_1_dataset[0]
# RMS = t_1_dataset[1]
# plt.loglog(spaces, RMS, "--sr")

# plt.xlabel("cells")
# plt.ylabel("root mean square error")

plt.savefig("gaussian_IC_t=1_M2_4_6_uncollided.pdf")
plt.show()


plt.figure(6)

""" plots RMS for plane IC all four cases M = 4
"""
plane_source_t_1_dataset = f[f"plane_IC/t=10/RMS/uncollided_moving_M_4"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "-^b")

plane_source_t_1_dataset = f[f"plane_IC/t=10/RMS/uncollided_static_M_4"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "--^b")

plane_source_t_1_dataset = f[f"plane_IC/t=10/RMS/no_uncollided_moving_M_4"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "-^r", mfc = "none")

plane_source_t_1_dataset = f[f"plane_IC/t=10/RMS/no_uncollided_static_M_4"]
spaces = plane_source_t_1_dataset[0]
RMS = plane_source_t_1_dataset[1]
plt.loglog(spaces, RMS, "--^r", mfc = "none")

plt.xlabel("cells")
plt.ylabel("root mean square error")
# plt.savefig("plane_IC_t=1.pdf")
plt.show()



plt.figure(7)
""" Plots RMS for MMS M=2, 4, 6

"""

t_1_dataset = f[f"MMS/t=1/RMS/no_uncollided_moving_M_2"]
spaces = t_1_dataset[0]
RMS = t_1_dataset[1]
plt.loglog(spaces, RMS, "-or", mfc = "none")

t_1_dataset = f[f"MMS/t=1/RMS/no_uncollided_moving_M_4"]
spaces = t_1_dataset[0]
RMS = t_1_dataset[1]
plt.loglog(spaces, RMS, "-^r", mfc = "none")

t_1_dataset = f[f"MMS/t=1/RMS/no_uncollided_moving_M_6"]
spaces = t_1_dataset[0]
RMS = t_1_dataset[1]
plt.loglog(spaces, RMS, "-sr", mfc = "none")




plt.xlabel("cells")
plt.ylabel("root mean square error")
plt.savefig("MMS_t=1.pdf")
plt.show()

f.close()
