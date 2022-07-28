#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 17 10:06:21 2022

@author: bennett
"""

import numpy as np
import scipy.interpolate as interpolate

from ..solver_functions.main_functions import plot_p1_su_olson_mathematica

def test_P1_against_mathematica(tfinal, xs, phi_sol, rad_or_mat):
    mathematica_rad = plot_p1_su_olson_mathematica()[0]
    mathematica_mat = plot_p1_su_olson_mathematica()[1]
    mathematica_xs = mathematica_rad[:,0]
    
    phi_interpolated = interpolate.interp1d(xs, phi_sol)
    # e_interpolated = interpolate(xs, e_sol)

    if rad_or_mat == "rad":
        diff_phi = phi_interpolated(mathematica_xs) - mathematica_rad[:,1]
    elif rad_or_mat == "mat":
        diff_phi = phi_interpolated(mathematica_xs) - mathematica_mat[:,1]    
    
    # np.testing.assert_allclose(diff_phi, np.zeros(len(mathematica_xs)), rtol = 1e-7, atol = 1e-7)
    
    
