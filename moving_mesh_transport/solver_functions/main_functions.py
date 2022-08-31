#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 16 07:02:38 2022

@author: bennett
"""
import numpy as np
import scipy.integrate as integrate
import quadpy
import matplotlib.pyplot as plt
from pathlib import Path
from ..solver_classes.functions import find_nodes


from ..solver_classes.build_problem import build
from ..solver_classes.matrices import G_L
from ..solver_classes.numerical_flux import LU_surf
from ..solver_classes.sources import source_class
from ..solver_classes.uncollided_solutions import uncollided_solution
from ..solver_classes.phi_class import scalar_flux
from ..solver_classes.mesh import mesh_class
from ..solver_classes.rhs_class import rhs_class
from ..solver_classes.make_phi import make_output
from ..solver_classes.radiative_transfer import T_function
from timeit import default_timer as timer
from .wavespeed_estimator import wavespeed_estimator
from .wave_loc_estimator import find_wave

"""
This file contains functions used by solver
"""



def parameter_function(major, N_spaces, Ms, count):
    if major == 'cells':
        M = Ms[0]
        N_space = N_spaces[count]
    elif major == 'Ms':
        N_space = N_spaces[1]
        M = Ms[count]
    return N_space, M


def s2_source_type_selector(sigma, x0, thermal_couple, source_type, weights):
    """ 
    changes the name of the source type in order to select the correct 
    benchmark. For S2 benchmarks 
    """
    # thick source s8 
    if source_type[5] == 1:
        if sigma == 300:
            source_array_rad = [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0]
            source_array_mat = [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0]
        elif sigma == 0.5:
            source_array_rad = [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0]
            source_array_mat = [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0]
    elif source_type[2] == 1:
        if x0 == 400:
            source_array_rad = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0]
            source_array_mat = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]
        elif x0 == 0.5:
            source_array_rad = [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0]
            source_array_mat = [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0]
    return source_array_rad, source_array_mat
    

def time_step_function(t_array):
    N = len(t_array)
    res = np.zeros(N-1)
    for i in range(N-1):
        res[i] = t_array[i+1]-t_array[i]
    return res

def plot_p1_su_olson_mathematica():
    data_folder = Path("moving_mesh_transport/benchmarks")
    benchmark_mat_file_path = data_folder / "S2SuOlMat_t_1..txt"
    benchmark_rad_file_path = data_folder / "S2SuOlRadt_1..txt"
    
    su_olson_rad = np.loadtxt(benchmark_rad_file_path)
    su_olson_mat = np.loadtxt(benchmark_mat_file_path)
    plt.plot(su_olson_rad[:,0],su_olson_rad[:,1], "xk" )
    plt.plot(su_olson_mat[:,0],su_olson_mat[:,1], "xk" )
    
    return [su_olson_rad, su_olson_mat]

def solve(tfinal, N_space, N_ang, M, x0, t0, sigma_t, sigma_s, t_nodes, scattering_ratio, source_type, 
          uncollided, moving, move_type, thermal_couple, temp_function, rt, at, e_initial, choose_xs, specified_xs, 
          weights, sigma, particle_v, edge_v, cv0, estimate_wavespeed, find_wave_loc, thick, mxstp, wave_loc_array, find_edges_tol):

    if weights == "gauss_lobatto":
        mus = quadpy.c1.gauss_lobatto(N_ang).points
        ws = quadpy.c1.gauss_lobatto(N_ang).weights
    elif weights == "gauss_legendre":
        mus = quadpy.c1.gauss_legendre(N_ang).points
        ws = quadpy.c1.gauss_legendre(N_ang).weights
    if N_ang == 2:
        print("mus =", mus)
    xs_quad = quadpy.c1.gauss_legendre(2*M+1).points
    ws_quad = quadpy.c1.gauss_legendre(2*M+1).weights
    t_quad = quadpy.c1.gauss_legendre(t_nodes).points
    t_ws = quadpy.c1.gauss_legendre(t_nodes).weights
    initialize = build(N_ang, N_space, M, tfinal, x0, t0, scattering_ratio, mus, ws, xs_quad,
                       ws_quad, sigma_t, sigma_s, source_type, uncollided, moving, move_type, t_quad, t_ws,
                       thermal_couple, temp_function, e_initial, sigma, particle_v, edge_v, cv0, thick, wave_loc_array)
                       
    initialize.make_IC()
    IC = initialize.IC

    if thermal_couple == 0:
        deg_freedom = N_ang*N_space*(M+1)
    elif thermal_couple == 1:
        deg_freedom = (N_ang+1)*N_space*(M+1)
    mesh = mesh_class(N_space, x0, tfinal, moving, move_type, source_type, edge_v, thick, wave_loc_array) 
    matrices = G_L(initialize)
    num_flux = LU_surf(initialize)
    source = source_class(initialize)
    uncollided_sol = uncollided_solution(initialize)
    flux = scalar_flux(initialize)
    rhs = rhs_class(initialize)
    transfer = T_function(initialize)
    
    def RHS(t, V):
        return rhs.call(t, V, mesh, matrices, num_flux, source, uncollided_sol, flux, transfer)
    
    start = timer()
    reshaped_IC = IC.reshape(deg_freedom)

    if estimate_wavespeed == False:
        tpnts = [tfinal]
    elif estimate_wavespeed == True:
        tpnts = np.linspace(0, tfinal, 25)
    
    sol = integrate.solve_ivp(RHS, [0.0,tfinal], reshaped_IC, method='DOP853', t_eval = tpnts , rtol = rt, atol = at, max_step = mxstp)
    end = timer()

    if estimate_wavespeed == True:
        wavespeed_array = wavespeed_estimator(sol, N_ang, N_space, ws, M, uncollided, mesh, 
                          uncollided_sol, thermal_couple, tfinal, x0)
    elif estimate_wavespeed == False:
        wavespeed_array = np.array([[0],[0], [0]])

    if find_wave_loc == True:
        wave_loc_finder = find_wave(N_ang, N_space, ws, M, uncollided, mesh, uncollided_sol, 
        thermal_couple, tfinal, x0, sol.t, find_edges_tol)
        left_edges, right_edges = wave_loc_finder.find_wave(sol)
    elif find_wave_loc == False:
        left_edges =  np.array([0])
        right_edges = np.array([0])


    if thermal_couple == 0:
        sol_last = sol.y[:,-1].reshape((N_ang,N_space,M+1))
    elif thermal_couple == 1:
        sol_last = sol.y[:,-1].reshape((N_ang+1,N_space,M+1))

    
    if sol.t.size > 1:
        timesteps = time_step_function(sol.t)
        print(np.max(timesteps), "max time step")
    
    mesh.move(tfinal)
    edges = mesh.edges
    
    if choose_xs == False:
        xs = find_nodes(edges, M)
        
    elif choose_xs == True:
        xs = specified_xs
        
    output = make_output(tfinal, N_ang, ws, xs, sol_last, M, edges, uncollided)
    phi = output.make_phi(uncollided_sol)
    if thermal_couple == 1:
        e = output.make_e()
    else:
        e = phi*0
    
    computation_time = end-start
    
    return xs, phi, e, computation_time, sol_last, ws, edges, wavespeed_array, tpnts, left_edges, right_edges



def problem_identifier():
    name_array = []

def plot_edges(edges, fign):
    plt.figure(fign)
    for ed in range(edges.size):
        plt.scatter(edges[ed], 0.0, s = 128, c = 'k', marker = "|")

def x0_function(x0, source_type, count):
        if source_type[3] or source_type[4] == 1:
            x0_new = x0[count]
        else:
            x0_new = x0[0]
        return x0_new