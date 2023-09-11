#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 16 07:02:38 2022

@author: bennett
"""
import numpy as np
import scipy.integrate as integrate
# import quadpy 
import matplotlib.pyplot as plt
from pathlib import Path
from ..solver_classes.functions import find_nodes
from ..solver_classes.functions import Pn

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
from ..solver_classes.opacity import sigma_integrator

from timeit import default_timer as timer
from .wavespeed_estimator import wavespeed_estimator
from .wave_loc_estimator import find_wave
from scipy.special import roots_legendre
import numpy.polynomial as poly
import scipy.special as sps
from functools import partial

"""
This file contains functions used by solver
"""



def parameter_function(major, N_spaces, Ms, count):
    if major == 'cells':
        M = Ms[count]
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

def solve(tfinal, N_space, N_ang, M, x0, t0, sigma_t, sigma_s, t_nodes, source_type, 
          uncollided, moving, move_type, thermal_couple, temp_function, rt, at, e_initial, choose_xs, specified_xs, 
          weights, sigma, particle_v, edge_v, cv0, estimate_wavespeed, find_wave_loc, thick, mxstp, wave_loc_array, 
          find_edges_tol, source_strength, move_factor, integrator, l, save_wave_loc, pad, leader_pad, xs_quad_order, 
          eval_times, eval_array, boundary_on, boundary_source_strength, boundary_source, sigma_func, Msigma,
          finite_domain, domain_width, fake_sedov_v0, test_dimensional_rhs, epsilon):

    # if weights == "gauss_lobatto":
    #     mus = quadpy.c1.gauss_lobatto(N_ang).points
    #     ws = quadpy.c1.gauss_lobatto(N_ang).weights
    # elif weights == "gauss_legendre":
    #     mus = quadpy.c1.gauss_legendre(N_ang).points
    #     ws = quadpy.c1.gauss_legendre(N_ang).weights
    # if N_ang == 2:
    mus, ws = quadrature(N_ang, weights, testing = True)
    #     print("mus =", mus)

    # xs_quad = quadpy.c1.gauss_legendre(2*M+1).points
    # ws_quad = quadpy.c1.gauss_legendre(2*M+1).weights

    xs_quad, ws_quad = quadrature(2*M+1, 'gauss_legendre')

    # t_quad = quadpy.c1.gauss_legendre(t_nodes).points
    t_quad, t_ws = quadrature(t_nodes, 'gauss_legendre')

    # t_ws = quadpy.c1.gauss_legendre(t_nodes).weights
    quad_thick_source, blank = quadrature(int(N_space/2+1), 'gauss_lobatto')
    quad_thick_edge, blank = quadrature(int(N_space/4+1), 'gauss_lobatto')
    # quad_thick_source = quadpy.c1.gauss_lobatto(int(N_space/2+1)).points
    # quad_thick_edge = quadpy.c1.gauss_lobatto(int(N_space/4+1)).points
    # quad_thick_source = ([quad_thick_source_inside, quad_thick_source_outside])
    
    initialize = build(N_ang, N_space, M, tfinal, x0, t0, mus, ws, xs_quad,
                       ws_quad, sigma_t, sigma_s, source_type, uncollided, moving, move_type, t_quad, t_ws,
                       thermal_couple, temp_function, e_initial, sigma, particle_v, edge_v, cv0, thick, 
                       wave_loc_array, source_strength, move_factor, l, save_wave_loc, pad, leader_pad, quad_thick_source,
                        quad_thick_edge, boundary_on, boundary_source_strength, boundary_source, sigma_func, Msigma,
                        finite_domain, domain_width, fake_sedov_v0, test_dimensional_rhs, epsilon)
                       
    initialize.make_IC()
    IC = initialize.IC

    if thermal_couple == 0:
        deg_freedom = N_ang*N_space*(M+1)
    elif thermal_couple == 1:
        deg_freedom = (N_ang+1)*N_space*(M+1)
    mesh = mesh_class(N_space, x0, tfinal, moving, move_type, source_type, edge_v, thick, move_factor,
                      wave_loc_array, pad, leader_pad, quad_thick_source, quad_thick_edge, finite_domain,
                      domain_width, fake_sedov_v0, boundary_on, t0) 
    matrices = G_L(initialize)
    num_flux = LU_surf(initialize)
    source = source_class(initialize)
    uncollided_sol = uncollided_solution(initialize)
    flux = scalar_flux(initialize)
    rhs = rhs_class(initialize)
    transfer = T_function(initialize)
    sigma_class = sigma_integrator(initialize)
    flux.load_AAA(sigma_class.AAA)
    
    def RHS(t, V):
        return rhs.call(t, V, mesh, matrices, num_flux, source, uncollided_sol, flux, transfer, sigma_class)
    
    start = timer()
    reshaped_IC = IC.reshape(deg_freedom)

    if estimate_wavespeed == False:
        tpnts = [tfinal]
    elif estimate_wavespeed == True:
        tpnts = np.linspace(0, tfinal, 10000)
    if eval_times == True:
        tpnts = eval_array
        print(tpnts, 'time points')

    
    sol = integrate.solve_ivp(RHS, [0.0,tfinal], reshaped_IC, method=integrator, t_eval = tpnts , rtol = rt, atol = at, max_step = mxstp)
    print(sol.y.shape,'sol y shape')
    print(eval_times, 'eval times')
    end = timer()
    print('solver finished')
    
    if save_wave_loc == True:
        print(save_wave_loc, 'save wave')
        wave_tpnts = rhs.times_list
        wave_xpnts = rhs.wave_loc_list
    else:
        wave_tpnts = np.array([0.0])
        wave_xpnts = np.array([0.0])

    # if estimate_wavespeed == True:
    #     wavespeed_array = wavespeed_estimator(sol, N_ang, N_space, ws, M, uncollided, mesh, 
    #                       uncollided_sol, thermal_couple, tfinal, x0)
    # elif estimate_wavespeed == False:
    wavespeed_array = np.zeros((1,1,1))

    if find_wave_loc == True:
        wave_loc_finder = find_wave(N_ang, N_space, ws, M, uncollided, mesh, uncollided_sol, 
        thermal_couple, tfinal, x0, sol.t, find_edges_tol, source_type, sigma_t)
        left_edges, right_edges, T_front_location = wave_loc_finder.find_wave(sol)
    elif find_wave_loc == False:
        left_edges =  np.zeros(1)
        right_edges = np.zeros(1)
        T_front_location = np.zeros(1)

    
    if thermal_couple == 0:
        extra_deg_freedom = 0
       
        sol_last = sol.y[:,-1].reshape((N_ang,N_space,M+1))
        if eval_times ==True:
            sol_array = sol.y.reshape((eval_array.size, N_ang,N_space,M+1)) 
    elif thermal_couple == 1:
        extra_deg_freedom = 1
        sol_last = sol.y[:,-1].reshape((N_ang+1,N_space,M+1))
        if eval_times == True:
            sol_array = sol.y.reshape((eval_array.size, N_ang+1,N_space,M+1)) 


    
    if sol.t.size > 1:
        timesteps = time_step_function(sol.t)
        print(np.max(timesteps), "max time step")
    
    mesh.move(tfinal)
    edges = mesh.edges
    
    if choose_xs == False:
        xs = find_nodes(edges, M)
        
    elif choose_xs == True:
        xs = specified_xs
    # print(xs, 'xs')
    if eval_times == False:
        output = make_output(tfinal, N_ang, ws, xs, sol_last, M, edges, uncollided)
        phi = output.make_phi(uncollided_sol)
        psi = output.psi_out # this is the collided psi
        exit_dist, exit_phi = output.get_exit_dist(uncollided_sol)
        xs_ret = xs
        if thermal_couple == 1:
            e = output.make_e()
        else:
            e = phi*0
    else:
        phi = np.zeros((eval_array.size, xs.size))
        psi = np.zeros((eval_array.size, N_ang, xs.size))
        exit_dist = np.zeros((eval_array.size, N_ang, 2))
        exit_phi = np.zeros((eval_array.size, 2))
        xs_ret = np.zeros((eval_array.size, xs.size))
        for it, tt in enumerate(eval_array):
            mesh.move(tt)
            edges = mesh.edges
            if choose_xs == False:
                xs = find_nodes(edges, M)
            elif choose_xs == True:
                xs = specified_xs
            output = make_output(tt, N_ang, ws, xs, sol.y[:,it].reshape((N_ang+extra_deg_freedom,N_space,M+1)), M, edges, uncollided)
            phi[it,:] = output.make_phi(uncollided_sol)
            psi[it, :, :] = output.psi_out # this is the collided psi
            exit_dist[it], exit_phi[it] = output.get_exit_dist(uncollided_sol)
            xs_ret[it] = xs
            if thermal_couple == 1:
                e = output.make_e()
            else:
                e = phi*0
    computation_time = end-start
    
    return xs_ret, phi, psi, exit_dist, exit_phi,  e, computation_time, sol_last, mus, ws, edges, wavespeed_array, tpnts, left_edges, right_edges, wave_tpnts, wave_xpnts, T_front_location, mus



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


def quadrature(n, name, testing = True):
    ws = np.zeros(n)
    xs = np.zeros(n)
    # roots, weights = roots_legendre(n-1)
    roots = np.zeros(n)
    if name == 'gauss_legendre':
        xs, ws = poly.legendre.leggauss(n)
    elif name == 'gauss_lobatto':
        if n > 1:
            # brackets = sps.legendre(n-1).weights[:, 0]
            xs_brackets, blanl = poly.legendre.leggauss(n-1)
            brackets = xs_brackets
        else:
            brackets = np.array([-1,1])
        for i in range(n-2):
            roots[i+1] = bisection(partial(eval_legendre_deriv, n-1),brackets[i], brackets[i+1])

    # mesh = np.linspace(-1, 1, 300)

    # plt.plot(mesh, eval_legendre_deriv(n-1, mesh))
    # plt.plot(roots, 0*roots, "o")


        xs = roots
        xs[0] = -1
        xs[-1] = 1
        for nn in range(1, n-1):
            inn = nn + 1
            ws[nn] = 2 / (n*(n-1)) / (Pn(n-1, np.array([roots[nn]]), -1.0, 1.0)[0])**2
            ws[0] = 2/ (n*(n-1))
            ws[-1] = 2/ (n*(n-1))
        # if testing == True:
        #     testxs = quadpy.c1.gauss_lobatto(n).points
        #     testws = quadpy.c1.gauss_lobatto(n).weights
        #     np.testing.assert_allclose(testxs, xs)
        #     np.testing.assert_allclose(testws, ws)


        # if testing == True:
        #     # testxs = quadpy.c1.gauss_legendre(n).points
        #     # testws = quadpy.c1.gauss_legendre(n).weights
        #     np.testing.assert_allclose(testxs, xs)
        #     np.testing.assert_allclose(testws, ws)
    return xs, ws





def bisection(f, a, b, tol=1e-14):
    assert np.sign(f(a)) != np.sign(f(b))
    while b-a > tol:
        m = a + (b-a)/2
        fm = f(m)
        if np.sign(f(a)) != np.sign(fm):
            b = m
        else:
            a = m
            
    return m

def eval_legendre_deriv(n, x):
    return (
        (x*sps.eval_legendre(n, x) - sps.eval_legendre(n-1, x))
        /
        ((x**2-1)/n))
