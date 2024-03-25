#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 10:05:00 2022

@author: bennett
"""
import sys
sys.path.append('/Users/bennett/Documents/Github/transport_benchmarks/')
sys.path.append('/Users/bennett/Documents/Github/transport_benchmarks/benchmarks')
from ..solver_classes.functions import *
from ..solver_classes.build_problem import build
from ..solver_classes.sources import source_class
from ..solver_classes.uncollided_solutions import uncollided_solution
from ..solver_functions.main_functions import quadrature
from ..solver_classes.functions import normTn
from benchmarks import integrate_greens as intg
import scipy.special as sp
import scipy.integrate as integrate
# from functions import *
import scipy.integrate as integrate
from numpy.testing import assert_allclose as npassert
from scipy.special import expi
import numpy as np
import quadpy
import math
from scipy.interpolate import interp1d
from numba import jit
from ..solver_classes.functions import quadrature

def test_normPns():
    # normalized
    for i in range(10):
         res = integrate.fixed_quad(lambda x: normPn(i,x,0.,2.)**2,0,2,n=i+1)[0]
         npassert(res, 1.0)
         # orthogonal 
         for j in range(10):
             if j  !=  i:
                 res = integrate.fixed_quad(lambda x: normPn(i,x,0.,2.) * normPn(j,x,0.,2.),0,2,n= i + j + 2)[0]
                 npassert(abs(res), 0.0, rtol = 1e-11, atol = 1e-11)
def test_expi():
    xs_neg = np.linspace(-25, -1)
    xs_pos = np.linspace(1, 25)
    for i in range(xs_neg.size):
        xn = xs_neg[i]
        xp = xs_pos[i]
        res = (numba_expi(xn) - expi(xn)) + (numba_expi(xp) - expi(xp)) 
        npassert(abs(res), 0.0, rtol = 1e-11, atol = 1e-11)
        
def test_node_finder_slab():
    edges = np.array([-1, 1])
    Ms = [1,2,3,4,5,6,7,8,9,10]
    geometry = {'slab': True, 'sphere': False}
    for M in Ms:
        nodes = find_nodes(edges, M, geometry)
        true_nodes = quadpy.c1.gauss_legendre(M)
        npassert(nodes, true_nodes.points, rtol = 1e-11, atol = 1e-11)

def test_node_finder_slab_2():
    aa = np.random.rand()*2-1
    bb = aa + np.random.rand()
    edges = np.array([aa, bb])
    Ms = [1,2,3,4,5,6,7,8,9,10]
    geometry = {'slab': True, 'sphere': False}
    for M in Ms:
        nodes = find_nodes(edges, M, geometry)
        basis_evaluated = normPn(M, nodes, aa, bb)
        for ix in range(nodes.size):
            assert(abs(basis_evaluated[ix]) <= 1e-10)

def test_node_finder_sphere():
    aa = np.random.rand()
    bb = aa + np.random.rand()
    edges = np.array([aa, bb])
    Ms = [1,2,3,4,5,6,7,8]
    geometry = {'slab': False, 'sphere': True}
    for M in Ms:
        nodes = find_nodes(edges, M, geometry)
        basis_evaluated = normTn(M, nodes, aa, bb)
        for ix in range(nodes.size):
            assert(abs(basis_evaluated[ix]) <= 1e-10)

    



def test_ortho_Tn():
    R = 3
    aalist = np.random.rand(4) * R
    bblist = aalist + np.random.rand(4) + 0.1
    for irun in range(4):
        aa = aalist[irun]
        bb = bblist[irun]
        assert(aa<bb)
        M = 4
        ortho_array = np.zeros((M+1, M+1))

        for i in range(M+1):
            for j in range(M+1): 
                integrand = lambda x: normTn(i, np.array([x]), aa, bb ) * weight_func_Tn(np.array([x]), aa, bb) * normTn(j, np.array([x]), aa, bb )
                ortho_array[i,j] = integrate.quad(integrand, aa, bb)[0]
        npassert(ortho_array, np.identity(M+1), atol = 1e-8)

def test_ortho_Tn_gauss_quad():
    R = 1
    aalist = np.random.rand(4) * R
    bblist = aalist + np.random.rand(4) + 0.1
    for irun in range(4):
        aa = aalist[irun]
        bb = bblist[irun]
        assert(aa<bb)
        M = 4
        ortho_array = np.zeros((M+1, M+1))
        ortho_array_gauss = np.zeros((M+1, M+1))
        xs_quad, ws_quad = quadrature(2*M+1, 'chebyshev')
        # xs_quad = quadpy.c1.chebyshev_gauss(250).points
        # ws_quad = quadpy.c1.gauss_chebyshev(250).weights
        aa = 0
        bb = 1
        for i in range(M+1):
            for j in range(M+1): 
                integrand = lambda x: normTn(i, x, aa, bb) * weight_func_Tn(x, aa, bb) * normTn(j, x, aa, bb)
                integrand2 = lambda x: normTn(i, x, aa, bb) * 2 * normTn(j, x, aa, bb)
                argument  = 0.5*(bb-aa)*xs_quad + (aa+bb)*0.5
                res = integrand2(argument)
                # integrandxi = lambda x: integrand(0.5*(bb-aa)*np.array([x]) + (aa+bb)*0.5)
                # print(res)
                # integral = integrate.quad(integrandxi, -1, 1)[0] * (bb-aa) * 0.5
                integral = 0.5*(bb-aa)*np.sum(ws_quad * res)

                # print(integral)
                ortho_array[i,j] = integral
        npassert(ortho_array, np.identity(M+1), atol = 1e-8)

def test_ortho_Pn():
    x1 = 7
    x2 = -7
    aalist = np.random.rand(4) * (x2-x1) - x1
    bblist = aalist + np.random.rand(4) + 0.1
    for irun in range(4):
        aa = aalist[irun]
        bb = bblist[irun]
        assert(aa<bb)
        M = 4
        ortho_array = np.zeros((M+1, M+1))
        for i in range(M+1):
            for j in range(M+1): 
                integrand = lambda x: normPn(i, x, aa, bb ) * normPn(j, x, aa, bb )
                ortho_array[i,j] = integrate.fixed_quad(integrand, aa, bb, n = 50 * M)[0]
        npassert(ortho_array, np.identity(M+1), atol = 1e-8)

def test_ortho_Pn_gauss_quad():
    R = 3
    aalist = np.random.rand(4) * R
    bblist = aalist + np.random.rand(4) + 0.1
    for irun in range(4):
        aa = aalist[irun]
        bb = bblist[irun]
        assert(aa<bb)
        M = 2
        ortho_array = np.zeros((M+1, M+1))
        ortho_array_gauss = np.zeros((M+1, M+1))
        # xs_quad, ws_quad = quadrature(2*M+1 + 300, 'gauss_legendre')
        xs_quad = quadpy.c1.gauss_legendre(150).points
        ws_quad = quadpy.c1.gauss_legendre(150).weights

        for i in range(M+1):
            for j in range(M+1): 
                integrand = lambda x: normPn(i, x, aa, bb )  * normPn(j, x, aa, bb)
                argument  = 0.5*(bb-aa)*xs_quad + (aa+bb)*0.5
                res = integrand(argument)
                # print(res)
                integral = 0.5*(bb-aa)*np.sum(ws_quad * res)
                # print(integral)
                ortho_array[i,j] = integral
        npassert(ortho_array, np.identity(M+1), atol = 1e-8)

def test_return_sol():
    test_func = lambda x: sp.jv(0, x) 
    M = 3
    cells = np.array([0, 0.25, 0.5, 0.75, 1.0])
    coeffs = np.zeros((cells.size, ))
    for ix in range(cells.size-1, M+1):
        aa = cells[ix]
        bb = cells[ix+1]
        for j in range(M+1):
            integrand = lambda x: Tn(i, np.array([x]), aa, bb ) * weight_func_Tn(np.array([x]), aa, bb) 
            coeffs[ix, j] = integrate.quad(integrand, aa, bb)[0]
    
def test_finite_diff():
    # testfunc = (1-x**2) * np.sin(x)
    pts = 5000
    xs = quadpy.c1.gauss_lobatto(pts).points
    sol_array = np.sin(xs) * (1-xs**2)
    deriv_array = xs*0
    true_derivative = (1 - xs**2)*np.cos(xs) - 2*xs*np.sin(xs)

    for ix in range(xs.size):
        deriv_array[ix] = finite_diff_uneven(xs, ix, sol_array , left = (ix == 0), right = (ix == xs.size-1))
    

    # plt.plot(xs, true_derivative)
    # plt.plot(xs, deriv_array-true_derivative, 'o', mfc = 'none')
    # plt.show()

    npassert(deriv_array, true_derivative, atol = 1e-6, rtol = 1e-6)


def test_Tn_edgevalues():
    cells = np.zeros(20)
    cells[0] = 0
    M = 6

    for ix in range(1, cells.size):
        cells[ix] = np.random.rand() + cells[ix-1]
    
    for ix in range(0, cells.size-1):
        a = cells[ix]
        b = cells[ix+1]
        h = math.sqrt(cells[ix+1] - cells[ix])

        edgeval = 1 / h / math.sqrt(math.pi)
        for j in range(M+1):
            if j == 0:
                B_left = edgeval
                B_right = edgeval
            else:
                B_right = math.sqrt(2) * edgeval
                if j%2 == 0: 
                    B_left = B_right
                else:
                    B_left = -B_right

            npassert(B_left, normTn(j, np.array([a]), a, b))
            npassert(B_right, normTn(j, np.array([b]), a, b))

           
           
def test_quadrature_chebyshev():
    nn = 50
    quad_test = np.polynomial.chebyshev.chebgauss(nn)
    my_quad_xs, my_quad_ws = quadrature(nn, 'chebyshev')
    npassert(quad_test[0], my_quad_xs)
    npassert(quad_test[1], my_quad_ws)

def test_quadrature_gss_lobatto():
    nn = 50
    quad_test_xs = quadpy.c1.gauss_lobatto(nn).points
    quad_test_ws = quadpy.c1.gauss_lobatto(nn).weights
    my_quad_xs, my_quad_ws = quadrature(nn, 'gauss_lobatto')
    npassert(quad_test_xs, my_quad_xs)
    npassert(quad_test_ws, my_quad_ws)



            




# def test_integrate_quad():
#     t = 0.5
    
#     edges = np.array([-t - 0.5, -0.5, 0.0, 0.5, t + 0.5])
    
#     N_ang = 128
#     N_space = 4
#     M = 5
#     tfinal = 5
#     x0 = 0.5
#     mus = quadpy.c1.gauss_lobatto(N_ang).points
#     ws = quadpy.c1.gauss_lobatto(N_ang).weights
#     xs_quad = quadpy.c1.gauss_legendre(M+2).points
#     ws_quad = quadpy.c1.gauss_legendre(M+2).weights
#     t_quad = quadpy.c1.gauss_legendre(100).points
#     t_ws = quadpy.c1.gauss_legendre(100).weights
#     sigma_t = np.ones(N_space)
#     sigma_s = np.ones(N_space)
#     move_type = np.array([1,0,0])
#     source_type = np.array([0,1,0,0,0,0])
#     uncollided = True
#     moving = True
    
#     initialize = build(N_ang, N_space, M, tfinal, x0, mus, ws, xs_quad, ws_quad, sigma_t, sigma_s, source_type,
#                  uncollided, moving, move_type, t_quad, t_ws)
    
#     source = source_class(initialize)
#     uncollided = uncollided_solution(initialize)
#     @jit
#     def func(x,t):
#         return np.cos(x)**2 + x**3 + 100*x
#     for ix in range(N_space):
#         for j in range(M+1):
#             xL = edges[ix]
#             xR = edges[ix+1]
#             source.integrate_quad(t, xL, xR, j, func)
#             normed_uncollided = lambda x: normPn(j, x, edges[ix], edges[ix+1]) * (np.cos(x)**2 + x**3 + 100*x)
#             res_2 = integrate.fixed_quad(normed_uncollided, edges[ix], edges[ix+1], n = M+2)[0]
#             npassert(source.S[j],res_2, rtol = 1e-11, atol = 1e-11)
    
    


# def test_Tn_expansion():
#     M = 4
#     a = 0.6
#     b = .9
#     xs = np.linspace(a, b)
#     res = np.zeros((M+1))

#     xs_quad, ws_quad = quadrature(2*2*M+1, 'chebyshev')
#     # xs_quad = quadpy.c1.gauss_legendre(2*M+1).points
#     # ws_quad = quadpy.c1.gauss_legendre(2*M+1).weight
#     # quadob = np.polynomial.chebyshev.chebgauss(10*2*M+1)
#     # xs_quad = quadob[0]
#     # ws_quad = quadob[1]
    

#     jn = 0
#     test_func = lambda x: sp.jv(jn, (a + b -2*x)/(a-b)) 

#     # Tn_basis = lambda x, n: eval_Tn_standard(n, 2*x**2-1) if n%2 == 0 else x * eval_Tn_standard(n, 2*x**2-1)
#         # if n == 0: 
#         #     return eval_Tn_standard(0,2*x**2-1)
#         # elif n % 2 == 0:
#         #     return eval_Tn_standard(n,2*x**2-1)
#         # elif n % 2 !=0: 
#         #     return x * eval_Tn_standard(n,2*x**2-1)
#     def basis(x, n):
#         return  eval_Tn_standard(n,2*x**2-1)

#     # Tn_basis = lambda x, n: eval_Tn_standard(n, x) 
#     # Tn_test1 = lambda x: eval_Tn_standard(2, x) * eval_Tn_standard(2, x) /np.sqrt(1-x**2)
#     # Tn_test2 = lambda x: eval_Tn_standard(3, x) * eval_Tn_standard(3, x) /np.sqrt(1-x**2)
#     # Tn_test3 = lambda x: eval_Tn_standard(2, x) * eval_Tn_standard(3, x) /np.sqrt(1-x**2)
    
#     # np.testing.assert_allclose(integrate.quad(Tn_test1,-1,1)[0], math.pi/2)
#     # np.testing.assert_allclose(integrate.quad(Tn_test2,-1,1)[0], math.pi/2)
#     # np.testing.assert_allclose(integrate.quad(Tn_test3,-1,1)[0], 0)

#     Tn_basis_test = lambda x: Tn_basis_shifted(x,3, a, b) * Tn_basis_shifted(x,3, a, b) * rho_shifted(x, a, b)
#     np.testing.assert_allclose(integrate.quad(Tn_basis_test, a, b)[0], math.pi/4)

#     for ic in range(M+1):
#         factor  = 4/math.pi
#         if ic == 0:
#             factor = 2/math.pi
#         else:
#             if ic % 2 != 0:
#                 factor = 4 / math.pi
                

        

#         # integrand = lambda x, n: test_func(x) * Tn_basis(x, n) / np.sqrt(1-x**2)
#         integrand = lambda x, n: test_func(x) * Tn_basis((a + b - 2*x) / (a-b), n)

#         integrand_vanilla = lambda x, n: test_func(x) * Tn_basis_shifted(x, n, a, b) * rho_shifted(x, a, b)     

#         res[ic] = factor * integrate.quad(integrand_vanilla, a, b, args = (ic))[0]

#         argument = (b-a)*0.5*xs_quad + (a+b)*0.5
#         # res[ic] =  factor * (b-a) * 0.5* np.sum(ws_quad * integrand(argument, ic))

#     # re-build solution
#     sol = 0*xs
#     # plt.figure(176)
#     # plt.plot(np.linspace(0, M, M+1), np.abs(res), 'o')
#     for j in range(M+1):
#         # sol += res[j] * eval_Tn(j, xs, a, b)
#         #  xnew = (b - xs) / (b-a)
#          xnew = (a + b - 2*xs) / (a-b)
#         #  sol += res[j] * basis(xs, j)
#          sol += res[j] * Tn_basis_shifted(xs, j,a ,b)

#     # plt.figure(177)
#     # plt.plot(xs, sol, label = 'interpolated')
#     # plt.plot(xs, test_func(xs), 'o', mfc = 'none',label = 'true')
#     # # plt.ylim(0,2)
#     # plt.xlabel('x')
#     # plt.title(f'j({jn},x)')
#     # plt.legend()
#     # plt.show()
#     np.testing.assert_allclose(sol, test_func(xs), atol =1e-7)

# def test_Tn():



#     xs = np.linspace(-1,1)
#     xshift = 2 * xs**2 -1
#     T0 = 1
#     T1 = xs*eval_Tn_standard(1, xshift)
#     T2 = eval_Tn_standard(2, xshift)
#     T3 = xs*eval_Tn_standard(3, xshift)
#     T4 = eval_Tn_standard(4, xshift)
#     np.testing.assert_allclose(T0, eval_Tn(0,xs, -1, 1))
#     np.testing.assert_allclose(T1, eval_Tn(1,xs, -1, 1))
#     # np.testing.assert_allclose(T2, eval_Tn(2,xs, -1, 1))



# def test_spherical_sol_interp():
# def get_coefficients(test_func, ic, a, b):
#     integrand_vanilla = lambda x, n: test_func(x,a,b) * Tn_basis_shifted(x, n, a, b) * rho_shifted(x, a, b)
#     res = factor_func(ic) * integrate.quad(integrand_vanilla, a, b, args = (ic))[0]
#     return res

# def factor_func(ic):
#     factor  = 4/math.pi
#     if ic == 0:
#         factor = 2/math.pi
#     else:
#         if ic % 2 != 0:
#             factor = 4 / math.pi
#     return factor

# def test_interpolate_point_source():
#     tf = 1.2
#     M = 4
#     xs = np.linspace(-tf, tf, 100)
# #    edges = np.array([0.0, 0.2, 0.5, 0.62, 0.84, 1.0032, 1.19])
#     edges = np.array([-1.19, 1.19])
#     sol_array = xs * 0
#     coefficient_array = np.zeros((edges.size-1, M+1))
#     sol_data = intg.point_source(tf, 500)
#     sol_data[1][0] = sol_data[1][1]
#     sol_data[1][-1] = sol_data[1][-2]
#     plt.close()
#     plt.close()
#     sol_interp = interp1d(sol_data[0], sol_data[1])

#     def test_func(xs, a, b):
#         return sol_interp(np.abs(xs))  #* Tn_basis_shifted(xs, 0, a, b)
#         a = 0
#         b = 1.2
#         # return sp.jv(3, (a + b -2*xs)/(a-b)) 
#         # return 1.0 + xs * 0
    
 


#     for icell in range(edges.size-1):
#         a = edges[icell]
#         b = edges[icell+1]
#         for j in range(M+1):
#             coefficient_array[icell, j] = get_coefficients(test_func, j, a, b)

#     for ix in range(xs.size):
#         icell = np.searchsorted(edges[:], xs[ix])
#         if (icell == 0):
#                 icell = 1
#         if (icell >= edges.size):
#                 icell = edges.size - 1
#         a = edges[icell-1]
#         b = edges[icell]
#         for j in range(M+1):
#             # if j == 0:
#             #     print(coefficient_array[icell, j])
#             # sol_array[ix] += Tn_basis_shifted(xs[ix], j, a, b) * coefficient_array[icell-1, j]
#             sol_array[ix] += Tn_basis_shifted(xs[ix], j, a, b) * coefficient_array[icell-1, j]

#     plt.figure(1)
#     for icell in range(edges.size-1):
#         a = edges[icell]
#         b = edges[icell+1]
#         leftind = np.argmin(np.abs(xs - a))
#         rightind = np.argmin(np.abs(xs - b))
#         plt.plot(xs[leftind:rightind], test_func(xs[leftind:rightind], a, b), c = 'b')
#     plt.plot(xs, sol_array, 'o', mfc = 'none', label = 'interpolated')
#     plt.legend()
#     plt.show()
#     # np.testing.assert_allclose(sol_array, test_func(xs))
    



# # def quadratic_TN()


# def Tn_basis(x, n):
#         return eval_Tn_standard(int(2 * n ), x)
    
# def Tn_basis_shifted(x, n, a, b):
#     return Tn_basis((a + b - 2*x) / (a-b), n) / math.sqrt(b-a)

# def rho_shifted(x, a, b):
#     xn = (a + b - 2*x) / (a-b)
#     return 1/np.sqrt(1-xn**2)

