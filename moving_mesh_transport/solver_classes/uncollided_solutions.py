#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 09:03:23 2022

@author: bennett
"""
import numpy as np
import math

from .build_problem import build
from .functions import normPn, numba_expi, uncollided_su_olson_s2, uncollided_square_s2
from .functions import uncollided_s2_gaussian, uncollided_s2_gaussian_thick
from.functions import heaviside_vector, heaviside_scalar, normTn
# from scipy.special import expi as expi2
import numba as nb
from numba import float64, int64, deferred_type
from numba.experimental import jitclass
###############################################################################
"""
I wonder if it would be faster for this function and other source functions to
make temp a self. variable and call it instead of returning it 
"""
###############################################################################

build_type = deferred_type()
build_type.define(build.class_type.instance_type)
params_default = nb.typed.Dict.empty(key_type=nb.typeof('par_1'),value_type=nb.typeof(1))

data = [("S", float64[:]),
        ("source_type", int64[:]),
        ("uncollided", int64),
        ("moving", int64),
        ("M", int64),
        ("x0", float64),
        ("t", float64),
        ("xs_quad", float64[:]),
        ("ws_quad", float64[:]),
        ("tfinal", float64),
        ("t0", float64),
        ("sqs_interval_1", float64),
        ("sqs_interval_2", float64),
        ("sqs_interval_3", float64),
        ("sqs_interval_4", float64),
        ("t_quad", float64[:]),
        ("t_ws", float64[:]),
        ("N_ang", int64),
        ("sigma", float64),
        ('source_strength', float64),
        ('uncollided_solution_return', float64[:]),
        ('v0', float64),
        ('ws', float64[:]),
        ('mus', float64[:]),
        ('t0_source', float64),
        ('geometry', nb.typeof(params_default)),
        ]
###############################################################################
@jitclass(data)
class uncollided_solution(object):
    def __init__(self, build):
        self.source_type =np.array(list(build.source_type), dtype = np.int64) 
        self.uncollided = build.uncollided
        self.x0 = build.x0
        self.M = build.M
        self.xs_quad = build.xs_quad
        self.ws_quad = build.ws_quad
        self.moving = build.moving
        self.tfinal = build.tfinal
        self.t0 = build.t0
        self.sqs_interval_1 = 0.0
        self.sqs_interval_2 = 0.0
        self.sqs_interval_3 = 0.0
        self.sqs_interval_4 = 0.0
        self.t_quad = build.t_quad
        self.t_ws = build.t_ws
        self.N_ang = build.N_ang
        self.sigma = build.sigma
        self.source_strength = build.source_strength
        self.uncollided_solution_return = np.zeros(1)
        self.v0 = build.fake_sedov_v0
        self.ws = build.ws
        self.mus = build.mus
        self.t0_source = build.t0
        self.geometry = build.geometry
###############################################################################
        
    def integrate_quad_gaussian_source(self, t, x, a, b, func):
        """ integrates the gaussian source over tau from 0 to min(t,t0)
        """
        argument = (b-a)/2 * self.t_quad + (a+b)/2
        return(b-a)/2 * np.sum(self.t_ws * func(argument, t, x)) 
    
    def integrate_quad_sedov_source(self, t, x, a, b, func):

        # argument = (b-a)/2 * self.t_quad + (a+b)/2
        argument = self.mus
        funcval = argument * 0
        for imu, xmu in enumerate(argument):
            funcval[imu] = func(xmu, t, x)
        # return (b-a)/2 * np.sum(self.t_ws * funcval) 
        return 2*np.sum(np.multiply(self.ws, funcval))
        
    def gaussian_source_integrand(self, tau, t, x):
        abx = abs(x)
        temp = tau*0
        for i in range(tau.size):
            tp = t - tau[i]
            
            if tp != 0:
                erf1 = math.erf((tp - abx)/self.sigma) 
                erf2 = math.erf((tp + abx)/self.sigma)
                temp[i] = math.exp(-tp) * (erf1 + erf2) / tp / 4.0
            else:
                temp[i] = 0.0
        return temp

    def fake_sedov_integrand(self, mu, t, x):
        c1 = 1.0
        v0 = self.v0
        x0 = -self.x0
        t0 =  (x0-x)/mu + t # time the particle is emitted
        x02 = 0.0
        sqrt_pi = math.sqrt(math.pi)
        kappa = 2
        rho0 = 0.2
        b2 = ((v0*x0) - t0*v0*mu)/(v0 + mu)
        b1 = max(x, b2)
        b4 = x0
        b3 =  min(x, b2)
        def t1(s):
            return (sqrt_pi*kappa*mu*math.erf((v0*(s - x0) + (c1 + s + t0*v0)*mu)/(kappa*mu)))/(2.*(v0 + mu))
        def t2(s):
            return rho0 * s

        mfp = t1(b1) - t1(b2) + t2(b3) - t2(b4) 
        return np.exp(-mfp / mu) * heaviside_scalar(mu - abs((x - x0)/ (t+1e-10))) * heaviside_scalar(abs(x0-x) - mu*(t-self.t0_source))
    
        
    def gaussian_source_uncollided_solution(self, xs, t):
        temp = xs*0
        sqrtpi = math.sqrt(math.pi)
        t_ceiling = min(t,self.t0)
        for ix in range(xs.size):
            x = xs[ix]
            result = self.integrate_quad_gaussian_source(t, x, 0.0, t_ceiling, self.gaussian_source_integrand)

            temp[ix] = result
        return temp * sqrtpi * self.sigma   
    
    
    
    
    
    # def sq_s_f1(self, t):
    #     # print(self.sqs_interval_2)
    #     arg1 = self.sqs_interval_2-t
    #     arg2 = 0.0 - t
    #     # print("################")
    #     # print(expi(arg2))
    #     print(numba_expi(-1))
    #     print(expi(0.0-1.0))
    #     # print("##################")
    #     return -self.x0 * expi(arg1) - (-self.x0 * expi(arg2))
        
    # def sq_s_f2(self, t, x):
    #     # print(self.sqs_interval_2)
    #     # print(self.sqs_interval_3)
    #     if self.sqs_interval_2 != t:
    #         eval_left = 0.5*((-self.x0 + abs(x)) * expi(self.sqs_interval_2-t) + math.exp(self.sqs_interval_2 - t))
    #     else:
    #         eval_left = 0.5
            
    #     if self.sqs_interval_3 != t:
    #         eval_right = 0.5*((-self.x0 + abs(x)) * expi(self.sqs_interval_3-t) + math.exp(self.sqs_interval_3 - t))
    #     else:
    #         eval_right = 0.5
    #     return eval_right-eval_left

    # def sq_s_f3(self, t):
    #     # print(self.sqs_interval_3)
    #     # print(self.sqs_interval_4)
    #     return math.exp(self.sqs_interval_4-t) - math.exp(self.sqs_interval_3-t)
    
    def square_source_uncollided_solution(self, xs, t):
        # i dont think thid is used
        temp = xs*0
        for ix in range(xs.size):
            x = xs[ix]
            temp[ix] = uncollided_square_s2(x, t, self.x0, self.t0)
            
            # t1 = 0.0
            # t2 = 0.0
            # t3 = 0.0
            # end = min(self.t0, t - abs(x) + self.x0)
            # self.sqs_interval_2 = min(end, t - self.x0 - abs(x))
            # if self.sqs_interval_2 < 0.0:
            #     self.sqs_interval_2 = 0.0
            # self.sqs_interval_3 = min(end, t - self.x0 + abs(x))
            # if self.sqs_interval_3 < 0.0:
            #     self.sqs_interval_3 = 0.0
            # self.sqs_interval_4 = end
            # if self.sqs_interval_4 < 0:
            #     self.sqs_interval_4 = 0.0
            # t1 = self.sq_s_f1(t)
            # t2 = self.sq_s_f2(t, x)
            # t3 = self.sq_s_f3(t)
            # temp[ix] = t1 + t2 + t3
            
        return temp
        
    def square_IC_uncollided_solution(self, xs, t):
        temp = xs*0
        for ix in range(xs.size):
            xx = xs[ix]
            # abxx = abs(xx)
            if (t <= self.x0) and (xx >= -self.x0 + t) and (xx <= self.x0 - t):
                temp[ix] = math.exp(-t)
            elif t > self.x0  and (-t + self.x0 <=  xx) and (t - self.x0 >= xx):
                temp[ix] = math.exp(-t) * self.x0 / (t + 1e-12)
            elif (xx < t + self.x0) and (xx > -t - self.x0):
                if (self.x0 - xx >= t) and (self.x0 + xx <= t):
                    temp[ix] = math.exp(-t)*(t + xx + self.x0)/(2.0 * t + 1e-12)
                elif (self.x0 - xx <= t) and (self.x0 + xx >= t):
                    temp[ix] = math.exp(-t)*(t - xx + self.x0)/(2.0 * t + 1e-12)
        return temp * self.source_strength
    
    def shell_IC_uncollided_solution(self, xs, t):
        temp = xs*0 
        sigma_abs = 1
        N01  = 4 * math.pi * self.x0**3 / 3
        N0 = 1
        n0 = N0 / (4. * math.pi / 3. * (self.x0 ** 3)) / (4. * math.pi)
        for ix, r in enumerate(xs):
            tt = t + 1e-12
            mu_crit = min(1., max(-1.,0.5*(tt/r+r/tt-self.x0**2/(r*tt))))
            r2 = r ** 2 + t ** 2 - 2 * mu_crit * r * t
            # if np.sqrt(r2) < self.x0: 
            temp[ix] = n0 * 2 * math.pi *  (1. - mu_crit ) * np.exp(-t * sigma_abs) 
        return temp


    def gaussian_IC_uncollided_solution(self, xs, t):
        temp = xs*0
        sqrtpi = math.sqrt(math.pi)
        for ix in range(xs.size):
            xx = xs[ix]
            temp[ix] = self.sigma * math.exp(-t) * sqrtpi * (math.erf((t-xx)/self.sigma) + math.erf((t+xx)/self.sigma))/(4.0 * t + 1e-12) * self.source_strength
        return temp
        
    def plane_IC_uncollided_solution(self, xs, t):
        temp = xs*0
        for ix in range(xs.size):
            if (-t <= xs[ix] <= t):
 
                temp[ix] = math.exp(-t)/(2*t+1e-12) 
                # elif self.geometry['sphere'] == True:
                #     eta = xs[ix] / t
                #     temp[ix] = math.exp(-t)/(2* math.pi * t**2 * math.sqrt(1-eta**2) + 1e-12 ) 
        return temp
    
    def plane_IC_uncollided_solution_integrated(self, t, xL, xR):
        #since the solution is piecewise constant, it is easy to integrate over
        # a cell
        return  (math.exp(-t)/(2*t+self.x0)*math.sqrt(xR-xL))
    
    def su_olson_s2_uncollided_solution(self, xs, t):
        temp = xs*0
        for ix in range(xs.size):
            temp[ix] = uncollided_su_olson_s2(xs[ix], t, self.x0, self.t0)
        return temp
    
    def gaussian_s2_thick_integrand(self, taus, t, x):
        temp = taus*0
        for it in range(taus.size):
            tau = taus[it]
            temp[it] = (math.exp(-t + tau)*(math.exp(-(math.sqrt(3)*t + 3*x - math.sqrt(3)*tau)**2/(9.*self.sigma**2)) + math.exp(-(-(math.sqrt(3)*t) + 3*x + math.sqrt(3)*tau)**2/(9.*self.sigma**2))))/2.
        return temp
    
    def gaussian_s2(self, xs, t):
        temp = xs*0
        for ix in range(xs.size):
            if self.sigma <= 5:
                temp[ix] = uncollided_s2_gaussian(xs[ix], t, self.sigma, self.t0)
            else:
                temp[ix] = self.integrate_quad_gaussian_source(t, xs[ix], 0.0, min(t,self.t0), self.gaussian_s2_thick_integrand)
        return temp
    
    def fake_sedov(self, xs, t):
        temp = xs*0
        for ix, xx in enumerate(xs):
            temp[ix] = self.integrate_quad_sedov_source(t, xx, -1, 1, self.fake_sedov_integrand)
        return temp
    
    def point_source(self, rhos, t):   
        temp = rhos*0
        for ix in range(rhos.size):
            rho = rhos[ix]
            if abs(rho-t) <= 1e-10:
                temp[ix] = math.exp(-t)/4/math.pi/t**2/rho
        return temp

    def uncollided_solution(self, xs, t):
        if self.uncollided == True:
            if self.geometry['slab'] == True:
                if self.source_type[0] == 1:
                    self.uncollided_solution_return = (self.plane_IC_uncollided_solution(xs, t) * self.source_strength)
                
                elif self.source_type[1] == 1:
                    self.uncollided_solution_return = (self.square_IC_uncollided_solution(xs, t) * self.source_strength)
            
                elif self.source_type[2] == 1:
                    if self.N_ang == 2:
                        self.uncollided_solution_return = (self.su_olson_s2_uncollided_solution(xs, t) * self.source_strength)
                    else:
                        self.uncollided_solution_return = (self.square_source_uncollided_solution(xs, t) * self.source_strength)
                    
                elif self.source_type[3] == 1:      
                    self.uncollided_solution_return =  (self.gaussian_IC_uncollided_solution(xs, t) * self.source_strength)        
                
                elif self.source_type[5] == 1:
                    if self.N_ang == 2:
                        self.uncollided_solution_return = (self.gaussian_s2(xs,t) * self.source_strength)
                    else:
                        self.uncollided_solution_return = (self.gaussian_source_uncollided_solution(xs, t) * self.source_strength)
                elif np.all(self.source_type==0):
                    self.uncollided_solution_return = self.fake_sedov(xs,t) * self.source_strength
            
            elif self.geometry['sphere'] == True:
                if self.source_type[0] == 1:
                    self.uncollided_solution_return = (self.point_source(xs, t) * self.source_strength)
                elif self.source_type[1] == 1:
                    self.uncollided_solution_return = self.shell_IC_uncollided_solution(xs, t) * self.source_strength
        else:
            self.uncollided_solution_return = np.zeros(xs.size)

        return self.uncollided_solution_return
        
# import quadpy
# import matplotlib.pyplot as plt

# M = 1
# N_ang = 256
# N_space = 4
# tfinal = 1.0
# x0 = 0.5
# move_type = np.array([1,0,0,0])
# sigma_t = np.ones(N_space)
# sigma_s = np.ones(N_space)
# time = True 
# plotting = True
# RMS = True
# uncollided = True
# moving = True
# ws = quadpy.c1.gauss_lobatto(N_ang).weights
# xs_quad = quadpy.c1.gauss_legendre(M+2).points
# ws_quad = quadpy.c1.gauss_legendre(M+2).weights
# mus = quadpy.c1.gauss_lobatto(N_ang).points
# source_type = np.array([0,0,1,0,0])
# initialize = build(N_ang, N_space, M, tfinal, x0, mus, ws, xs_quad, ws_quad, sigma_t, sigma_s, source_type, uncollided, moving, move_type, time, plotting, RMS)

# # source = source_class(initialize)
# uncollided_sol = uncollided_solution(initialize)
# # print(source.uncollided_solution(np.ones(1)*0.7, 1.0))
# xs = np.linspace(0,1.5,100)
# phi = uncollided_sol.uncollided_solution(xs,1.0)
# plt.plot(xs,phi)
# print(numba_expi(-1))
# print(numba_expi(-1.0))
# for i in range(xs.size):
#     print(expi(xs[i])- expi2(xs[i]))
# print(source.f1(1,0,0.5))

# print(0.5*((-0.5 + abs(0)) * expi(0.2-1) + math.exp(0.2 - 1)))
