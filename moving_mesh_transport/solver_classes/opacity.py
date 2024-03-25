import numpy as np
import math
from .build_problem import build

from numba.experimental import jitclass
from numba import int64, float64, deferred_type, prange
from .functions import Pn, normPn
from numba import types, typed
import numba as nb

build_type = deferred_type()
build_type.define(build.class_type.instance_type)
kv_ty = (types.int64, types.unicode_type)
params_default = nb.typed.Dict.empty(key_type=nb.typeof('par_1'),value_type=nb.typeof(1))

data = [('N_ang', int64), 
        ('N_space', int64),
        ('M', int64),
        ('sigma_t', float64),
        ('sigma_s', float64),
        ('sigma_a', float64),
        ('mus', float64[:]),
        ('ws', float64[:]),
        ('x0', float64),
        ("xL", float64),
        ("xR", float64),
        ('sigma_func', nb.typeof(params_default)),
        ('Msigma', int64),
        ('AAA', float64[:,:,:]),
        ('xs_quad', float64[:]),
        ('ws_quad', float64[:]),
        ('edges', float64[:]),
        ('std', float64), 
        ('cs', float64[:,:]), 
        ('VV', float64[:]),
        ('VP', float64[:]),
        ('moving', float64),
        ('sigma_v', float64), 
        ('fake_sedov_v0', float64),

        ]


@ jitclass(data)
class sigma_integrator():
    def __init__(self, build):
        self.sigma_t = build.sigma_t
        self.sigma_s = build.sigma_s
        print(self.sigma_s,'sigma_s')
        self.sigma_a = self.sigma_t - self.sigma_s
        print(self.sigma_a,'sigma_a')
        self.sigma_func = build.sigma_func
        self.M = build.M
        self.Msigma = build.Msigma
        self.xs_quad = build.xs_quad
        self.ws_quad = build.ws_quad
        self.std = 2
        self.N_space = build.N_space
        self.edges = np.zeros(self.N_space + 1)
        self.cs = np.zeros((self.N_space, self.Msigma+ 1))
        self.VV = np.zeros(self.M+1)
        self.VP = np.zeros(self.M+1)
        self.AAA = np.zeros((self.M+1, self.M + 1, self.Msigma + 1))
        self.moving = False
        # if self.sigma_func['fake_sedov'] == True:
        #     self.moving = True
        # self.sigma_v = 0.005
        self.sigma_v = build.fake_sedov_v0

        # initialize integrals of Basis Legendre polynomials
        self.create_integral_matrices()

    def integrate_quad(self, a, b, i, j, k):
        argument = (b-a)/2 * self.xs_quad + (a+b)/2
        fact = np.sqrt(2*i + 1) * np.sqrt(2*j + 1) * np.sqrt(2*k + 1) / 2
        self.AAA[i,j,k] = fact * (b-a)/2 * np.sum(self.ws_quad *  Pn(i, argument, a, b) * Pn(j, argument, a, b) * Pn(k, argument, a, b))
    
    def integrate_moments(self, a, b, j, k, t):
        argument = (b-a)/2 * self.xs_quad + (a+b)/2
        self.cs[k, j] = (b-a)/2 * np.sum(self.ws_quad * self.sigma_function(argument, t) * normPn(j, argument, a, b))

    def both_even_or_odd(self, i, j, k):
        if i % 2 == 0:
            if (j + k) % 2 == 0:
                return True
            else:
                return False
        if i % 2 != 0:
            if (j + k) % 2 != 0:
                return True
            else:
                return False

    def create_integral_matrices(self):
        """
        creates a matrix with every integral over [-1,1] of the three normalized Legendre polynomials of order
        i, j, k. Each entry must be divided by sqrt(xR-xL) 
        """
        for i in range(self.M + 1):
            for j in range(self.M + 1):
                for k in range(self.Msigma + 1):
                    if (j + k >= i) and (self.both_even_or_odd(i, j, k)):
                        self.integrate_quad(-1, 1, i, j, k)
        # print(self.AAA)
    
    def sigma_moments(self, edges, t):
        for i in range(self.N_space):
            if (edges[i] != self.edges[i]) or (edges[i+1] != self.edges[i+1]) or self.moving == True :
                for j in range(self.Msigma + 1):
                    self.integrate_moments(edges[i], edges[i+1], j, i, t)
        self.edges = edges
    
    def xi2(self, x, t, x0, c1, v0tilde):
        return -x - c1 - v0tilde*(t)

    def heaviside(self,x):
        if x < 0.0:
            return 0.0
        else:
            return 1.0

    def heaviside_vector(self, x):
        return_array = np.ones(x.size)
        for ix, xx in enumerate(x):
            if xx < 0:
                return_array[ix] = 0.0
        return return_array

    def sigma_function(self, x, t):

        if self.sigma_func['constant'] == 1:
            return x * 0 + 1.0
        elif self.sigma_func['gaussian'] == 1:
            return np.exp(- x**2 /(2* self.std**2))  # probably shouldn't have sigma_a here
            # return x * 0 + 1.0
        elif self.sigma_func['siewert1'] == 1: # siewert with omega_0 = 1, s = 1
            return np.exp(-x - 2.5)
        elif self.sigma_func['siewert2'] == 1:
            return np.exp(-x/100000000000)
        elif self.sigma_func['fake_sedov'] == 1:
            # return np.exp(-(x- self.sigma_v * t)**2/(2*self.std**2))
            c1 = 1
            xi2x = self.xi2(x, t, 0, c1, self.sigma_v)
            rho2 = 0.2
            res = np.exp(-xi2x**2/self.std**2) * self.heaviside_vector(-xi2x - c1) + rho2*self.heaviside_vector(xi2x + c1)
            # vec_test = self.heaviside_vector(-xi2x - c1)
            # found = False
            # index = 0
            # if np.any(vec_test == 0):
            #     while found == False and index < x.size:
            #         if vec_test[index] == 1:
            #             found == True
            #             print(x[index], 'location of shock', t, 't')
            #             print(vec_test)
            #             print(-self.sigma_v*t - x[index])
            #             print("#--- --- --- --- --- --- ---#")
            #         index += 1


            return res
    
    def make_vectors(self, edges, u, space):

        # self.sigma_moments(edges) # take moments of the opacity
        xL = edges[space]
        xR = edges[space+1]
        dx = math.sqrt(xR-xL)
        if self.sigma_func['constant'] == True:
            self.VV = u * self.sigma_t
        else:
            for i in range(self.M + 1):
                for j in range(self.M + 1):
                    for k in range(self.Msigma + 1):
                        self.VV[i] +=   (self.sigma_a + self.sigma_s) * self.cs[space, k] * u[j] * self.AAA[i, j, k] / dx




    





