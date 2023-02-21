import numpy as np

from .build_problem import build

from numba.experimental import jitclass
from numba import int64, float64, deferred_type, prange
from .functions import Pn, normPn

build_type = deferred_type()
build_type.define(build.class_type.instance_type)
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
        ('sigma_func', int64[:]),
        ('Msigma', int64),
        ('AAA', float64[:,:,:]),
        ('xs_quad', float64[:]),
        ('ws_quad', float64[:]),
        ('edges', float64[:]),
        ('std', float64), 
        ('cs', float64[:,:]), 
        ('VV', float64[:]),
        ('VP', float64[:])
        ]


@ jitclass(data)
class sigma_integrator():
    def __init__(self, build):
        self.sigma_t = build.sigma_t
        self.sigma_s = build.sigma_s
        self.sigma_a = self.sigma_t - self.sigma_s
        self.sigma_func = build.sigma_func
        self.M = build.M
        self.Msigma = build.Msigma
        self.xs_quad = build.xs_quad
        self.ws_quad = build.ws_quad
        self.std = 0.5
        self.N_space = build.N_space
        self.edges = np.zeros(self.N_space + 1)
        self.cs = np.zeros((self.N_space, self.Msigma+ 1))
        self.VV = np.zeros(self.M+1)
        self.VP = np.zeros(self.M+1)
        self.AAA = np.zeros((self.M+1, self.M + 1, self.Msigma + 1))

        # initialize integrals of Basis Legendre polynomials
        self.create_integral_matrices()

    def integrate_quad(self, a, b, i, j, k):
        argument = (b-a)/2 * self.xs_quad + (a+b)/2
        fact = np.sqrt(2*i + 1) * np.sqrt(2*j + 1) * np.sqrt(2*k + 1) / 2
        self.AAA[i,j,k] = fact * (b-a)/2 * np.sum(self.ws_quad *  Pn(i, argument, a, b) * Pn(j, argument, a, b) * Pn(k, argument, a, b))
    
    def integrate_moments(self, a, b, j, k):
        argument = (b-a)/2 * self.xs_quad + (a+b)/2
        self.cs[k, j] = (b-a)/2 * np.sum(self.ws_quad * self.sigma_function(argument) * normPn(j, argument, a, b))

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
                    # if (j + k >= i) and (self.both_even_or_odd(i, j, k)):
                    self.integrate_quad(-1, 1, i, j, k)
        print(self.AAA)
    
    def sigma_moments(self, edges):
        for i in range(self.N_space):
            if (edges[i] != self.edges[i]) or (edges[i+1] != self.edges[i+1]):
                for j in range(self.Msigma + 1):
                    self.integrate_moments(edges[i], edges[i+1], j, i)
                print(self.cs, 'cs')
        self.edges = edges
        

    def sigma_function(self, x):
        if self.sigma_func[0] == 1:
            return x * 0 + 1.0
        elif self.sigma_func[1] == 1:
            # return np.exp(- x**2 /(2* self.std**2)) * self.sigma_a
            return x * 0 + 1.0
    
    def make_vectors(self, edges, u, space):
        VV = u * 0
        self.sigma_moments(edges) # take moments of the opacity
        xL = edges[space]
        xR = edges[space+1]
        dx = np.sqrt(xR-xL)
        
        for i in range(self.M + 1):
            for j in range(self.M + 1):
                for k in range(self.Msigma + 1):
                    VV[i] += self.cs[space, k] * u[j] * self.AAA[i, j, k] / dx
        return VV, self.VP



    





