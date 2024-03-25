import matplotlib.pyplot as plt
import numpy as np


def xi2(x, t, x0, c1, v0tilde):
    return -x - c1 - v0tilde*(t)

def heaviside_vector(x):
    return_array = np.ones(x.size)
    for ix, xx in enumerate(x):
        if xx < 0:
            return_array[ix] = 0.0
    return return_array

def opacity(x,t,std, sigma_v):
   # return np.exp(-(x- self.sigma_v * t)**2/(2*self.std**2))
        c1 = 1
        xi2x = xi2(x, t, 0, c1, sigma_v)
        rho2 = 0.2
        res = np.exp(-xi2x**2/std**2) * heaviside_vector(-xi2x - c1) + rho2*heaviside_vector(xi2x + c1)
        return res



def plotter(t):
    # plt.ion()
    xs = np.linspace(-5,5,500)
    std = 2
    sigma_v = 0.05
    y = opacity(xs, t, std, sigma_v)
    # plt.plot(xs,y)
    # plt.show()