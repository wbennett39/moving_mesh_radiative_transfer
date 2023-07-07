import numpy as np
import scipy.special as spec
import math


def exact_sol(alpha = 0, beta = 4, npnts = 100, omega = 0, lambda_ = 0, mu = 1):
    m = omega*(1-mu)
    k = omega * (1+lambda_)
    n = (4 + alpha - beta) / beta
    p = 2-k-m+(1-m)*n
    print('p = ', p)
    if 2 - k - m >0:
        l = (1-m)/ (2-k-m)
    elif 2- k - m <0:
        l = -1/n -(d-m)/ (2-k-m)
    if n > 0:
        xi0 = (p*abs(2-k-m)**(n+1) / (n*spec.beta(l, 1/n+1)**n))**(1/p) 

        xi = np.linspace(0, xi0, npnts)

        temp = (n * (xi0**(2-k-m)-xi**(2-k-m)) / (p * (2-k-m)))**(1/n)
    elif omega == k == m == 0:
        xi = np.linspace(0, 4, npnts)

        temp = 0.5 * (1/math.sqrt(math.pi)) * np.exp(-xi**2/4)

    else:
        print('case not yet implemented')
        assert 0
    return xi, temp, m, k, n, p

