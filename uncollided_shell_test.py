import numpy as np 
import math
import matplotlib.pyplot as plt
x0 = 0.5


def shell_IC_uncollided_solution(xs, t):
        temp = xs*0 
        sigma_abs = 0.0
        N01  = 4 * math.pi * x0**3 / 3
        N0 = 1 /N01
        n0 = N0 / (4. * np.pi / 3. * (x0 ** 3)) / (4. * np.pi) /2
        # n0 =   N0 * 3 / 4 / math.pi /x0**3 /4/math.pi
        for ix, r in enumerate(xs):
            tt = t + 1e-12
            mu_crit = min(1., max(-1.,0.5*(tt/r + r/tt - x0**2/(r*tt))))
            r2 = r ** 2 + t ** 2 - 2 * mu_crit * r * t
            # if np.sqrt(r2) < x0: 
            temp[ix] =  n0 * 2 * math.pi *  (1. - mu_crit ) * np.exp(-t * sigma_abs)
        return temp



tlist = [0.1, 0.5, 1,2,3,4,5,6,7,8,9,10,12,12,13,14]
for tf in tlist:
# tf = 1.0
    xs = np.linspace(0.000001,tf+x0)
    uncol = shell_IC_uncollided_solution(xs, tf)
    plt.plot(xs, uncol, '--', label = f't={tf}')
plt.legend()
plt.xlim(0,7)
plt.show()