import sys
sys.path.append('/Users/bennett/Documents/Github/transport_benchmarks/')
import matplotlib.pyplot as plt
from moving_mesh_transport.plots.plot_functions.show import show
from moving_mesh_transport.loading_and_saving.load_solution import load_sol
import numpy as np
import math
from scipy.special import erf
import scipy.integrate as integrate
from tqdm import tqdm
from benchmarks import integrate_greens as intg
from scipy.interpolate import interp1d as interp


def RMSE(l1, l2):
    diff = (l1-l2)**2
    return np.sqrt(np.mean(diff))

def shell_source_convergence():
    tfinal_list = [1.0, 5.0]
    spaces_list = np.array([[3, 5, 10, 15], [10 , 10, 10, 10]])
    RMSE_list = np.zeros((spaces_list.shape))
    Ms_list = [3,3]
    res = intg.shell_source(1.0, 5)
    plt.close()
    loader = load_sol('transport', 'square_IC', 'transport', s2 = False, file_name = 'run_data.hdf5')

    for it, tt in enumerate(tfinal_list):
        for ix, space in enumerate(spaces_list[it]):
            loader.call_sol(tt, Ms_list[it], 0.5, space, 'rad', False, True, 1.0)
            phi = loader.phi
            xs = loader.xs
            bench = intg.shell_source(tt, 150, choose_xs = True, xpnts = xs)
            plt.close()
            interp_bench = interp(bench[0], bench[1] + bench[2])
            interp_bench_uncol = interp(bench[0],  bench[2])
            RMSE_list[it, ix] = RMSE(phi, bench[1] + bench[2])
            plt.figure(12)
            # plt.plot(xs, interp_bench(xs), 'k-')
            plt.plot(xs, np.abs(phi-interp_bench(xs)), 'o', mfc = 'none')
            plt.figure(13)
            plt.plot(xs, shell_IC_uncollided_solution(xs, tt) - interp_bench_uncol(xs), 'o', mfc = 'none')
            # plt.plot(bench[0], bench[2], 'k-')
        plt.figure(it)
        plt.loglog(spaces_list[it, :], RMSE_list[it, :], '-o')
        # plt.loglog(np.array(spaces_list[it]), np.array(spaces_list[it])**(-1.0))
        plt.show()


def shell_IC_uncollided_solution(xs, t):
        temp = xs*0 
        x0 = 0.5
        sigma_abs = 1
        N01  = 4 * math.pi * x0**3 / 3
        N0 = 1
        n0 = N0 / (4. * math.pi / 3. * (x0 ** 3)) / (4. * math.pi)
        for ix, r in enumerate(xs):
            tt = t + 1e-12
            mu_crit = min(1., max(-1.,0.5*(tt/r+r/tt-x0**2/(r*tt))))
            r2 = r ** 2 + t ** 2 - 2 * mu_crit * r * t
            # if np.sqrt(r2) < self.x0: 
            temp[ix] = n0 * 2 * math.pi *  (1. - mu_crit ) * np.exp(-t * sigma_abs) 
        return temp







       
       
