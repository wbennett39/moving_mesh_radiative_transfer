import numpy as np
import matplotlib.pyplot as plt

from .loading_and_saving.load_solution import load_sol
from .plots.plot_functions import show
from .solver_classes import make_phi

def plot(tfinal, M,  N_space, problem_name, source_name, rad_or_transport, c, s2,
 cv0, x0_or_sigma , mat_or_rad, uncollided, moving, fign, name):
    data = load_sol(problem_name, source_name, rad_or_transport, c, s2, cv0)
    data.call_sol(tfinal, M, x0_or_sigma, N_space, mat_or_rad, uncollided, moving)
    N_ang = np.shape(data.coeff_mat)[0]
    print('loading solutions for', N_space, 'spaces', N_ang, 'angles', M, 'M')
    
    plt.ion()
    plt.figure(fign)
    plt.plot(data.xs, data.phi, 'k-')
    show(name)
    plt.show()
    plt.close()

    return data.xs, data.phi

# plot(1,4,16,'transport', 'square_s', 'transport', 1.0, False, 0.0, 0.5, 'rad', True, True )

def plot_thin_nonlinear_problems():
    tfinal_list = [1,10,20]


    plot(tfinal, M,  N_space, problem_name, source_name, rad_or_transport, c, s2,
 cv0, x0_or_sigma , mat_or_rad, uncollided, moving, fign, name)

