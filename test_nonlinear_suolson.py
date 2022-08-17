from moving_mesh_transport.benchmarks import integrate_greens as intg
from moving_mesh_transport.plots import plotting_script as plotter
from moving_mesh_transport import solver
import matplotlib.pyplot as plt
from moving_mesh_transport.plots.plot_square_s_times import main as plot_square_s_times

from moving_mesh_transport.solver_functions.run_functions import run
from numpy import genfromtxt

xs_mcclarren = genfromtxt('x.csv', delimiter=',')[1:,0]
T_mcclarren = genfromtxt('T.csv', delimiter=',')[1]
phi_mcclarren = genfromtxt('phi.csv', delimiter=',')[1]


run = run()
run.load('rad_transfer_const_cv')
print('loaded input script for nonlinear rad transfer')
run.square_source()
xs = run.xs
phi = run.phi
e = run.e

plt.figure(25)
plt.ion()
plt.plot(xs_mcclarren, phi_mcclarren, 'k-', label = 'phi from Dr. McClarren')
plt.plot(xs_mcclarren, T_mcclarren, 'r--', label = 'T from Dr. McClarren')
plt.plot(xs, phi, '-o', mfc='none', label ='my solution phi')
plt.plot(xs, e, '-^', mfc='none', label ='my solution e')

# plt.plot(xs, e/0.3, '-^', mfc='none', label ='my solution T')
plt.legend()
plt.show()