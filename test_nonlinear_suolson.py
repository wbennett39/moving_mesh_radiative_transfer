from moving_mesh_transport.benchmarks import integrate_greens as intg
from moving_mesh_transport.plots import plotting_script as plotter
from moving_mesh_transport import solver
import matplotlib.pyplot as plt
from moving_mesh_transport.plots.plot_square_s_times import main as plot_square_s_times

from moving_mesh_transport.solver_functions.run_functions import run
from numpy import genfromtxt

# xs_mcclarren = genfromtxt('x_cv1.csv', delimiter=',')[1:]
# T_mcclarren = genfromtxt('T_cv1.csv', delimiter=',')[1:]
# phi_mcclarren = genfromtxt('phi_cv1.csv', delimiter=',')[1:]

xs_mcclarren = genfromtxt('x.csv', delimiter=',')[1:]
T_mcclarren = genfromtxt('T.csv', delimiter=',')[1:]
phi_mcclarren = genfromtxt('phi.csv', delimiter=',')[1:]


run = run()
run.load('rad_transfer_const_cv')
print('loaded input script for nonlinear rad transfer')
run.square_source(True, True)
xs = run.xs
phi = run.phi
e = run.e
a = 0.01372
c = 299.98
cv = 0.3
cvbar = cv/a
print(a*c, 'ac')
plt.figure(25)
plt.ion()
plt.plot(xs_mcclarren, phi_mcclarren, 'k-', label = 'phi from Dr. McClarren')
plt.plot(xs_mcclarren, T_mcclarren, 'r--', label = 'T from Dr. McClarren')
plt.plot(xs, phi*a*c, '-o', mfc='none', label ='my solution phi')
plt.plot(xs, e/cvbar, '-^', mfc='none', label ='my solution T')


# plt.plot(xs, e/0.3, '-^', mfc='none', label ='my solution T')
T_rat = max(e/cvbar)/max(T_mcclarren)
phi_rat = max(phi*a*c)/max(phi_mcclarren)
print(T_rat,  'my T over McClarrens T')
print(1/T_rat,  'McClarrens T over my T')
print(phi_rat,  'my PHI over McClarrens phi')
print(1/phi_rat,  'McClarrens phi over my phi')

print(max(phi), 'my phi max')
print(max(e/cvbar), 'my T max')
print(max(T_mcclarren), 'McClarren T max')
print(max(phi_mcclarren), 'McClarren phi max')

plt.legend()
plt.show()