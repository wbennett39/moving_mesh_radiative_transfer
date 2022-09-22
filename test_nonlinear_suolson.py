from moving_mesh_transport.benchmarks import integrate_greens as intg
from moving_mesh_transport.plots import plotting_script as plotter
from moving_mesh_transport import solver
import matplotlib.pyplot as plt
from moving_mesh_transport.plots.plot_square_s_times import main as plot_square_s_times

from moving_mesh_transport.solver_functions.run_functions import run
from numpy import genfromtxt
from scipy.interpolate import interp1d
import numpy as np

# xs_mcclarren = genfromtxt('x_cv1.csv', delimiter=',')[1:]
# T_mcclarren = genfromtxt('T_cv1.csv', delimiter=',')[1:]
# phi_mcclarren = genfromtxt('phi_cv1.csv', delimiter=',')[1:]

xs_mcclarren = genfromtxt('x_cv1.csv', delimiter=',')[1:]
T_mcclarren = genfromtxt('T_cv1.csv', delimiter=',')[1:]
phi_mcclarren = genfromtxt('phi_cv1.csv', delimiter=',')[1:]

xs_thick_mcclarren = genfromtxt('x_cv1_thick.csv', delimiter=',')[1:]
T_thick_mcclarren = genfromtxt('T_cv1_thick.csv', delimiter=',')[1:]
phi_thick_mcclarren = genfromtxt('phi_cv1_thick.csv', delimiter=',')[1:]

T_mcclar_interp = interp1d(xs_mcclarren, T_mcclarren, kind = 'cubic')
phi_mcclar_interp = interp1d(xs_mcclarren, phi_mcclarren, kind = 'cubic')

run = run()
run.load('rad_transfer_const_cv')
print('loaded input script for nonlinear rad transfer')
run.square_source(True, True)
xs = run.xs
phi = run.phi
e = run.e
run.load('rad_transfer_const_cv_thick')
print('loaded input script for nonlinear rad transfer')
run.square_source(False, False)
xs_thick = run.xs
phi_thick = run.phi
e_thick = run.e





a = 0.01372
c = 299.98
cv = 0.03
cvbar = cv/a
print(a*c, 'ac')
plt.figure(25)
plt.ion()
plt.plot(xs, phi*a*c, 'o', mfc='none', label ='my solution phi')
plt.plot(xs, e/cvbar, '^', mfc='none', label ='my solution T')
plt.plot(xs_mcclarren, phi_mcclarren, 'k-', label = 'phi from Dr. McClarren')
plt.plot(xs_mcclarren, T_mcclarren, 'r--', label = 'T from Dr. McClarren')

plt.xlim(-2,2)
plt.show()

plt.figure(26)
plt.ion()
plt.plot(xs_thick, phi_thick*a*c, 'o', mfc='none', label ='my solution phi')
plt.plot(xs_thick, e_thick/cvbar, '^', mfc='none', label ='my solution T')

plt.plot(xs_thick_mcclarren, phi_thick_mcclarren, 'k-', label = 'phi from Dr. McClarren')
plt.plot(xs_thick_mcclarren, T_thick_mcclarren, 'r--', label = 'T from Dr. McClarren')

plt.xlim(350, 450)
plt.show()

# plt.plot(xs, e/0.3, '-^', mfc='none', label ='my solution T')

print("### ### ### ### ### ### ### ###")
RMSE = np.sqrt(np.mean((phi*a*c- phi_mcclar_interp(xs))**2))
print(RMSE)
print("### ### ### ### ### ### ### ###")

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