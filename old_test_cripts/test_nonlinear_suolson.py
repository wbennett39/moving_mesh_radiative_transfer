# from moving_mesh_transport.benchmarks import integrate_greens as intg
# from moving_mesh_transport.plots import plotting_script as plotter
# from moving_mesh_transport import solver
# import matplotlib.pyplot as plt
# from moving_mesh_transport.plots.plot_square_s_times import main as plot_square_s_times
# from moving_mesh_transport.solution_plotter import plot as sol_plot

# from moving_mesh_transport.solver_functions.run_functions import run
# from numpy import genfromtxt
# from scipy.interpolate import interp1d
# import numpy as np

# a = 0.01372
# c = 299.98
# cv = 0.03
# cvbar = cv/a
# cvbar2 = cv/a/800


# # load sweep 2d results

# # sigma = 1, t = 1/c x0 = 0.5
# xs_mcclarren = genfromtxt('sweep2d_results/x_cv1.csv', delimiter=',')[1:]
# T_mcclarren = genfromtxt('sweep2d_results/T_cv1.csv', delimiter=',')[1:]
# phi_mcclarren = genfromtxt('sweep2d_results/phi_cv1.csv', delimiter=',')[1:]


# # sigma = 1, t = 31.6228/c x0 = 0.5
# xs_mcclarren_31 = genfromtxt('sweep2d_results/x_cv31.csv', delimiter=',')[1:]
# T_mcclarren_31 = genfromtxt('sweep2d_results/T_cv31.csv', delimiter=',')[1:]
# phi_mcclarren_31 = genfromtxt('sweep2d_results/phi_cv31.csv', delimiter=',')[1:]


# # sigma = 1, t = 100/c, x0 = 400.0
# xs_thick_mcclarren = genfromtxt('sweep2d_results/x_cv1_thick.csv', delimiter=',')[1:]
# T_thick_mcclarren = genfromtxt('sweep2d_results/T_cv1_thick.csv', delimiter=',')[1:]
# phi_thick_mcclarren = genfromtxt('sweep2d_results/phi_cv1_thick.csv', delimiter=',')[1:]

# gaussian_xs_thick_mc = genfromtxt('sweep2d_results/gaussian_x_cv1_thick.csv', delimiter=',')[1:]
# gaussian_T_thick_mc = genfromtxt('sweep2d_results/gaussian_T_cv1_thick.csv', delimiter=',')[1:]
# gaussian_phi_thick_mc = genfromtxt('sweep2d_results/gaussian_phi_cv1_thick.csv', delimiter=',')[1:]

# # sigma = 800, t = 0.125/c, x0 = 0.5
# xs_thick_opacity_mcclarren = genfromtxt('sweep2d_results/x_cv1_thick_sigma400.csv', delimiter=',')[1:] * 800.0
# T_thick_opacity_mcclarren = genfromtxt('sweep2d_results/T_cv1_thick_sigma400.csv', delimiter=',')[1:] 
# phi_thick_opacity_mcclarren = genfromtxt('sweep2d_results/phi_cv1_thick_sigma400.csv', delimiter=',')[1:]


# T_mcclar_interp = interp1d(xs_mcclarren, T_mcclarren, kind = 'cubic')
# phi_mcclar_interp = interp1d(xs_mcclarren, phi_mcclarren, kind = 'cubic')

# # load my results


# thin_nonlinear_plot = sol_plot(1.0, 5,  256, 'rad_transfer_const_cv', 'square_s', 'rad', 0.0, False, 
#             0.03, 0.5, 'rad', True, True, 1, 'garbage_plot.pdf', mkr1 = 'ob', mkr2 = '^b')

# thin_nonlinear_plot_31 = sol_plot(31.6228, 4,  32, 'transfer_const_cv=0.03', 'square_s', 'rad', 0.0, False, 
#             0.03, 0.5, 'rad', True, True, 1, 'garbage_plot.pdf', mkr1 = 'ob', mkr2 = '^b')


# # thick_source_nonlinear_plot = sol_plot(100.0, 10,  32, 'rad_transfer_const_cv_thick', 'square_s', 'rad', 0.0, False, 
# #             0.03, 400.0 , 'rad', False, False, 2, 'garbage_plot.pdf',mkr1 = 'ob', mkr2 = '^b' )

# # thick_opacity_nonlinear_plot = sol_plot(0.125, 0,  64, 'rad_transfer_const_cv_thick', 'square_s', 'rad', 0.0, False, 
# #             0.03, 0.5 , 'rad', False, False, 2, 'garbage_plot.pdf', mkr1 = 'ob', mkr2 = '^b')

# # gaussian_thick_nonlinear_plot = sol_plot(100.0, 15,  16, 'rad_transfer_const_cv_thick', 'gaussian_s', 'rad', 0.0, False, 
# #             0.03, 300.0 , 'rad', False, False, 2, 'garbage_plot.pdf', mkr1 = 'ob', mkr2 = '^b')

# plt.close()
# plt.close()
# plt.close()
# plt.close()

# thin_xs, thin_nonlinear_phi = thin_nonlinear_plot.plot() 
# thin_xs_31, thin_nonlinear_phi_31 = thin_nonlinear_plot_31.plot() 
# # thick_xs, thick_nonlinear_phi = thick_source_nonlinear_plot.plot()
# # thick_xs_gaussian, thick_nonlinear_phi_gaussian = gaussian_thick_nonlinear_plot.plot()
# # thick_opacity_xs, thick_opacity_nonlinear_phi = thick_opacity_nonlinear_plot.plot()

# thin_nonlinear_e = thin_nonlinear_plot.e

# thin_nonlinear_e_31 = thin_nonlinear_plot_31.e
# # thick_nonlinear_e = thick_source_nonlinear_plot.e
# # thick_gaussian_nonlinear_e = gaussian_thick_nonlinear_plot.e
# # thick_opacity_e = thick_opacity_nonlinear_plot.e 

# thin_phi = thin_nonlinear_phi * a * c

# thin_phi_31 = thin_nonlinear_phi_31 * a * c
# # thick_phi = thick_nonlinear_phi * a * c
# # thick_opacity_phi = thick_opacity_nonlinear_phi * a * c
# # thick_nonlinear_phi_gaussian = thick_nonlinear_phi_gaussian * a * c


# thin_T = thin_nonlinear_e / cvbar
# thin_T_31 = thin_nonlinear_e_31 / cvbar
# # thick_T_gaussian = thick_gaussian_nonlinear_e / cvbar
# # thick_T = thick_nonlinear_e / cvbar
# # thick_opacity_T = thick_opacity_e / cvbar

# # thick_opacity_xs = thick_opacity_xs * 800.0
# # thick_xs_gaussian = thick_xs_gaussian / 800.0


# # run = run()
# # run.load('rad_transfer_const_cv')
# # print('loaded input script for nonlinear rad transfer')
# # run.square_source(True, True)
# # xs = run.xs
# # phi = run.phi
# # e = run.e
# # run.load('rad_transfer_const_cv_thick')
# # print('loaded input script for nonlinear rad transfer')
# # run.square_source(False, False)
# # xs_thick = run.xs
# # phi_thick = run.phi
# # e_thick = run.e



# print(a*c, 'ac')
# plt.figure(25)
# plt.ion()
# plt.title('t = 1 cv = 0.03')
# plt.plot(thin_xs, thin_phi, 'o', mfc='none', label ='my solution phi')
# plt.plot(thin_xs, thin_T, '^', mfc='none', label ='my solution T')
# plt.plot(xs_mcclarren, phi_mcclarren, 'k-', label = 'phi sweep 2d')
# plt.plot(xs_mcclarren, T_mcclarren, 'r--', label = 'T sweep 2d')

# plt.xlim(-2,2)
# plt.legend()
# plt.show()

# plt.figure(26)
# plt.ion()
# plt.title('t = 31.6228 cv = 0.03')
# plt.plot(thin_xs_31, thin_phi_31, 'o', mfc='none', label ='my solution phi')
# plt.plot(thin_xs_31, thin_T_31, '^', mfc='none', label ='my solution T')
# plt.plot(xs_mcclarren_31, phi_mcclarren_31, 'k-', label = 'phi sweep 2d')
# plt.plot(xs_mcclarren_31, T_mcclarren_31, 'r--', label = 'T sweep 2d')


# # plt.plot(thick_xs, thick_phi, 'o', mfc='none', label ='my solution phi')
# # plt.plot(thick_xs, thick_T, '^', mfc='none', label ='my solution T')

# # plt.plot(thick_opacity_xs, thick_opacity_phi, 'o', mfc='none', label ='my solution phi, sigma = 800')
# # plt.plot(thick_opacity_xs, thick_opacity_T, '^', mfc='none', label ='my solution T, sigma = 800')


# # plt.plot(xs_thick_mcclarren, phi_thick_mcclarren, 'k-', label = 'phi sweep 2d ')
# # plt.plot(xs_thick_mcclarren, T_thick_mcclarren, 'r--', label = 'T sweep 2d')

# # plt.plot(xs_thick_opacity_mcclarren, phi_thick_opacity_mcclarren, 'b-', label = 'phi sweep 2d sigma = 800 ')
# # plt.plot(xs_thick_opacity_mcclarren, T_thick_opacity_mcclarren, 'g--', label = 'T sweep 2d sigma = 800')

# # plt.xlim(350, 450)
# # plt.show()

# # plt.figure(27)

# # plt.plot(thick_xs_gaussian, thick_nonlinear_phi_gaussian, 'o', mfc = 'none', label= 'my code phi')
# # plt.plot(thick_xs_gaussian, thick_T_gaussian, 'o', mfc = 'none', label= 'my code T')

# # plt.plot(gaussian_xs_thick_mc, gaussian_phi_thick_mc, 'k-', label = 'phi sweep 2d ')
# # plt.plot(gaussian_xs_thick_mc, gaussian_T_thick_mc, 'r--', label = 'T sweep 2d')
# # plt.xlim(-2,2)
# # plt.show()



# # plt.plot(xs, e/0.3, '-^', mfc='none', label ='my solution T')

# # print("### ### ### ### ### ### ### ###")
# # RMSE = np.sqrt(np.mean((phi*a*c- phi_mcclar_interp(xs))**2))
# # print(RMSE)
# # print("### ### ### ### ### ### ### ###")

# # T_rat = max(e/cvbar)/max(T_mcclarren)
# # phi_rat = max(phi*a*c)/max(phi_mcclarren)
# # print(T_rat,  'my T over McClarrens T')
# # print(1/T_rat,  'McClarrens T over my T')
# # print(phi_rat,  'my PHI over McClarrens phi')
# # print(1/phi_rat,  'McClarrens phi over my phi')

# # print(max(phi), 'my phi max')
# # print(max(e/cvbar), 'my T max')
# # print(max(T_mcclarren), 'McClarren T max')
# # print(max(phi_mcclarren), 'McClarren phi max')

# plt.legend()
# plt.show()