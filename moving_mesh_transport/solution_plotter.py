# import numpy as np
# import matplotlib.pyplot as plt

# from .loading_and_saving.load_solution import load_sol
# from .plots.plot_functions.show import show
# from .plots.plotting_script import plot_coefficients
# from .solver_classes import make_phi, build_problem
# from .solver_classes.uncollided_solutions import uncollided_solution
# from .loading_and_saving.load_bench import load_bench
# from pathlib import Path
# from scipy.interpolate import interp1d
# import quadpy
# import csv
# from .misc_functions import source_name_finder

# def trunc(values, decs=0):
#     return np.trunc(values*10**decs)/(10**decs)

# def sparsify(xlist, ylist, num = 50):
#     print('calling sparsify')
#     def every_other(xlist, ylist):
#         start = np.argmin(np.abs(xlist))
#         xlist_new = np.linspace(xlist[0], xlist[-1])
#         yint = interp1d(xlist, ylist, 'cubic')
#         ylist_new = yint(xlist_new)
#         return np.array(xlist_new), np.array(ylist_new)

#     xlist = np.array(xlist)
#     ylist = np.array(ylist)
    
#     while xlist.size > num:
#         assert(xlist.size == ylist.size)
#         xlist, ylist = every_other(xlist, ylist)
#     return xlist, ylist




# su_tfinal_list = np.array([0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228,  100.0])
# su_xs_list = np.array([0.01, 0.1, 0.17783, 0.31623, 0.45, 0.5, 0.56234, 0.75, 1.0, 1.33352, 1.77828, 3.16228, 5.62341, 10.0, 17.78279])
        
# class plot:
#     def __init__(self,tfinal, M,  N_space, problem_name, source_name, rad_or_transport, c, s2,
#                 cv0, x0_or_sigma , mat_or_rad, uncollided, moving, fign, name, mkr1='k-', mkr2='k--', mkr3 = 'k:', mkr4 = '-.k', mkr5 = ":k", mkr6 = "-*k", file_name = 'run_data_crc_dec15-3.hdf5'):
#         self.tfinal = tfinal
#         self.M = M
#         self.problem_name = problem_name
#         self.source_name = source_name
#         self.rad_or_transport = rad_or_transport
#         self.c = c
#         self.s2 = s2
#         self.cv0=cv0
#         self.x0_or_sigma = x0_or_sigma
#         self.mat_or_rad = mat_or_rad
#         self.uncollided = uncollided
#         self.moving=moving
#         self.fign = fign
#         self.name = name
#         self.N_space = N_space
#         self.mkr1 = mkr1
#         self.mkr2 = mkr2
#         self.mkr3 = mkr3
#         self.mkr4 = mkr4
#         self.mkr5 = mkr5
#         self.mkr6 = mkr6


#         if self.problem_name == 'su_olson_thick_s2':
#             self.file_name = 'run_data_crc_dec7-4.hdf5'
#         else:
#             if (self.source_name == 'square_s' or self.x0_or_sigma == 0.375):
#                 self.file_name = file_name
#             # elif self.problem_name == 'transfer_const_cv=0.03':
#             #     self.file_name = 'run_data_crc_nov15.hdf5'
#             elif self.tfinal >= 10.0 and self.source_name == 'gaussian_s' and self.cv0 == 0.03:
#                 self.file_name = 'run_data_crc_dec7-4.hdf5'
#             else:
#                 self.file_name = 'run_data_crc_nov23.hdf5'


#     def load(self):
#         data = load_sol(self.problem_name, self.source_name, self.rad_or_transport, self.c, self.s2, self.cv0, file_name = self.file_name)
#         data.call_sol(self.tfinal, self.M, self.x0_or_sigma, self.N_space, self.mat_or_rad, self.uncollided, self.moving)
#         self.N_ang = np.shape(data.coeff_mat)[0] -1
#         self.coeff_mat = data.coeff_mat
#         self.ws = data.ws
#         self.edges = data.edges
#         self.e = data.e
#         self.phi = data.phi
#         self.xs = data.xs

#         print('loading solutions for', self.N_space, 'spaces',self.N_ang, 'angles', self.M, 'M')
        
#         return data.xs, data.phi
    
#     def plot(self, uncollided_class):
#         plt.ion()
#         plt.figure(self.fign)

#         middle = np.argmin(np.abs(self.xs))
#         xs_plot = self.xs[middle:]
#         print(xs_plot[-1], 'edge')
#         phi_plot = self.phi[middle:]
#         e_plot = self.e[middle:]
    
#         if self.s2 == True:
#             xs_plot = - xs_plot
#             # marker1 = self.mkr4
#             # marker2 = self.mkr5
#             # marker3 = self.mkr6
#         # else:
#         marker1 = self.mkr1
#         marker2 = self.mkr2
#         marker3 = self.mkr3
#         if self.problem_name in ['transfer_const_cv=0.03', 'transfer_const_cv=0.03_s2', 'transfer_const_cv=0.03_thick', 'transfer_const_cv=0.03_thick_s2', 'transfer_const_cv=0.03_s2']:
#             a = 0.01372
#             cvbar = 0.03/a
#             T = e_plot / cvbar
#             plt.plot(xs_plot, np.power(phi_plot,.25), marker1, mfc = 'none')
#             plt.plot(xs_plot, T, marker2, mfc = 'none')
#             maxT = max(T)
#             maxphi = max(phi_plot**.25)
#             height = max(maxT, maxphi)
#             uncollided_sol = uncollided_class.uncollided_solution(xs_plot, self.tfinal)
#             if max(uncollided_sol > 1e-8):
#                 plt.plot(xs_plot, np.power(uncollided_sol, 0.25), 'k-.', mfc = 'none')

#         elif self.problem_name in ['su_olson', 'su_olson_s2', 'su_olson_thick', 'su_olson_thick_s2']:
#             T = np.power(e_plot, .25)
#             maxe = max(e_plot)
#             maxphi = max(phi_plot)
#             height = max(maxe, maxphi)
#             plt.plot(xs_plot, phi_plot, marker1, mfc = 'none')
#             plt.plot(xs_plot, e_plot, marker2, mfc = 'none')
#             # sparse_xs = np.linspace(0, xs_plot[-1], 25)
#             uncollided_sol = uncollided_class.uncollided_solution(xs_plot, self.tfinal)
#             if max(uncollided_sol > 1e-6):
#                 plt.plot(xs_plot, uncollided_sol, 'k-.', mfc = 'none')
#         # plt.plot(xs_plot, T, marker3, mfc = 'none')
#         if self.s2 == False:
#             plt.xlim(-xs_plot[-1] - 0.05, xs_plot[-1] + 0.05)
#             left = -xs_plot[-1]/3
#             if self.problem_name == 'su_olson_thick':
#                 plt.xlim(-1.5, 1.5)
#                 left = -1.5/3
#             elif self.problem_name == 'su_olson':
#                 if self.source_name == 'square_s':
#                     if self.tfinal == 100.0:
#                         plt.xlim(-25,25)
#                         left = -25/3
#             # elif self.problem_name == 'transfer_const_cv=0.03':
#             #     if self.source_name == 'gaussian_s':
#             #         if self.tfinal == 100:
#             #             plt.xlim(-16,16)
#             #             left = -16/3

            
#             # if self.tfinal == 100.0 or self.tfinal == 30.0:
            
#             right = -left
#             txtheight = height * 1.05
#             plt.text(left, txtheight, r'$S_2$', fontsize = 'xx-large', horizontalalignment = 'center' )
#             plt.text(right, txtheight, 'Transport', fontsize = 'xx-large', horizontalalignment = 'center' )
#             plt.axvline(x = 0, color = 'tab:gray')
        
#         # diffusion results

#         if self.problem_name in ['su_olson_thick'] and self.source_name == 'square_s' and self.s2 == False:
#             xvals = []
#             T_vals = []
#             with open(f'diffusion_results/t_diff_t={self.tfinal}.csv', 'r') as file:
#                 reader = csv.reader(file)
#                 for i, line in enumerate(reader):
#                     if i > 0:
#                         xvals.append(float(line[1]))
#                         T_vals.append(float(line[0]))
#             xvals, T_vals = sparsify(xvals, T_vals)
#             plt.plot(xvals, T_vals, 'k^', mfc = 'none')
        
#         elif self.problem_name in ['su_olson_thick'] and self.source_name == 'gaussian_s' and self.s2 == False:
#             xvals = []
#             T_vals = []
#             with open(f'diffusion_results/t_diff_gauss_t={self.tfinal}.csv', 'r') as file:
#                 reader = csv.reader(file)
#                 for i, line in enumerate(reader):
#                     if i > 0:
#                         xvals.append(float(line[1]))
#                         T_vals.append(float(line[0]))
#             xvals, T_vals = sparsify(xvals, T_vals)
#             plt.plot(xvals, T_vals, 'k^', mfc = 'none')

#         elif self.problem_name in ['transfer_const_cv=0.03_thick'] and self.source_name == 'gaussian_s' and self.s2 == False:
#             xvals = []
#             T_vals = []
#             with open(f'diffusion_results/t_diff_gauss_nonlin_t={self.tfinal}.csv', 'r') as file:
#                 reader = csv.reader(file)
#                 for i, line in enumerate(reader):
#                     if i > 0:
#                         xvals.append(float(line[1]))
#                         T_vals.append(float(line[0]))
#             xvals, T_vals = sparsify(xvals, T_vals)
#             plt.plot(xvals, T_vals, 'k^', mfc = 'none')


#             # data_phi_trunc = trunc(data_phi, 15)



#         if self.s2 == True:
#             show(self.name)
#             plt.show()
#             plt.close()

# # plot(1,4,16,'transport', 'square_s', 'transport', 1.0, False, 0.0, 0.5, 'rad', True, True )

# def plot_thin_nonlinear_problems(M=12, N_space = 64, problem_name = 'rad_transfer_const_cv', rad_or_transport = 'rad', 
#                                 c = 0.0, s2 = False, cv0=0.03, x0_or_sigma = 0.5, mat_or_rad ='rad', uncollided = True, moving = False):
#     tfinal_list = [0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228,100.0]
#     source_name_list = ['gaussian_s']
#     fign = 1
#     delim = '_'
#     path = Path("moving_mesh_transport")
#     for tfinal in tfinal_list:
#         for source_name in source_name_list:
            
#             string = problem_name + delim + str(tfinal) + delim + source_name + delim + 'x0='+ str(x0_or_sigma) + delim + 'cv0=' + str(cv0) 

#             name = str(path / 'plots' / 'solution_plots') + '/' + string

#             plotter = plot(tfinal, M,  N_space, problem_name, source_name, rad_or_transport, c, s2, 
#             cv0, x0_or_sigma , mat_or_rad, uncollided, moving, fign, name)

#             fign +=1

#         plotter.plot()

# def plot_thin_nonlinear_problems_s2(M=10, N_space = 128, problem_name = 'transfer_const_cv=0.03_s2', rad_or_transport = 'rad', 
#                                 c = 0.0, s2 = False, cv0=0.03, x0_or_sigma = 0.5, mat_or_rad ='rad', uncollided = True, moving = False):
#     # tfinal_list = [0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228,100.0]
#     tfinal_list = [10.0, 31.6228,100.0]

#     source_name_list = ['gaussian_s']
#     fign = 1
#     delim = '_'
#     path = Path("moving_mesh_transport")
#     for tfinal in tfinal_list:
#         for source_name in source_name_list:
            
#             string = problem_name + delim + str(tfinal) + delim + source_name + delim + 'x0='+ str(x0_or_sigma) + delim + 'cv0=' + str(cv0) 

#             name = str(path / 'plots' / 'solution_plots') + '/' + string

#             plotter = plot(tfinal, M,  N_space, problem_name, source_name, rad_or_transport, c, s2, 
#             cv0, x0_or_sigma , mat_or_rad, uncollided, moving, fign, name)

#             fign +=1

#         plotter.plot()


# def plot_thick_nonlinear_problems(M=10, N_spaces = [64, 64, 128], problem_name = 'rad_transfer_const_cv_thick', rad_or_transport = 'rad', 
#                                 c = 0.0, s2 = False, cv0=0.03, x0_or_sigma = 0.375, mat_or_rad ='rad', uncollided = False, moving = False):
    
#     tfinal_list = [0.3, 3.0, 30.0]
#     source_name_list = ['gaussian_s']
#     fign = 1
#     delim = '_'
#     path = Path("moving_mesh_transport")
#     for count, tfinal in  enumerate(tfinal_list):
            
#         string = problem_name + delim + str(tfinal) + delim + source_name_list[0] + delim + 'x0='+ str(x0_or_sigma) + delim + 'cv0=' + str(cv0) 

#         name = str(path / 'plots' / 'solution_plots') + '/' + string

#         plotter = plot(tfinal, M,  N_spaces[count], problem_name, source_name_list[0], rad_or_transport, c, s2, 
#         cv0, x0_or_sigma , mat_or_rad, uncollided, moving, fign, name)

#         fign +=1

#         plotter.plot()

  


# def plot_thick_nonlinear_problems_s2(M=10, N_spaces = [64, 64, 128], problem_name = 'rad_transfer_const_cv_thick_s2', rad_or_transport = 'rad', 
#                                 c = 0.0, s2 = False, cv0=0.03, x0_or_sigma = 400.0, mat_or_rad ='rad', uncollided = False, moving = False):
    
#     tfinal_list = [0.3, 3.0, 3.0]
#     source_name_list = ['gaussian_s']
#     fign = 1
#     delim = '_'
#     path = Path("moving_mesh_transport")
#     for count, tfinal in  enumerate(tfinal_list):
            
#         string = problem_name + delim + str(tfinal) + delim + source_name_list[0] + delim + 'x0='+ str(x0_or_sigma) + delim + 'cv0=' + str(cv0) 

#         name = str(path / 'plots' / 'solution_plots') + '/' + string

#         plotter = plot(tfinal, M,  N_spaces[count], problem_name, source_name_list[1], rad_or_transport, c, s2, 
#         cv0, x0_or_sigma , mat_or_rad, uncollided, moving, fign, name)

#         fign +=1

#         plotter.plot()

#     for count, tfinal in  enumerate(tfinal_list):
            
#         string = problem_name + delim + str(tfinal) + delim + source_name_list[0] + delim + 'x0='+ str(x0_or_sigma) + delim + 'cv0=' + str(cv0) 

#         name = str(path / 'plots' / 'solution_plots') + '/' + string

#         plotter = plot(tfinal, M,  N_spaces[count + 3], problem_name, source_name_list[0], rad_or_transport, c, s2, 
#         cv0, x0_or_sigma , mat_or_rad, uncollided, moving, fign, name)

#         fign +=1

#         plotter.plot()
    

# def plot_su_olson(M=4, N_space = 128, problem_name = 'su_olson', rad_or_transport = 'rad', 
#                                 c = 0.0, s2 = False, cv0=0.0, x0_or_sigma = 0.5, mat_or_rad ='rad', uncollided = True, moving = True):
#     # tfinal_list = [0.1, 1.0,5.0,10.0,31.6228,100.0]
    
#     source_name_list = ['square_s']
#     fign = 1
#     delim = '_'
#     path = Path("moving_mesh_transport")
#     for tfinal in tfinal_list:
#         for source_name in source_name_list:
            
#             string = problem_name + delim + str(tfinal) + delim + source_name + delim + 'x0='+ str(x0_or_sigma) + delim + 'cv0=' + str(cv0) 

#             name = str(path / 'plots' / 'solution_plots') + '/' + string
#             print(problem_name, 'problem name')
#             plotter = plot(tfinal, M,  N_space, problem_name, source_name, rad_or_transport, c, s2, 
#             cv0, x0_or_sigma , mat_or_rad, uncollided, moving, fign, name)

            

#             plotter.plot()

#             fign +=1
# def plot_su_olson_s2(M=4, N_space = 128, problem_name = 'su_olson_s2', rad_or_transport = 'rad', 
#                                 c = 0.0, s2 = False, cv0=0.0, x0_or_sigma = 0.5, mat_or_rad ='rad', uncollided = True, moving = True):
#     # tfinal_list = [0.1, 1.0,5.0,10.0,31.6228,100.0]
#     tfinal_list = [0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228]
#     source_name_list = ['square_s']
#     fign = 1
#     delim = '_'
#     path = Path("moving_mesh_transport")
#     for tfinal in tfinal_list:
#         for source_name in source_name_list:
            
#             string = problem_name + delim + str(tfinal) + delim + source_name + delim + 'x0='+ str(x0_or_sigma) + delim + 'cv0=' + str(cv0) 

#             name = str(path / 'plots' / 'solution_plots') + '/' + string
#             print(problem_name, 'problem name')
#             plotter = plot(tfinal, M,  N_space, problem_name, source_name, rad_or_transport, c, s2, 
#             cv0, x0_or_sigma , mat_or_rad, uncollided, moving, fign, name)

            

#             plotter.plot()

#             fign +=1
# def plot_su_olson_gaussian(M=12, N_space = 16, problem_name = 'su_olson', rad_or_transport = 'rad', 
#                                 c = 0.0, s2 = False, cv0=0.0, x0_or_sigma = 0.5, mat_or_rad ='rad', uncollided = True, moving = False):
#     # tfinal_list = [0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228, 100.0]
#     tfinal_list = [31.6228]

#     source_name_list = ['gaussian_s']
#     fign = 1
#     delim = '_'
#     path = Path("moving_mesh_transport")
#     for tfinal in tfinal_list:
#         for source_name in source_name_list:
            
#             string = problem_name + delim + str(tfinal) + delim + source_name + delim + 'x0='+ str(x0_or_sigma) + delim + 'cv0=' + str(cv0) 

#             name = str(path / 'plots' / 'solution_plots') + '/' + string
#             print(problem_name, 'problem name')
#             plotter = plot(tfinal, M,  N_space, problem_name, source_name, rad_or_transport, c, s2, 
#             cv0, x0_or_sigma , mat_or_rad, uncollided, moving, fign, name)

            

#             plotter.plot()

#             fign +=1


# def plot_thick_suolson_problems(M=10, N_spaces = [64, 64, 128, 32, 32, 32], problem_name = 'su_olson_thick', rad_or_transport = 'rad', 
#                                 c = 0.0, s2 = False, cv0=0.00, x0_or_sigma = 0.5, mat_or_rad ='rad', uncollided = False, moving = False):
    
#     tfinal_list = [0.3,3.0,30.0]
#     source_name_list = ['square_s']
#     fign = 1
#     delim = '_'
#     path = Path("moving_mesh_transport")
#     for count, tfinal in  enumerate(tfinal_list):
            
#         string = problem_name + delim + str(tfinal) + delim + source_name_list[0] + delim + 'x0='+ str(x0_or_sigma) + delim + 'cv0=' + str(cv0) 

#         name = str(path / 'plots' / 'solution_plots') + '/' + string

#         plotter = plot(tfinal, M,  N_spaces[count], problem_name, source_name_list[0], rad_or_transport, c, s2, 
#         cv0, x0_or_sigma , mat_or_rad, uncollided, moving, fign, name)


#         fign +=1

#         plotter.plot()
#         plotter = plot(tfinal, M,  N_spaces[count], problem_name, source_name_list[0], rad_or_transport, c, s2, 
#         cv0, x0_or_sigma , 'mat', uncollided, moving, fign, name)
#         plotter.plot()

#     # for count, tfinal in  enumerate(tfinal_list):
            
#     #     string = problem_name + delim + str(tfinal) + delim + source_name_list[0] + delim + 'x0='+ str(x0_or_sigma) + delim + 'cv0=' + str(cv0) 

#     #     name = str(path / 'plots' / 'solution_plots') + '/' + string

#     #     plotter = plot(tfinal, M,  N_spaces[count + 3], problem_name, source_name_list[0], rad_or_transport, c, s2, 
#     #     cv0, x0_or_sigma , mat_or_rad, uncollided, moving, fign, name)

#     #     fign +=1

#     #     plotter.plot()

# def make_tables_su_olson(Ms=[10], N_spaces = [32], problem_name = 'su_olson', rad_or_transport = 'rad', 
#                                 c = 0.0, s2 = False, cv0=0.0, x0_or_sigma = 0.5, mat_or_rad ='rad', filenames = ['su_olson_phi.csv','su_olson_e.csv'], source_name_list = ['square_s'], uncollided = True, moving = True):

#     # xs_list = 
#     # source_name_list = ['square_s']\
#     decimals = 6
#     xs_list = su_xs_list

#     tfinal_list = su_tfinal_list
#     fign = 1
#     delim = '_'
#     path = Path("moving_mesh_transport")
#     csv_file = filenames[0]
#     if s2 == True:
#         plus = 3
#     else:
#         plus = 1

#     data_phi = np.zeros((xs_list.size+plus, tfinal_list.size+1))
#     data_e = np.zeros((xs_list.size+plus, tfinal_list.size+1))
#     if s2 == True:
#         data_phi[0,1:] = tfinal_list
#         data_e[0,1:] = tfinal_list
#         data_phi[1:-2,0] = xs_list
#         data_e[1:-2,0] = xs_list
#         # data_phi[-1, 0] = 'RMSE'
#         # data_e[-1, 0] = 'RMSE'
#     else:
#         data_phi[0,1:] = tfinal_list
#         data_e[0,1:] = tfinal_list
#         data_phi[1:,0] = xs_list
#         data_e[1:,0] = xs_list

#     # data_phi = np.zeros((N_space, tfinal_list.size))
#     # data_e = np.zeros((N_space, tfinal_list.size))

#     for count, tfinal in enumerate(tfinal_list):
#         print('t=',tfinal)
#         for source_name in source_name_list:
#             if s2 == True:
#                 if source_name_list[0] == 'square_s': 
#                     benchmark_phi =  load_bench([0,0,0,0,0,0,0,0,1], tfinal, 0.5, 0.0, False)
#                     benchmark_e =  load_bench([0,0,0,0,0,0,0,0,0,1], tfinal, 0.5, 0.0, False)
#                     if tfinal > 10.0:
#                         uncollided = False

#             string = problem_name + delim + str(tfinal) + delim + source_name + delim + 'x0='+ str(x0_or_sigma) + delim + 'cv0=' + str(cv0) 

#             name = str(path / 'plots' / 'solution_plots') + '/' + string
#             print(problem_name, 'problem name')
#             M = Ms[count]
#             N_space = N_spaces[count]


#             xs_quad = quadpy.c1.gauss_legendre(2*M+1).points
#             ws_quad = quadpy.c1.gauss_legendre(2*M+1).weights
#             t_quad = quadpy.c1.gauss_legendre(40).points
#             t_ws = quadpy.c1.gauss_legendre(40).weights


#             l = 1.0

#             plotter = plot(tfinal, M,  N_space, problem_name, source_name, rad_or_transport, c, s2, 
#             cv0, x0_or_sigma , mat_or_rad, uncollided, moving, fign, name)
#             xs, phi = plotter.load()
 

#             quick_build = build_problem.build(N_ang = plotter.N_ang, N_space= N_space, M=M, tfinal = tfinal, x0 = 0.5, 
#             t0 = 10.0, mus = np.array([0.0,0.0 ]), ws = plotter.ws, xs_quad = xs_quad, ws_quad= ws_quad,  sigma_t = 1.0, 
#             sigma_s = 0.0, source_type = np.array([0,0,1,0,0,0,0,0,0,0,0,0,0,0,0]), uncollided = uncollided, moving = moving,
#             move_type = np.array([1]), t_quad = t_quad, t_ws = t_ws, thermal_couple = 1.0, temp_function =  np.array([1,0]), 
#             e_initial = 0.0, sigma =  0, particle_v = 1.0, edge_v = 1, cv0 = 0, thick = False, wave_loc_array = np.zeros((1,1,1)),
#             source_strength = 1.0, move_factor = 1.0, l = 1.0, save_wave_loc= 0, pad =0, leader_pad = 0, quad_thick_source = np.zeros(3),
#             quad_thick_edge = np.zeros(3), boundary_on = np.array([0]), boundary_source_strength = 0, boundary_source = 0, sigma_func = np.array([1]),  Msigma = 0)


#             uncollided_class = uncollided_solution(quick_build)
#             plotter.plot(uncollided_class)

       
            




#             # print(plotter.N_ang, 'angle')
#             # print(plotter.edges, 'edges')


#             output_maker = make_phi.make_output(tfinal, plotter.N_ang, plotter.ws, xs_list, plotter.coeff_mat, M, plotter.edges, uncollided)

#             phi_new = output_maker.make_phi(uncollided_class)
#             e_new = output_maker.make_e()

#             # from scipy.interpolate import interp1d
#             # interp_phi = interp1d(xs, phi_new, kind = 'cubic')
#             # RMSE = np.sqrt(np.mean(phi - interp_phi(np.abs(xs)))**2)
#             # print(RMSE, 'RMSE', '------', tfinal)
#             if s2 == True:
#                 data_phi[1:-2, count+1] = trunc(phi_new, decimals)
#                 data_e[1:-2, count+1] = trunc(e_new, decimals)
                
#             else:
#                 data_phi[1:, count+1] = trunc(phi_new, decimals)
#                 data_e[1:, count+1] = trunc(e_new, decimals)

            
#             if s2 == True and problem_name =='su_olson_s2':
#                 res = [ele for ele in plotter.xs if ele >= 1e-10]
#                 xs_new = np.sort(res) * np.ones(int(plotter.xs.size/2))
#                 xs_new = np.sort(np.append(xs_new, xs_list))
#                 # xs_new = xs_list
       
#                 # xs_new = np.sort(xs_list) * np.ones(xs_list.size)
#                 # xs_points = np.array(xs_list)
#                 xs_points = np.copy(xs_new)
#                 # xs_new = plotter.xs
#                 output_maker = make_phi.make_output(tfinal, plotter.N_ang, plotter.ws, xs_new, plotter.coeff_mat, M, plotter.edges, uncollided)
#                 phi_new = output_maker.make_phi(uncollided_class)
#                 e_new = output_maker.make_e()

#                 output_maker2 = make_phi.make_output(tfinal, plotter.N_ang, plotter.ws, xs_new, plotter.coeff_mat, M, plotter.edges, uncollided)
#                 # if len(xs_new) == len(set(xs_new)):
#                 #     print('no dupes')
#                 # else:
#                 #     print('dupeys!')
#                 phi_new2_intp = interp1d(xs_new, output_maker2.make_phi(uncollided_class), kind = 'cubic')
#                 e_new2_intp = interp1d(xs_new, output_maker2.make_e(), kind = 'cubic')

#                 phi_new2 = xs_points * 0 
#                 e_new2 = xs_points * 0
#                 bench_array = xs_points * 0
#                 for ix in range(xs_points.size):
#                     if xs_points[ix] <= xs_new[-1]:
#                         phi_new2[ix] = phi_new2_intp(xs_points[ix])
#                         e_new2[ix] = e_new2_intp(xs_points[ix])
                        
#                         # if xs_points[ix] == 1.0 and tfinal == 1.0:
#                         #     print(phi_new2[ix], 'phi')
#                         #     print(benchmark_phi(np.array([1.0]))[0], 'bench')
#                         #     print(phi_new2[ix] - benchmark_phi(np.array([1.0]))[0], 'error')

#                 bench_phi_xnew = benchmark_phi(np.abs(xs_new))[0]
#                 bench_phi_xlist = benchmark_phi(np.abs(xs_points))[0]


#                 phi_RMS = np.sqrt(np.mean((phi_new - bench_phi_xnew )**2))
#                 # print(bench_phi_xnew, 'bench phi')
#                 # print(phi_new2, 'phi sol')
#                 phi_RMS_2 = np.sqrt(np.mean((phi_new2 - bench_phi_xlist)**2))
#                 # print(phi_RMS_2, 'phi error here')
#                 # print(phi_new2)
                

#                 e_RMS = np.sqrt(np.mean((e_new - benchmark_e(np.abs(xs_new))[0])**2))
#                 e_RMS_2 = np.sqrt(np.mean((e_new2 - benchmark_e(np.abs(xs_points))[0])**2))
#                 error_e2 = e_new2 - benchmark_e(np.abs(xs_points))[0]
#                 max_error_e2 = max(error_e2)
                
#                 index_of_max_2 = np.argmin(np.abs(e_new2 - benchmark_e(np.abs(xs_points))[0]- max_error_e2))
#                 x_of_max_2 = xs_points[index_of_max_2]

#                 error_e = e_new - benchmark_e(np.abs(xs_new))[0]
#                 max_error_e = max(error_e)
                
#                 index_of_max = np.argmin(np.abs(e_new - benchmark_e(np.abs(xs_new))[0]- max_error_e))
#                 x_of_max = xs_new[index_of_max]
#                 if max_error_e2 > 1e-5:
#                     print('---   ---   ---   ---   ---   ---   ---')
#                     print(' ')
#                     print(' ')
#                     print(' ')
                    
#                     print(x_of_max_2, 'x', max_error_e2, 'error')
#                     print('error e2', error_e2)
#                     print('full')
#                     print(x_of_max, 'x', max_error_e, 'error')
#                     print('error e2', error_e)
#                     print(' ')
#                     print(' ')
#                     print(xs_new, 'xnew')
#                     print(' ')
#                     print('---   ---   ---   ---   ---   ---   ---')

#                 data_e[-2,count+1] = '{:0.3e}'.format(e_RMS)
#                 data_phi[-2, count+1] = '{:0.3e}'.format(phi_RMS)

#                 data_e[-1,count+1] = '{:0.3e}'.format(e_RMS_2)
#                 data_phi[-1, count+1] = '{:0.3e}'.format(phi_RMS_2)
#                 print('--- --- --- --- --- --- --- ---')
#                 print(phi_RMS, 'rmse')
#                 print(e_RMS, 'rmse e')
#                 print('--- --- --- --- --- --- --- ---')

#                 if tfinal == 1.0:
#                     plt.figure(27)
#                     plt.plot(xs_new, e_new, 'bo', mfc = 'none')
#                     plt.plot(xs_points, e_new2, 'rx')
#                     plt.show()



#             # plt.figure(fign)
#             # plt.plot(xs_list, phi_new)
#             # plt.plot(xs, phi, 'o')
#             # plt.show()
#             fign+=1
    
#     with open('tables/' + filenames[0], 'w') as file:
#         writer = csv.writer(file)
#         # data_phi_trunc = trunc(data_phi, 15)
#         writer.writerows(data_phi)

#     with open('tables/' + filenames[1], 'w') as file:
#         writer = csv.writer(file)
#         # data_e_trunc = trunc(data_e, 15)
#         writer.writerows(data_e)

# def make_tables_gaussian_thin(Ms=[10], N_spaces = [32], problem_name = 'su_olson', rad_or_transport = 'rad', 
#                                 c = 0.0, s2 = False, cv0=0.0, x0_or_sigma = 0.5, mat_or_rad ='rad', filenames = ['su_olson_phi.csv','su_olson_e.csv'], source_name_list = ['gaussian_s'], uncollided = True, moving = False):

#     # xs_list = 
#     # source_name_list = ['square_s']
#     xs_list = su_xs_list
#     tfinal_list = su_tfinal_list
#     decimals = 6
#     fign = 1
#     delim = '_'
#     path = Path("moving_mesh_transport")
#     csv_file = filenames[0]
#     if s2 == True:
#         plus = 3
#     else:
#         plus = 1

#     data_phi = np.zeros((xs_list.size+plus, tfinal_list.size+1))
#     data_e = np.zeros((xs_list.size+plus, tfinal_list.size+1))
#     if s2 == True:
#         data_phi[0,1:] = tfinal_list
#         data_e[0,1:] = tfinal_list
#         data_phi[1:-2,0] = xs_list
#         data_e[1:-2,0] = xs_list
#         # data_phi[-1, 0] = 'RMSE'
#         # data_e[-1, 0] = 'RMSE'
#     else:
#         data_phi[0,1:] = tfinal_list
#         data_e[0,1:] = tfinal_list
#         data_phi[1:,0] = xs_list
#         data_e[1:,0] = xs_list


#     # data_phi = np.zeros((N_space, tfinal_list.size))
#     # data_e = np.zeros((N_space, tfinal_list.size))

#     for count, tfinal in enumerate(tfinal_list):
#         print('t=',tfinal)
#         if (tfinal == 31.6228 or tfinal == 100.0) and s2 == True and problem_name == 'su_olson_s2':
#             moving = True
#         for source_name in source_name_list:
#             if s2 == True:
#                 if source_name_list[0] == 'gaussian_s': 
#                     benchmark_phi= load_bench([0,0,0,0,0,0,0,0,0,0,1], tfinal, 0.5, 0.0, False)
#                     benchmark_e =  load_bench([0,0,0,0,0,0,0,0,0,0,0,1], tfinal, 0.5, 0.0, False)
            
#             string = problem_name + delim + str(tfinal) + delim + source_name + delim + 'x0='+ str(x0_or_sigma) + delim + 'cv0=' + str(cv0) 

#             name = str(path / 'plots' / 'solution_plots') + '/' + string
#             print(problem_name, 'problem name')
#             M = Ms[count]
#             N_space = N_spaces[count]
#             plotter = plot(tfinal, M,  N_space, problem_name, source_name, rad_or_transport, c, s2, 
#             cv0, x0_or_sigma , mat_or_rad, uncollided, moving, fign, name)

#             xs_quad = quadpy.c1.gauss_legendre(2*M+1).points
#             ws_quad = quadpy.c1.gauss_legendre(2*M+1).weights
#             t_quad = quadpy.c1.gauss_legendre(40).points
#             t_ws = quadpy.c1.gauss_legendre(40).weights

#             xs, phi = plotter.load()
#             l = 1.0

#             # quick_build = build_problem.build(plotter.N_ang, N_space, M, tfinal, 0.5, 10.0, 1.0, np.array([0.0]), plotter.ws, xs_quad, ws_quad,  np.array([1.0]),  np.array([1.0]), 
#             # np.array([0,0,0,0,0,1,0,0,0,0,0,0,0,0,0]), uncollided, moving,  np.array([1]), t_quad, t_ws, 1.0, np.array([1,0]), 0.0, 0.5, 1.0, 
#             # 1, 0, False, np.zeros((1,1,1)), 1.0, 1.0, 1.0, 0, 0, 0, np.zeros(3), np.zeros(3))


#             quick_build = build_problem.build(N_ang = plotter.N_ang, N_space= N_space, M=M, tfinal = tfinal, x0 = 0.5, 
#             t0 = 10.0, mus = np.array([0.0,0.0 ]), ws = plotter.ws, xs_quad = xs_quad, ws_quad= ws_quad,  sigma_t = 1.0, 
#             sigma_s = 0.0, source_type = np.array([0,0,0,0,0,1,0,0,0,0,0,0,0,0,0]), uncollided = uncollided, moving = moving,
#             move_type = np.array([1]), t_quad = t_quad, t_ws = t_ws, thermal_couple = 1.0, temp_function =  np.array([1,0]), 
#             e_initial = 0.0, sigma =  0.5, particle_v = 1.0, edge_v = 1, cv0 = 0, thick = False, wave_loc_array = np.zeros((1,1,1)),
#             source_strength = 1.0, move_factor = 1.0, l = 1.0, save_wave_loc= 0, pad =0, leader_pad = 0, quad_thick_source = np.zeros(3),
#             quad_thick_edge = np.zeros(3), boundary_on = np.array([0]), boundary_source_strength = 0, boundary_source = 0, sigma_func = np.array([1]),  Msigma = 0)


#             uncollided_class = uncollided_solution(quick_build)
#             plotter.plot(uncollided_class)



#             output_maker = make_phi.make_output(tfinal, plotter.N_ang, plotter.ws, xs_list, plotter.coeff_mat, M, plotter.edges, uncollided)

#             phi_new = output_maker.make_phi(uncollided_class)
#             e_new = output_maker.make_e()

#             # print(plotter.N_ang, 'angle')
#             # print(plotter.edges, 'edges')
#             if s2 == True:
#                 data_phi[1:-2, count+1] = trunc(phi_new, decimals)
#                 data_e[1:-2, count+1] = trunc(e_new, decimals)
                
#             else:
#                 data_phi[1:, count+1] = trunc(phi_new, decimals)
#                 data_e[1:, count+1] = trunc(e_new, decimals)

            
#             if s2 == True and problem_name =='su_olson_s2':
#                 # xs_new = np.linspace(plotter.edges[0], plotter.edges[-1], 100)
#                 xs_new = plotter.xs * np.ones(plotter.xs.size)
#                 xs_points = xs_list
#                 # xs_new = plotter.xs
#                 output_maker = make_phi.make_output(tfinal, plotter.N_ang, plotter.ws, xs_new, plotter.coeff_mat, M, plotter.edges, uncollided)
#                 phi_new = output_maker.make_phi(uncollided_class)
#                 e_new = output_maker.make_e()

#                 output_maker2 = make_phi.make_output(tfinal, plotter.N_ang, plotter.ws, xs_list, plotter.coeff_mat, M, plotter.edges, uncollided)
#                 phi_new2 = interp1d(xs_list, output_maker2.make_phi(uncollided_class), kind = 'cubic')(xs_points)
#                 e_new2 = interp1d( xs_list, output_maker2.make_e(), kind = 'cubic')(xs_points)

#                 phi_RMS = np.sqrt(np.mean((phi_new - benchmark_phi(np.abs(xs_new))[0])**2))
#                 phi_RMS_2 = np.sqrt(np.mean((phi_new2 - benchmark_phi(np.abs(xs_points))[0])**2))

#                 e_RMS = np.sqrt(np.mean((e_new - benchmark_e(np.abs(xs_new))[0])**2))
#                 e_RMS_2 = np.sqrt(np.mean((e_new2 - benchmark_e(np.abs(xs_points))[0])**2))

#                 data_e[-2,count+1] = '{:0.3e}'.format(e_RMS)
#                 data_phi[-2, count+1] = '{:0.3e}'.format(phi_RMS)

#                 data_e[-1,count+1] = '{:0.3e}'.format(e_RMS_2)
#                 data_phi[-1, count+1] = '{:0.3e}'.format(phi_RMS_2)
#                 print('--- --- --- --- --- --- --- ---')
#                 print(phi_RMS, 'rmse')
#                 print(e_RMS, 'rmse e')
#                 print('--- --- --- --- --- --- --- ---')

#             # from scipy.interpolate import interp1d
#             # interp_phi = interp1d(xs, phi_new, kind = 'cubic')
#             # RMSE = np.sqrt(np.mean(phi - interp_phi(np.abs(xs)))**2)
#             # print(RMSE, 'RMSE', '------', tfinal)

#             # plt.plot(xs_list, phi_new)
#             # plt.plot(xs, phi, 'o')
#             # plt.show()
#             fign+=1
    
#     with open('tables/' + filenames[0], 'w') as file:
#         writer = csv.writer(file)
#         # data_phi_trunc = trunc(data_phi, 6)
#         writer.writerows(data_phi)

#     with open('tables/' + filenames[1], 'w') as file:
#         writer = csv.writer(file)
#         # data_e_trunc = trunc(data_e, 6)
#         writer.writerows(data_e)



# def make_tables_thick_problems(Ms=[10], N_spaces = [32], problem_name = 'su_olson', rad_or_transport = 'rad', 
#                                 c = 0.0, s2 = False, cv0=0.0, x0_or_sigma = 0.375, mat_or_rad ='rad', filenames = ['su_olson_phi.csv','su_olson_e.csv'], source_name_list = ['gaussian_s'], uncollided = True, moving = False):

#     # xs_list = 
#     # source_name_list = ['square_s']
#     npnts = 20
#     su_square_edge = 1.1
#     su_gauss_edge = 1.6
#     nl_gauss_edge = 1.4
#     if problem_name in ['su_olson_thick', 'su_olson_thick_s2']:
#         if source_name_list[0] == 'square_s':
#             edge_point = su_square_edge
#         elif source_name_list[0] == 'gaussian_s':
#             edge_point = su_gauss_edge
#     elif problem_name in ['transfer_const_cv=0.03_thick', 'transfer_const_cv=0.03_thick_s2']:
#         edge_point = nl_gauss_edge
    
#     xs_list = np.linspace(0, edge_point, npnts)
#     xs_list = np.round(xs_list, 4)


#     tfinal_list = np.array([0.3, 3.0, 30.0])
#     decimals = 6
#     fign = 1
#     delim = '_'
#     path = Path("moving_mesh_transport")
#     csv_file = filenames[0]


#     data_phi = np.zeros((xs_list.size+1, tfinal_list.size+1))
#     data_e = np.zeros((xs_list.size+1, tfinal_list.size+1))
#     data_phi[0,1:] = tfinal_list
#     data_e[0,1:] = tfinal_list
#     data_phi[1:,0] = xs_list
#     data_e[1:,0] = xs_list


#     # data_phi = np.zeros((N_space, tfinal_list.size))
#     # data_e = np.zeros((N_space, tfinal_list.size))

#     for count, tfinal in enumerate(tfinal_list):
#         print('t=',tfinal)
#         for source_name in source_name_list:
            
#             string = problem_name + delim + str(tfinal) + delim + source_name + delim + 'x0='+ str(x0_or_sigma) + delim + 'cv0=' + str(cv0) 

#             name = str(path / 'plots' / 'solution_plots') + '/' + string
#             print(problem_name, 'problem name')
#             M = Ms[count]
#             N_space = N_spaces[count]
#             plotter = plot(tfinal, M,  N_space, problem_name, source_name, rad_or_transport, c, s2, 
#             cv0, x0_or_sigma , mat_or_rad, uncollided, moving, fign, name)

#             xs_quad = quadpy.c1.gauss_legendre(2*M+1).points
#             ws_quad = quadpy.c1.gauss_legendre(2*M+1).weights
#             t_quad = quadpy.c1.gauss_legendre(40).points
#             t_ws = quadpy.c1.gauss_legendre(40).weights

#             xs, phi = plotter.load()
#             l = 1.0
      
#             # quick_build = build_problem.build(plotter.N_ang, N_space, M, tfinal, 0.5, 10.0, 1.0, np.array([0.0]), plotter.ws, xs_quad, ws_quad,  np.array([1.0]),  np.array([1.0]), 
#             # np.array([0,0,1,0,0,0,0,0,0,0,0,0,0,0,0]), uncollided, moving,  np.array([1]), t_quad, t_ws, 1.0, np.array([1,0]), 0.0, 0.375, 1.0, 
#             # 1, 0, False, np.zeros((1,1,1)), 1.0, 1.0, 1.0, 0, 0, 0, np.zeros(3), np.zeros(3))
#             quick_build = build_problem.build(N_ang = plotter.N_ang, N_space= N_space, M=M, tfinal = tfinal, x0 = 0.5, 
#             t0 = 10.0, mus = np.array([0.0,0.0 ]), ws = plotter.ws, xs_quad = xs_quad, ws_quad= ws_quad,  sigma_t = 1.0, 
#             sigma_s = 0.0, source_type = np.array([0,0,1,0,0,0,0,0,0,0,0,0,0,0,0]), uncollided = uncollided, moving = moving,
#             move_type = np.array([1]), t_quad = t_quad, t_ws = t_ws, thermal_couple = 1.0, temp_function =  np.array([1,0]), 
#             e_initial = 0.0, sigma =  0.5, particle_v = 1.0, edge_v = 1, cv0 = 0, thick = True, wave_loc_array = np.zeros((1,1,1)),
#             source_strength = 1.0, move_factor = 1.0, l = 1.0, save_wave_loc= 0, pad =0, leader_pad = 0, quad_thick_source = np.zeros(3),
#             quad_thick_edge = np.zeros(3), boundary_on = np.array([0]), boundary_source_strength = 0, boundary_source = 0, sigma_func = np.array([1]),  Msigma = 0)


#             uncollided_class = uncollided_solution(quick_build)
#             plotter.plot(uncollided_class)

#             output_maker = make_phi.make_output(tfinal, plotter.N_ang, plotter.ws, xs_list, plotter.coeff_mat, M, plotter.edges, uncollided)

#             phi_new = output_maker.make_phi(uncollided_class)
#             e_new = output_maker.make_e()

#             # print(plotter.N_ang, 'angle')
#             # print(plotter.edges, 'edges')

                

#             data_phi[1:, count+1] = trunc(phi_new, decimals)
#             data_e[1:, count+1] = trunc(e_new, decimals)

            

#             # from scipy.interpolate import interp1d
#             # interp_phi = interp1d(xs, phi_new, kind = 'cubic')
#             # RMSE = np.sqrt(np.mean(phi - interp_phi(np.abs(xs)))**2)
#             # print(RMSE, 'RMSE', '------', tfinal)

#             # plt.plot(xs_list, phi_new)
#             # plt.plot(xs, phi, 'o')
#             # plt.show()
#             fign+=1
    
#     with open('tables/' + filenames[0], 'w') as file:
#         writer = csv.writer(file)
#         # data_phi_trunc = trunc(data_phi, 6)
#         writer.writerows(data_phi)

#     with open('tables/' + filenames[1], 'w') as file:
#         writer = csv.writer(file)
#         # data_e_trunc = trunc(data_e, 6)
#         writer.writerows(data_e)

# # def make_tables_const_cv_thin(M=10, N_space = 32, problem_name = 'rad_transfer_const_cv', rad_or_transport = 'rad', 
# #                                 c = 0.0, s2 = False, cv0=0.03, x0_or_sigma = 0.5, mat_or_rad ='rad', filenames = ['square_const_cv_phi.csv','square_const_cv_e.csv'], source_name_list = ['square_s'], uncollided = True, moving = True):

# #     tfinal_list = np.array([0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228,  100.0])
# #     xs_list = np.array([0.01, 0.1, 0.17783, 0.31623, 0.45, 0.5, 0.56234, 0.75, 1.0, 1.33352, 1.77828, 3.16228, 5.62341, 10.0, 17.78279])
# #     # xs_list = 
# #     # source_name_list = ['square_s']
# #     fign = 1
# #     delim = '_'
# #     path = Path("moving_mesh_transport")
# #     csv_file = 'su_olson.csv'
# #     data_phi = np.zeros((xs_list.size+1, tfinal_list.size+1))
# #     data_e = np.zeros((xs_list.size+1, tfinal_list.size+1))
# #     data_phi[0,1:] = tfinal_list
# #     data_e[0,1:] = tfinal_list
# #     data_phi[1:,0] = xs_list
# #     data_e[1:,0] = xs_list

# #     # data_phi = np.zeros((N_space, tfinal_list.size))
# #     # data_e = np.zeros((N_space, tfinal_list.size))

# #     for count, tfinal in enumerate(tfinal_list):
# #         print('t=',tfinal)
# #         for source_name in source_name_list:
            
# #             string = problem_name + delim + str(tfinal) + delim + source_name + delim + 'x0='+ str(x0_or_sigma) + delim + 'cv0=' + str(cv0) 

# #             name = str(path / 'plots' / 'solution_plots') + '/' + string
# #             print(problem_name, 'problem name')
# #             plotter = plot(tfinal, M,  N_space, problem_name, source_name, rad_or_transport, c, s2, 
# #             cv0, x0_or_sigma , mat_or_rad, uncollided, moving, fign, name)

# #             xs_quad = quadpy.c1.gauss_legendre(2*M+1).points
# #             ws_quad = quadpy.c1.gauss_legendre(2*M+1).weights
# #             t_quad = quadpy.c1.gauss_legendre(40).points
# #             t_ws = quadpy.c1.gauss_legendre(40).weights

# #             xs, phi = plotter.plot()

# #             quick_build = build_problem.build(plotter.N_ang, N_space, M, tfinal, 0.5, 10.0, 1.0, np.array([0.0]), plotter.ws, xs_quad, ws_quad, np.array([1.0]), np.array([0.0]), 
# #             np.array([0,0,1,0,0,0,0,0,0,0,0,0,0,0,0]), uncollided, moving, np.array([1]), t_quad, t_ws, 1.0, np.array([1,0]), 0.0, 0, 1.0, 
# #             1, 0, False, np.zeros((1,1,1)), 1.0, 1.0)

# #             uncollided_class = uncollided_solution(quick_build)



# #             # print(plotter.N_ang, 'angle')
# #             # print(plotter.edges, 'edges')


# #             output_maker = make_phi.make_output(tfinal, plotter.N_ang, plotter.ws, xs_list, plotter.coeff_mat, M, plotter.edges, uncollided)

# #             phi_new = output_maker.make_phi(uncollided_class)
# #             e_new = output_maker.make_e()

# #             # from scipy.interpolate import interp1d
# #             # interp_phi = interp1d(xs, phi_new, kind = 'cubic')
# #             # RMSE = np.sqrt(np.mean(phi - interp_phi(np.abs(xs)))**2)
# #             # print(RMSE, 'RMSE', '------', tfinal)
# #             data_phi[1:, count+1] = phi_new
# #             data_e[1:, count+1] = e_new
# #             # plt.figure(fign)
# #             # plt.plot(xs_list, phi_new)
# #             # plt.plot(xs, phi, 'o')
# #             # plt.show()
# #             fign+=1
    
# #     with open(filenames[0], 'w') as file:
# #         writer = csv.writer(file)
# #         data_phi_trunc = trunc(data_phi, 10)
# #         writer.writerows(data_phi_trunc)

# #     with open(filenames[1], 'w') as file:
# #         writer = csv.writer(file)
# #         data_e_trunc = trunc(data_e, 10)
# #         writer.writerows(data_e_trunc)


# def plot_coeffs_nov23_crc():
#     # plot_coefficients(tfinals = [0.3, 3.0, 30.0],  M=10, source_name = 'gaussian_s',  N_spaces = [128], 
#     # problem_name = 'transfer_const_cv=0.03_thick_s2', rad_or_transport ='transfer', x0_or_sigma = 0.375,
#     # c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = False, s2 = False, moving = False, line = '-',
#     # legend = True, fign = 1)

#     # plt.close()
#     # plt.close()
#     # plt.close()
#     # plt.close()
    
#     # plot_coefficients(tfinals = [0.3, 3.0, 30.0],  M=10, source_name = 'gaussian_s',  N_spaces = [128], 
#     # problem_name = 'su_olson_thick', rad_or_transport ='transfer', x0_or_sigma = 0.375,
#     # c = 0.0, cv0=0.0,mat_or_rad = 'rad', uncollided = False, s2 = False, moving = False, line = '-',
#     # legend = True, fign = 1)


#     # plt.close()
#     # plt.close()
#     # plt.close()
#     # plt.close()

#     # plot_coefficients(tfinals = [0.3, 3.0, 30.0],  M=12, source_name = 'gaussian_s',  N_spaces = [64], 
#     # problem_name = 'su_olson_thick_s2', rad_or_transport ='transfer', x0_or_sigma = 0.375,
#     # c = 0.0, cv0=0.0,mat_or_rad = 'rad', uncollided = False, s2 = False, moving = False, line = '-',
#     # legend = True, fign = 1)


#     # plt.close()
#     # plt.close()
#     # plt.close()
#     # plt.close()
#     # plt.close()

#     # plot_coefficients(tfinals = [0.3, 3.0, 30.0],  M=10, source_name = 'square_s',  N_spaces = [128], 
#     # problem_name = 'su_olson_thick', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     # c = 0.0, cv0=0.0,mat_or_rad = 'rad', uncollided = False, s2 = False, moving = False, line = '-',
#     # legend = True, fign = 1)


#     # plt.close()
#     # plt.close()
#     # plt.close()
#     # plt.close()

#     # plot_coefficients(tfinals = [0.3, 3.0, 30.0],  M=10, source_name = 'square_s',  N_spaces = [128], 
#     # problem_name = 'su_olson_thick_s2', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     # c = 0.0, cv0=0.0,mat_or_rad = 'rad', uncollided = False, s2 = False, moving = False, line = '-',
#     # legend = True, fign = 1)


#     # plt.close()
#     # plt.close()
#     # plt.close()
#     # plt.close()
#     # plt.close()

#     ########### THIN PROBLEMS ##########

#     # Gaussian nonlinear

#     plot_coefficients(tfinals = [0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228, 100.0],  Ms=[12,12,12,12,10,10,10], source_name = 'gaussian_s',  N_spaces = [64,64,64,64,64,64,64], 
#     problem_name = 'transfer_const_cv=0.03', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = False, line = '-',
#     legend = True, fign = 1)


#     plt.close()
#     plt.close()
#     plt.close()
#     plt.close()
    
#     plot_coefficients(tfinals =[0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228, 100.0],  Ms=[12,12,12,12,10,10,10], source_name = 'gaussian_s',  N_spaces = [64,64,64,64,128,128,128], 
#     problem_name = 'transfer_const_cv=0.03_s2', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = False, line = '-',
#     legend = True, fign = 1)


    
#     plt.close()
#     plt.close()
#     plt.close()
#     plt.close()

#     # # SU-OLSON
#     # plot_coefficients(tfinals = [0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228, 100.0],  Ms=[12,12,12,12,12,12,12], source_name = 'gaussian_s',   N_spaces = [64,64,64,64,64,64,64], 
#     # problem_name = 'su_olson', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     # c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = False, line = '-',
#     # legend = True, fign = 1)

#     # plt.close()
#     # plt.close()
#     # plt.close()
#     # plt.close()

#     # plot_coefficients(tfinals = [0.1, 0.31623, 1.0, 3.16228, 10.0],  Ms=[12,12,12,12,12,12,12], source_name = 'gaussian_s',   N_spaces = [64,64,64,64,64,32,8], 
#     # problem_name = 'su_olson_s2', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     # c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = False, line = '-',
#     # legend = True, fign = 1)

#     # plot_coefficients(tfinals = [31.6228, 100.0],  Ms=[8,12], source_name = 'gaussian_s',   N_spaces = [64,32], 
#     # problem_name = 'su_olson_s2', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     # c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = ':',
#     # legend = True, fign = 1)


#     # plt.close()
#     # plt.close()
#     # plt.close()
#     # plt.close()

#     # plot_coefficients(tfinals = [0.1, 0.31623, 1.0, 3.16228],  Ms=[6,6,6,6,6,6,6], source_name = 'square_s',   N_spaces = [128,128,128,128,128,128,128], 
#     # problem_name = 'su_olson_s2', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     # c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
#     # legend = True, fign = 1)

#     # plt.close()
#     # plt.close()
#     # plt.close()
#     # plt.close()
    

# def plot_coeffs_nov28_crc():

#     print('su olson')
#     plot_coefficients(tfinals = [0.1, 0.31623, 1.0, 3.16228],  Ms=[3,4,5,6], source_name = 'square_s',   N_spaces = [256,256,256,256], 
#     problem_name = 'su_olson', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
#     legend = True, fign = 1)

#     plot_coefficients(tfinals = [10.0],  Ms=[6], source_name = 'square_s',   N_spaces = [128], 
#     problem_name = 'su_olson', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
#     legend = True, fign = 1)

#     plot_coefficients(tfinals = [31.6228, 100.0],  Ms=[10,10], source_name = 'square_s',   N_spaces = [32,32], 
#     problem_name = 'su_olson', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
#     legend = True, fign = 1)

#     plt.close()
#     plt.close()
#     plt.close()
#     plt.close()
#     plt.close()
#     plt.close()
#     plt.close()
#     plt.close()
    
#     print('su olson s2')

#     plot_coefficients(tfinals = [0.1, 0.31623, 1.0, 3.16228],  Ms=[4,4,4,4], source_name = 'square_s',   N_spaces = [128,128,128,128], 
#     problem_name = 'su_olson_s2', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
#     legend = True, fign = 1)

#     plot_coefficients(tfinals = [10.0],  Ms=[6], source_name = 'square_s',   N_spaces = [128], 
#     problem_name = 'su_olson_s2', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
#     legend = True, fign = 1)

#     plot_coefficients(tfinals = [31.6228, 100.0],  Ms=[8,8], source_name = 'square_s',   N_spaces = [32,32], 
#     problem_name = 'su_olson_s2', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = False, s2 = False, moving = True, line = '-',
#     legend = True, fign = 1)

#     plt.close()
#     plt.close()
#     plt.close()
#     plt.close()
#     plt.close()
#     plt.close()
#     plt.close()
#     plt.close()

#     print('const cv')
#     plot_coefficients(tfinals = [0.1, 0.31623, 1.0, 3.16228],  Ms=[3,4,5,6], source_name = 'square_s',   N_spaces = [256,256,256,256], 
#     problem_name = 'transfer_const_cv=0.03', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
#     legend = True, fign = 1)

#     plot_coefficients(tfinals = [10.0],  Ms=[8], source_name = 'square_s',   N_spaces = [128], 
#     problem_name = 'transfer_const_cv=0.03', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
#     legend = True, fign = 1)

#     plot_coefficients(tfinals = [31.6228, 100.0],  Ms=[10,8], source_name = 'square_s',   N_spaces = [32,64], 
#     problem_name = 'transfer_const_cv=0.03', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
#     legend = True, fign = 1)

#     plt.close()
#     plt.close()
#     plt.close()
#     plt.close()

#     print('const cv s2')
#     plot_coefficients(tfinals = [0.1, 0.31623, 1.0, 3.16228],  Ms=[4,4,4,4], source_name = 'square_s',   N_spaces = [128,128,128,128], 
#     problem_name = 'transfer_const_cv=0.03_s2', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
#     legend = True, fign = 1)

#     plot_coefficients(tfinals = [10.0],  Ms=[8], source_name = 'square_s',   N_spaces = [128], 
#     problem_name = 'transfer_const_cv=0.03_s2', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
#     legend = True, fign = 1)

#     plot_coefficients(tfinals = [31.6228, 100.0],  Ms=[10,10], source_name = 'square_s',   N_spaces = [32,32], 
#     problem_name = 'transfer_const_cv=0.03_s2', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = False, s2 = False, moving = True, line = '-',
#     legend = True, fign = 1)

#     plt.close()
#     plt.close()
#     plt.close()
#     plt.close()

#     # thick sources

#     print(' --- --- --- --- --- --- ')
#     print('su olson thick')

#     plot_coefficients(tfinals = [0.3, 3.0, 30.0],  Ms=[10,10,10], source_name = 'square_s',  N_spaces = [128,128,128], 
#     problem_name = 'su_olson_thick', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     c = 0.0, cv0=0.0,mat_or_rad = 'rad', uncollided = False, s2 = False, moving = False, line = '-',
#     legend = True, fign = 1)


#     plt.close()
#     plt.close()
#     plt.close()
#     plt.close()

#     print(' --- --- --- --- --- --- ')
#     print('su olson thick s2 square')

#     plot_coefficients(tfinals = [0.3, 3.0, 30.0],  Ms=[10,10,10], source_name = 'square_s',  N_spaces = [128,128,128], 
#     problem_name = 'su_olson_thick_s2', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     c = 0.0, cv0=0.0,mat_or_rad = 'rad', uncollided = False, s2 = False, moving = False, line = '-',
#     legend = True, fign = 1)


#     plt.close()
#     plt.close()
#     plt.close()
#     plt.close()

#     print(' --- --- --- --- --- --- ')
#     print('su olson thick gaussian')
#     plot_coefficients(tfinals = [0.3, 3.0, 30.0],  Ms=[10,10,10], source_name = 'gaussian_s',  N_spaces = [128,128,128], 
#     problem_name = 'su_olson_thick', rad_or_transport ='transfer', x0_or_sigma = 0.375,
#     c = 0.0, cv0=0.0,mat_or_rad = 'rad', uncollided = False, s2 = False, moving = False, line = '-',
#     legend = True, fign = 1)


#     plt.close()
#     plt.close()
#     plt.close()
#     plt.close()


#     print(' --- --- --- --- --- --- ')
#     print('su olson thick s2 gaussian')

#     plot_coefficients(tfinals = [0.3, 3.0, 30.0],  Ms=[10,10,10], source_name = 'gaussian_s',  N_spaces = [128,128,128], 
#     problem_name = 'su_olson_thick_s2', rad_or_transport ='transfer', x0_or_sigma = 0.375,
#     c = 0.0, cv0=0.0,mat_or_rad = 'rad', uncollided = False, s2 = False, moving = False, line = '-',
#     legend = True, fign = 1)


#     plt.close()
#     plt.close()
#     plt.close()
#     plt.close()

#     print(' --- --- --- --- --- --- ')
#     print('cv const thick gaussian')
#     plot_coefficients(tfinals = [0.3, 3.0, 30.0],  Ms=[10,10,10], source_name = 'gaussian_s',  N_spaces = [128,128,128], 
#     problem_name = 'transfer_const_cv=0.03_thick', rad_or_transport ='transfer', x0_or_sigma = 0.375,
#     c = 0.0, cv0=0.0,mat_or_rad = 'rad', uncollided = False, s2 = False, moving = False, line = '-',
#     legend = True, fign = 1)

#     plt.close()
#     plt.close()
#     plt.close()
#     plt.close()


#     print('cv const thick s2 gaussian')
#     plot_coefficients(tfinals = [0.3, 3.0, 30.0],  Ms=[10,10,10], source_name = 'gaussian_s',  N_spaces = [128,128,128], 
#     problem_name = 'transfer_const_cv=0.03_thick_s2', rad_or_transport ='transfer', x0_or_sigma = 0.375,
#     c = 0.0, cv0=0.0,mat_or_rad = 'rad', uncollided = False, s2 = False, moving = False, line = '-',
#     legend = True, fign = 1)

#     plt.close()
#     plt.close()
#     plt.close()
#     plt.close()


#     plot_coefficients(tfinals = [0.3, 3.0],  Ms=[10,10,10], source_name = 'square_s',  N_spaces = [128,128,128], 
#     problem_name = 'su_olson_thick', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     c = 0.0, cv0=0.0,mat_or_rad = 'rad', uncollided = False, s2 = False, moving = False, line = '-',
#     legend = True, fign = 1)


#     plt.close()
#     plt.close()
#     plt.close()
#     plt.close()

#     plot_coefficients(tfinals = [0.3, 3.0, 30.0],  Ms=[10,10,10], source_name = 'gaussian_s',  N_spaces = [128,128,128], 
#     problem_name = 'transfer_const_cv=0.03_thick', rad_or_transport ='transfer', x0_or_sigma = 0.375,
#     c = 0.0, cv0=0.0,mat_or_rad = 'rad', uncollided = False, s2 = False, moving = False, line = '-',
#     legend = True, fign = 1)


#     plt.close()
#     plt.close()
#     plt.close()
#     plt.close()


#     plot_coefficients(tfinals = [0.3, 3.0, 30.0],  Ms=[10,10,10], source_name = 'gaussian_s',  N_spaces = [128,128,128], 
#     problem_name = 'transfer_const_cv=0.03_thick_s2', rad_or_transport ='transfer', x0_or_sigma = 0.375,
#     c = 0.0, cv0=0.0,mat_or_rad = 'rad', uncollided = False, s2 = False, moving = False, line = '-',
#     legend = True, fign = 1)


#     plt.close()
#     plt.close()
#     plt.close()
#     plt.close()




#     # plot_coefficients(tfinals = [0.1, 0.31623, 1.0, 3.16228],  M=3, source_name = 'square_s',  N_spaces = [256], 
#     # problem_name = 'su_olson_s2', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     # c = 0.0, cv0=0.0,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
#     # legend = True, fign = 1)

#     # plt.close()
#     # plt.close()
#     # plt.close()
#     # plt.close()

#     # # NONLINEAR SQUARE

#     # plot_coefficients(tfinals = [0.1, 0.31623, 1.0, 3.16228],  Ms=[6,6,6,6,6,6,6], source_name = 'square_s',   N_spaces = [128,128,128,128,128,128,128], 
#     # problem_name = 'transfer_const_cv=0.03', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     # c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
#     # legend = True, fign = 1)

#     # plt.close()
#     # plt.close()
#     # plt.close()
#     # plt.close()

#     # plt.close()
#     # plt.close()
#     # plt.close()
#     # plt.close()

#     # plot_coefficients(tfinals = [0.1, 0.31623, 1.0,  3.16228],  M=4, source_name = 'square_s',  N_spaces = [128], 
#     # problem_name = 'transfer_const_cv=0.03', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     # c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
#     # legend = True, fign = 1)

#     # plt.close()
#     # plt.close()
#     # plt.close()
#     # plt.close()
# def plot_coeffs_nov31_crc():
     
#     plot_coefficients(tfinals = [10.0],  Ms=[4], source_name = 'square_s',  N_spaces = [128], 
#     problem_name = 'su_olson', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
#     legend = True, fign = 1)

#     plt.close()
#     plt.close()
#     plt.close()

#     plot_coefficients(tfinals = [10.0],  Ms=[8], source_name = 'square_s',   N_spaces = [64], 
#     problem_name = 'transfer_const_cv=0.03', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
#     legend = True, fign = 1)

#     plot_coefficients(tfinals = [31.6228, 100.0],  Ms=[8,8], source_name = 'square_s',   N_spaces = [64,64], 
#     problem_name = 'transfer_const_cv=0.03', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
#     legend = True, fign = 1)

#     plt.close()
#     plt.close()
#     plt.close()
#     plt.close()





# def plot_coeffs_all_local():
#     return 0
    


#     # plot_coefficients(tfinals = [0.3, 3.0, 30.0],  M=10, source_name = 'gaussian_s',  N_spaces = [128], 
#     # problem_name = 'su_olson_thick_s2', rad_or_transport ='transfer', x0_or_sigma = 0.375,
#     # c = 0.0, cv0=0.0,mat_or_rad = 'rad', uncollided = False, s2 = True, moving = False, line = '-',
#     # legend = True, fign = 1)


#     # plt.close()
#     # plt.close()
#     # plt.close()
#     # plt.close()

#     # plot_coefficients(tfinals = [0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228, 100.0],  M=8, source_name = 'square_s',  N_spaces = [32], 
#     # problem_name = 'su_olson', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     # c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
#     # legend = True, fign = 1)

#     # plt.close()
#     # plt.close()
#     # plt.close()
#     # plt.close()
    
#     # plot_coefficients(tfinals = [0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228, 100.0],  M=8, source_name = 'square_s',  N_spaces = [32], 
#     # problem_name = 'su_olson_s2', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     # c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
#     # legend = True, fign = 1)

#     # plt.close()
#     # plt.close()
#     # plt.close()
#     # plt.close()


#     # plot_coefficients(tfinals = [0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228, 100.0],  M=8, source_name = 'gaussian_s',  N_spaces = [16], 
#     # problem_name = 'su_olson', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     # c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
#     # legend = True, fign = 1)

#     # plt.close()
#     # plt.close()
#     # plt.close()
#     # plt.close()
    
#     # plot_coefficients(tfinals = [0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228, 100.0],  M=8, source_name = 'gaussian_s',  N_spaces = [16], 
#     # problem_name = 'su_olson_s2', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     # c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
#     # legend = True, fign = 1)

#     # plt.close()
#     # plt.close()
#     # plt.close()
#     # plt.close()

#     # plot_coefficients(tfinals = [0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228, 100.0],  M=8, source_name = 'square_s',  N_spaces = [32], 
#     # problem_name = 'transfer_const_cv=0.03', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     # c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
#     # legend = True, fign = 1)

#     # plt.close()
#     # plt.close()
#     # plt.close()
#     # plt.close()
    
#     # plot_coefficients(tfinals = [0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228, 100.0],  M=8, source_name = 'square_s',  N_spaces = [32], 
#     # problem_name = 'transfer_const_cv=0.03', rad_or_transport ='transfer', x0_or_sigma = 0.5,
#     # c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
#     # legend = True, fign = 1)

#     # plt.close()
#     # plt.close()
#     # plt.close()
#     # plt.close()

  
 
 

