import numpy as np
import matplotlib.pyplot as plt

from .loading_and_saving.load_solution import load_sol
from .plots.plot_functions.show import show
from .plots.plotting_script import plot_coefficients
from .solver_classes import make_phi, build_problem
from .solver_classes.uncollided_solutions import uncollided_solution
from pathlib import Path
import quadpy
import csv

def trunc(values, decs=0):
    return np.trunc(values*10**decs)/(10**decs)
    
class plot:
    def __init__(self,tfinal, M,  N_space, problem_name, source_name, rad_or_transport, c, s2,
                cv0, x0_or_sigma , mat_or_rad, uncollided, moving, fign, name, mkr1='k-', mkr2='k--', mkr3 = 'k:'):
        self.tfinal = tfinal
        self.M = M
        self.problem_name = problem_name
        self.source_name = source_name
        self.rad_or_transport = rad_or_transport
        self.c = c
        self.s2 = s2
        self.cv0=cv0
        self.x0_or_sigma = x0_or_sigma
        self.mat_or_rad = mat_or_rad
        self.uncollided = uncollided
        self.moving=moving
        self.fign = fign
        self.name = name
        self.N_space = N_space
        self.mkr1 = mkr1
        self.mkr2 = mkr2
        self.mkr3 = mkr3
        

    def plot(self):
        data = load_sol(self.problem_name, self.source_name, self.rad_or_transport, self.c, self.s2, self.cv0)
        data.call_sol(self.tfinal, self.M, self.x0_or_sigma, self.N_space, self.mat_or_rad, self.uncollided, self.moving)
        self.N_ang = np.shape(data.coeff_mat)[0] -1
        self.coeff_mat = data.coeff_mat
        self.ws = data.ws
        self.edges = data.edges
        self.e = data.e

        print('loading solutions for', self.N_space, 'spaces',self.N_ang, 'angles', self.M, 'M')
        
        plt.ion()
        plt.figure(self.fign)
        plt.plot(data.xs, data.phi, self.mkr1)
        plt.plot(data.xs, data.e, self.mkr2)

        if self.problem_name in ['rad_transfer_const_cv', 'rad_transfer_const_cv_s2', 'rad_transfer_const_cv_thick', 'rad_transfer_const_cv_thick_s2', 'transfer_const_cv=0.03_s2']:
            T = data.e / 0.03

        elif self.problem_name in ['su_olson', 'su_olson_s2', 'su_olson_thick', 'su_olson_thick_s2']:
            T = np.power(data.e,.25)
        
        plt.plot(data.xs, T, self.mkr3)

        show(self.name)
        plt.show()
        plt.close()

        return data.xs, data.phi

# plot(1,4,16,'transport', 'square_s', 'transport', 1.0, False, 0.0, 0.5, 'rad', True, True )

def plot_thin_nonlinear_problems(M=12, N_space = 64, problem_name = 'rad_transfer_const_cv', rad_or_transport = 'rad', 
                                c = 0.0, s2 = False, cv0=0.03, x0_or_sigma = 0.5, mat_or_rad ='rad', uncollided = True, moving = False):
    tfinal_list = [0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228,100.0]
    source_name_list = ['gaussian_s']
    fign = 1
    delim = '_'
    path = Path("moving_mesh_transport")
    for tfinal in tfinal_list:
        for source_name in source_name_list:
            
            string = problem_name + delim + str(tfinal) + delim + source_name + delim + 'x0='+ str(x0_or_sigma) + delim + 'cv0=' + str(cv0) 

            name = str(path / 'plots' / 'solution_plots') + '/' + string

            plotter = plot(tfinal, M,  N_space, problem_name, source_name, rad_or_transport, c, s2, 
            cv0, x0_or_sigma , mat_or_rad, uncollided, moving, fign, name)

            fign +=1

        plotter.plot()

def plot_thin_nonlinear_problems_s2(M=10, N_space = 128, problem_name = 'transfer_const_cv=0.03_s2', rad_or_transport = 'rad', 
                                c = 0.0, s2 = False, cv0=0.03, x0_or_sigma = 0.5, mat_or_rad ='rad', uncollided = True, moving = False):
    # tfinal_list = [0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228,100.0]
    tfinal_list = [10.0, 31.6228,100.0]

    source_name_list = ['gaussian_s']
    fign = 1
    delim = '_'
    path = Path("moving_mesh_transport")
    for tfinal in tfinal_list:
        for source_name in source_name_list:
            
            string = problem_name + delim + str(tfinal) + delim + source_name + delim + 'x0='+ str(x0_or_sigma) + delim + 'cv0=' + str(cv0) 

            name = str(path / 'plots' / 'solution_plots') + '/' + string

            plotter = plot(tfinal, M,  N_space, problem_name, source_name, rad_or_transport, c, s2, 
            cv0, x0_or_sigma , mat_or_rad, uncollided, moving, fign, name)

            fign +=1

        plotter.plot()


def plot_thick_nonlinear_problems(M=10, N_spaces = [64, 64, 128], problem_name = 'rad_transfer_const_cv_thick', rad_or_transport = 'rad', 
                                c = 0.0, s2 = False, cv0=0.03, x0_or_sigma = 0.375, mat_or_rad ='rad', uncollided = False, moving = False):
    
    tfinal_list = [0.3, 3.0, 30.0]
    source_name_list = ['gaussian_s']
    fign = 1
    delim = '_'
    path = Path("moving_mesh_transport")
    for count, tfinal in  enumerate(tfinal_list):
            
        string = problem_name + delim + str(tfinal) + delim + source_name_list[0] + delim + 'x0='+ str(x0_or_sigma) + delim + 'cv0=' + str(cv0) 

        name = str(path / 'plots' / 'solution_plots') + '/' + string

        plotter = plot(tfinal, M,  N_spaces[count], problem_name, source_name_list[0], rad_or_transport, c, s2, 
        cv0, x0_or_sigma , mat_or_rad, uncollided, moving, fign, name)

        fign +=1

        plotter.plot()

  


def plot_thick_nonlinear_problems_s2(M=10, N_spaces = [64, 64, 128], problem_name = 'rad_transfer_const_cv_thick_s2', rad_or_transport = 'rad', 
                                c = 0.0, s2 = False, cv0=0.03, x0_or_sigma = 400.0, mat_or_rad ='rad', uncollided = False, moving = False):
    
    tfinal_list = [0.3, 3.0, 3.0]
    source_name_list = ['gaussian_s']
    fign = 1
    delim = '_'
    path = Path("moving_mesh_transport")
    for count, tfinal in  enumerate(tfinal_list):
            
        string = problem_name + delim + str(tfinal) + delim + source_name_list[0] + delim + 'x0='+ str(x0_or_sigma) + delim + 'cv0=' + str(cv0) 

        name = str(path / 'plots' / 'solution_plots') + '/' + string

        plotter = plot(tfinal, M,  N_spaces[count], problem_name, source_name_list[1], rad_or_transport, c, s2, 
        cv0, x0_or_sigma , mat_or_rad, uncollided, moving, fign, name)

        fign +=1

        plotter.plot()

    for count, tfinal in  enumerate(tfinal_list):
            
        string = problem_name + delim + str(tfinal) + delim + source_name_list[0] + delim + 'x0='+ str(x0_or_sigma) + delim + 'cv0=' + str(cv0) 

        name = str(path / 'plots' / 'solution_plots') + '/' + string

        plotter = plot(tfinal, M,  N_spaces[count + 3], problem_name, source_name_list[0], rad_or_transport, c, s2, 
        cv0, x0_or_sigma , mat_or_rad, uncollided, moving, fign, name)

        fign +=1

        plotter.plot()
    

def plot_su_olson(M=4, N_space = 128, problem_name = 'su_olson', rad_or_transport = 'rad', 
                                c = 0.0, s2 = False, cv0=0.0, x0_or_sigma = 0.5, mat_or_rad ='rad', uncollided = True, moving = True):
    # tfinal_list = [0.1, 1.0,5.0,10.0,31.6228,100.0]
    tfinal_list = [0.1, 0.31623, 1.0]
    source_name_list = ['square_s']
    fign = 1
    delim = '_'
    path = Path("moving_mesh_transport")
    for tfinal in tfinal_list:
        for source_name in source_name_list:
            
            string = problem_name + delim + str(tfinal) + delim + source_name + delim + 'x0='+ str(x0_or_sigma) + delim + 'cv0=' + str(cv0) 

            name = str(path / 'plots' / 'solution_plots') + '/' + string
            print(problem_name, 'problem name')
            plotter = plot(tfinal, M,  N_space, problem_name, source_name, rad_or_transport, c, s2, 
            cv0, x0_or_sigma , mat_or_rad, uncollided, moving, fign, name)

            

            plotter.plot()

            fign +=1
def plot_su_olson_s2(M=4, N_space = 128, problem_name = 'su_olson_s2', rad_or_transport = 'rad', 
                                c = 0.0, s2 = False, cv0=0.0, x0_or_sigma = 0.5, mat_or_rad ='rad', uncollided = True, moving = True):
    # tfinal_list = [0.1, 1.0,5.0,10.0,31.6228,100.0]
    tfinal_list = [0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228]
    source_name_list = ['square_s']
    fign = 1
    delim = '_'
    path = Path("moving_mesh_transport")
    for tfinal in tfinal_list:
        for source_name in source_name_list:
            
            string = problem_name + delim + str(tfinal) + delim + source_name + delim + 'x0='+ str(x0_or_sigma) + delim + 'cv0=' + str(cv0) 

            name = str(path / 'plots' / 'solution_plots') + '/' + string
            print(problem_name, 'problem name')
            plotter = plot(tfinal, M,  N_space, problem_name, source_name, rad_or_transport, c, s2, 
            cv0, x0_or_sigma , mat_or_rad, uncollided, moving, fign, name)

            

            plotter.plot()

            fign +=1
def plot_su_olson_gaussian(M=12, N_space = 64, problem_name = 'su_olson', rad_or_transport = 'rad', 
                                c = 0.0, s2 = False, cv0=0.0, x0_or_sigma = 0.5, mat_or_rad ='rad', uncollided = True, moving = False):
    tfinal_list = [0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228, 100.0]
    source_name_list = ['gaussian_s']
    fign = 1
    delim = '_'
    path = Path("moving_mesh_transport")
    for tfinal in tfinal_list:
        for source_name in source_name_list:
            
            string = problem_name + delim + str(tfinal) + delim + source_name + delim + 'x0='+ str(x0_or_sigma) + delim + 'cv0=' + str(cv0) 

            name = str(path / 'plots' / 'solution_plots') + '/' + string
            print(problem_name, 'problem name')
            plotter = plot(tfinal, M,  N_space, problem_name, source_name, rad_or_transport, c, s2, 
            cv0, x0_or_sigma , mat_or_rad, uncollided, moving, fign, name)

            

            plotter.plot()

            fign +=1


def plot_thick_suolson_problems(M=10, N_spaces = [64, 64, 128, 32, 32, 32], problem_name = 'su_olson_thick', rad_or_transport = 'rad', 
                                c = 0.0, s2 = False, cv0=0.00, x0_or_sigma = 0.5, mat_or_rad ='rad', uncollided = False, moving = False):
    
    tfinal_list = [0.3,3.0,30.0]
    source_name_list = ['square_s']
    fign = 1
    delim = '_'
    path = Path("moving_mesh_transport")
    for count, tfinal in  enumerate(tfinal_list):
            
        string = problem_name + delim + str(tfinal) + delim + source_name_list[0] + delim + 'x0='+ str(x0_or_sigma) + delim + 'cv0=' + str(cv0) 

        name = str(path / 'plots' / 'solution_plots') + '/' + string

        plotter = plot(tfinal, M,  N_spaces[count], problem_name, source_name_list[0], rad_or_transport, c, s2, 
        cv0, x0_or_sigma , mat_or_rad, uncollided, moving, fign, name)


        fign +=1

        plotter.plot()
        plotter = plot(tfinal, M,  N_spaces[count], problem_name, source_name_list[0], rad_or_transport, c, s2, 
        cv0, x0_or_sigma , 'mat', uncollided, moving, fign, name)
        plotter.plot()

    # for count, tfinal in  enumerate(tfinal_list):
            
    #     string = problem_name + delim + str(tfinal) + delim + source_name_list[0] + delim + 'x0='+ str(x0_or_sigma) + delim + 'cv0=' + str(cv0) 

    #     name = str(path / 'plots' / 'solution_plots') + '/' + string

    #     plotter = plot(tfinal, M,  N_spaces[count + 3], problem_name, source_name_list[0], rad_or_transport, c, s2, 
    #     cv0, x0_or_sigma , mat_or_rad, uncollided, moving, fign, name)

    #     fign +=1

    #     plotter.plot()

def make_tables_su_olson(M=10, N_space = 32, problem_name = 'su_olson', rad_or_transport = 'rad', 
                                c = 0.0, s2 = False, cv0=0.0, x0_or_sigma = 0.5, mat_or_rad ='rad', filenames = ['su_olson_phi.csv','su_olson_e.csv'], source_name_list = ['square_s'], uncollided = True, moving = True):
    tfinal_list = np.array([0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228,  100.0])
    xs_list = np.array([0.01, 0.1, 0.17783, 0.31623, 0.45, 0.5, 0.56234, 0.75, 1.0, 1.33352, 1.77828, 3.16228, 5.62341, 10.0, 17.78279])
    # xs_list = 
    # source_name_list = ['square_s']
    fign = 1
    delim = '_'
    path = Path("moving_mesh_transport")
    csv_file = 'su_olson.csv'
    data_phi = np.zeros((xs_list.size+1, tfinal_list.size+1))
    data_e = np.zeros((xs_list.size+1, tfinal_list.size+1))
    data_phi[0,1:] = tfinal_list
    data_e[0,1:] = tfinal_list
    data_phi[1:,0] = xs_list
    data_e[1:,0] = xs_list

    # data_phi = np.zeros((N_space, tfinal_list.size))
    # data_e = np.zeros((N_space, tfinal_list.size))

    for count, tfinal in enumerate(tfinal_list):
        print('t=',tfinal)
        for source_name in source_name_list:
            
            string = problem_name + delim + str(tfinal) + delim + source_name + delim + 'x0='+ str(x0_or_sigma) + delim + 'cv0=' + str(cv0) 

            name = str(path / 'plots' / 'solution_plots') + '/' + string
            print(problem_name, 'problem name')
            plotter = plot(tfinal, M,  N_space, problem_name, source_name, rad_or_transport, c, s2, 
            cv0, x0_or_sigma , mat_or_rad, uncollided, moving, fign, name)

            xs_quad = quadpy.c1.gauss_legendre(2*M+1).points
            ws_quad = quadpy.c1.gauss_legendre(2*M+1).weights
            t_quad = quadpy.c1.gauss_legendre(40).points
            t_ws = quadpy.c1.gauss_legendre(40).weights

            xs, phi = plotter.plot()

            quick_build = build_problem.build(plotter.N_ang, N_space, M, tfinal, 0.5, 10.0, 1.0, np.array([0.0]), plotter.ws, xs_quad, ws_quad, np.array([1.0]), np.array([0.0]), 
            np.array([0,0,1,0,0,0,0,0,0,0,0,0,0,0,0]), uncollided, moving, np.array([1]), t_quad, t_ws, 1.0, np.array([1,0]), 0.0, 0, 1.0, 
            1, 0, False, np.zeros((1,1,1)), 1.0, 1.0)

            uncollided_class = uncollided_solution(quick_build)



            print(plotter.N_ang, 'angle')
            print(plotter.edges, 'edges')


            output_maker = make_phi.make_output(tfinal, plotter.N_ang, plotter.ws, xs_list, plotter.coeff_mat, M, plotter.edges, uncollided)

            phi_new = output_maker.make_phi(uncollided_class)
            e_new = output_maker.make_e()

            # from scipy.interpolate import interp1d
            # interp_phi = interp1d(xs, phi_new, kind = 'cubic')
            # RMSE = np.sqrt(np.mean(phi - interp_phi(np.abs(xs)))**2)
            # print(RMSE, 'RMSE', '------', tfinal)
            data_phi[1:, count+1] = phi_new
            data_e[1:, count+1] = e_new
            # plt.figure(fign)
            # plt.plot(xs_list, phi_new)
            # plt.plot(xs, phi, 'o')
            # plt.show()
            fign+=1
    
    with open(filenames[0], 'w') as file:
        writer = csv.writer(file)
        data_phi_trunc = trunc(data_phi, 10)
        writer.writerows(data_phi_trunc)

    with open(filenames[1], 'w') as file:
        writer = csv.writer(file)
        data_e_trunc = trunc(data_e, 10)
        writer.writerows(data_e_trunc)

def make_tables_const_cv_thin(M=10, N_space = 32, problem_name = 'rad_transfer_const_cv', rad_or_transport = 'rad', 
                                c = 0.0, s2 = False, cv0=0.03, x0_or_sigma = 0.5, mat_or_rad ='rad', filenames = ['square_const_cv_phi.csv','square_const_cv_e.csv'], source_name_list = ['square_s'], uncollided = True, moving = True):

    tfinal_list = np.array([0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228,  100.0])
    xs_list = np.array([0.01, 0.1, 0.17783, 0.31623, 0.45, 0.5, 0.56234, 0.75, 1.0, 1.33352, 1.77828, 3.16228, 5.62341, 10.0, 17.78279])
    # xs_list = 
    # source_name_list = ['square_s']
    fign = 1
    delim = '_'
    path = Path("moving_mesh_transport")
    csv_file = 'su_olson.csv'
    data_phi = np.zeros((xs_list.size+1, tfinal_list.size+1))
    data_e = np.zeros((xs_list.size+1, tfinal_list.size+1))
    data_phi[0,1:] = tfinal_list
    data_e[0,1:] = tfinal_list
    data_phi[1:,0] = xs_list
    data_e[1:,0] = xs_list

    # data_phi = np.zeros((N_space, tfinal_list.size))
    # data_e = np.zeros((N_space, tfinal_list.size))

    for count, tfinal in enumerate(tfinal_list):
        print('t=',tfinal)
        for source_name in source_name_list:
            
            string = problem_name + delim + str(tfinal) + delim + source_name + delim + 'x0='+ str(x0_or_sigma) + delim + 'cv0=' + str(cv0) 

            name = str(path / 'plots' / 'solution_plots') + '/' + string
            print(problem_name, 'problem name')
            plotter = plot(tfinal, M,  N_space, problem_name, source_name, rad_or_transport, c, s2, 
            cv0, x0_or_sigma , mat_or_rad, uncollided, moving, fign, name)

            xs_quad = quadpy.c1.gauss_legendre(2*M+1).points
            ws_quad = quadpy.c1.gauss_legendre(2*M+1).weights
            t_quad = quadpy.c1.gauss_legendre(40).points
            t_ws = quadpy.c1.gauss_legendre(40).weights

            xs, phi = plotter.plot()

            quick_build = build_problem.build(plotter.N_ang, N_space, M, tfinal, 0.5, 10.0, 1.0, np.array([0.0]), plotter.ws, xs_quad, ws_quad, np.array([1.0]), np.array([0.0]), 
            np.array([0,0,1,0,0,0,0,0,0,0,0,0,0,0,0]), uncollided, moving, np.array([1]), t_quad, t_ws, 1.0, np.array([1,0]), 0.0, 0, 1.0, 
            1, 0, False, np.zeros((1,1,1)), 1.0, 1.0)

            uncollided_class = uncollided_solution(quick_build)



            print(plotter.N_ang, 'angle')
            print(plotter.edges, 'edges')


            output_maker = make_phi.make_output(tfinal, plotter.N_ang, plotter.ws, xs_list, plotter.coeff_mat, M, plotter.edges, uncollided)

            phi_new = output_maker.make_phi(uncollided_class)
            e_new = output_maker.make_e()

            # from scipy.interpolate import interp1d
            # interp_phi = interp1d(xs, phi_new, kind = 'cubic')
            # RMSE = np.sqrt(np.mean(phi - interp_phi(np.abs(xs)))**2)
            # print(RMSE, 'RMSE', '------', tfinal)
            data_phi[1:, count+1] = phi_new
            data_e[1:, count+1] = e_new
            # plt.figure(fign)
            # plt.plot(xs_list, phi_new)
            # plt.plot(xs, phi, 'o')
            # plt.show()
            fign+=1
    
    with open(filenames[0], 'w') as file:
        writer = csv.writer(file)
        data_phi_trunc = trunc(data_phi, 10)
        writer.writerows(data_phi_trunc)

    with open(filenames[1], 'w') as file:
        writer = csv.writer(file)
        data_e_trunc = trunc(data_e, 10)
        writer.writerows(data_e_trunc)


def plot_coeffs_all_crc():
    # plot_coefficients(tfinals = [0.3, 3.0, 30.0],  M=10, source_name = 'gaussian_s',  N_spaces = [128], 
    # problem_name = 'transfer_const_cv=0.03_thick_s2', rad_or_transport ='transfer', x0_or_sigma = 0.375,
    # c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = False, s2 = False, moving = False, line = '-',
    # legend = True, fign = 1)

    # plt.close()
    # plt.close()
    # plt.close()
    # plt.close()
    
    # plot_coefficients(tfinals = [0.3, 3.0, 30.0],  M=10, source_name = 'gaussian_s',  N_spaces = [128], 
    # problem_name = 'su_olson_thick', rad_or_transport ='transfer', x0_or_sigma = 0.375,
    # c = 0.0, cv0=0.0,mat_or_rad = 'rad', uncollided = False, s2 = False, moving = False, line = '-',
    # legend = True, fign = 1)


    # plt.close()
    # plt.close()
    # plt.close()
    # plt.close()

    # plot_coefficients(tfinals = [0.3, 3.0, 30.0],  M=12, source_name = 'gaussian_s',  N_spaces = [64], 
    # problem_name = 'su_olson_thick_s2', rad_or_transport ='transfer', x0_or_sigma = 0.375,
    # c = 0.0, cv0=0.0,mat_or_rad = 'rad', uncollided = False, s2 = False, moving = False, line = '-',
    # legend = True, fign = 1)


    # plt.close()
    # plt.close()
    # plt.close()
    # plt.close()
    # plt.close()

    # plot_coefficients(tfinals = [0.3, 3.0, 30.0],  M=10, source_name = 'square_s',  N_spaces = [128], 
    # problem_name = 'su_olson_thick', rad_or_transport ='transfer', x0_or_sigma = 0.5,
    # c = 0.0, cv0=0.0,mat_or_rad = 'rad', uncollided = False, s2 = False, moving = False, line = '-',
    # legend = True, fign = 1)


    # plt.close()
    # plt.close()
    # plt.close()
    # plt.close()

    # plot_coefficients(tfinals = [0.3, 3.0, 30.0],  M=10, source_name = 'square_s',  N_spaces = [128], 
    # problem_name = 'su_olson_thick_s2', rad_or_transport ='transfer', x0_or_sigma = 0.5,
    # c = 0.0, cv0=0.0,mat_or_rad = 'rad', uncollided = False, s2 = False, moving = False, line = '-',
    # legend = True, fign = 1)


    # plt.close()
    # plt.close()
    # plt.close()
    # plt.close()
    # plt.close()

    ########### THIN PROBLEMS ##########

    # Gaussian nonlinear

    plot_coefficients(tfinals = [0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228, 100.0],  Ms=[12,12,12,12,10,10,10], source_name = 'gaussian_s',  N_spaces = [64,64,64,64,128,128,128], 
    problem_name = 'transfer_const_cv=0.03', rad_or_transport ='transfer', x0_or_sigma = 0.5,
    c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = False, line = '-',
    legend = True, fign = 1)


    plt.close()
    plt.close()
    plt.close()
    plt.close()
    
    plot_coefficients(tfinals =[0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228, 100.0],  Ms=[12,12,12,12,10,10,10], source_name = 'gaussian_s',  N_spaces = [64,64,64,64,128,128,128], 
    problem_name = 'transfer_const_cv=0.03_s2', rad_or_transport ='transfer', x0_or_sigma = 0.5,
    c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = False, line = '-',
    legend = True, fign = 1)


    
    plt.close()
    plt.close()
    plt.close()
    plt.close()

    # SU-OLSON

    plot_coefficients(tfinals = [0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228, 100.0],  Ms=[12,12,12,12,12,12,12], source_name = 'gaussian_s',   N_spaces = [64,64,64,64,64,64,64], 
    problem_name = 'su_olson', rad_or_transport ='transfer', x0_or_sigma = 0.5,
    c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = False, line = '-',
    legend = True, fign = 1)

    plt.close()
    plt.close()
    plt.close()
    plt.close()

    plot_coefficients(tfinals = [0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228, 100.0],  Ms=[12,12,12,12,12,12,12], source_name = 'gaussian_s',   N_spaces = [64,64,64,64,64,64,64], 
    problem_name = 'su_olson_s2', rad_or_transport ='transfer', x0_or_sigma = 0.5,
    c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = False, line = '-',
    legend = True, fign = 1)

    # plot_coefficients(tfinals = [0.1, 0.31623, 1.0, 3.16228],  M=4, source_name = 'square_s',  N_spaces = [128], 
    # problem_name = 'su_olson', rad_or_transport ='transfer', x0_or_sigma = 0.5,
    # c = 0.0, cv0=0.0,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
    # legend = True, fign = 1)

    # plt.close()
    # plt.close()
    # plt.close()
    # plt.close()

    # plot_coefficients(tfinals = [0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228],  M=4, source_name = 'square_s',  N_spaces = [128], 
    # problem_name = 'su_olson_s2', rad_or_transport ='transfer', x0_or_sigma = 0.5,
    # c = 0.0, cv0=0.0,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
    # legend = True, fign = 1)

    # plt.close()
    # plt.close()
    # plt.close()
    # plt.close()

    # # NONLINEAR SQUARE

    # plot_coefficients(tfinals = [0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228],  M=4, source_name = 'square_s',  N_spaces = [128], 
    # problem_name = 'transfer_const_cv=0.03_s2', rad_or_transport ='transfer', x0_or_sigma = 0.5,
    # c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
    # legend = True, fign = 1)

    # plt.close()
    # plt.close()
    # plt.close()
    # plt.close()

    # plot_coefficients(tfinals = [0.1, 0.31623, 1.0,  3.16228],  M=4, source_name = 'square_s',  N_spaces = [128], 
    # problem_name = 'transfer_const_cv=0.03', rad_or_transport ='transfer', x0_or_sigma = 0.5,
    # c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
    # legend = True, fign = 1)

    # plt.close()
    # plt.close()
    # plt.close()
    # plt.close()




def plot_coeffs_all_local():
    plot_coefficients(tfinals = [0.3, 3.0, 30.0],  M=8, source_name = 'gaussian_s',  N_spaces = [64], 
    problem_name = 'transfer_const_cv=0.03_thick', rad_or_transport ='transfer', x0_or_sigma = 0.375,
    c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = False, s2 = False, moving = False, line = '-',
    legend = True, fign = 1)

    plt.close()
    plt.close()
    plt.close()
    plt.close()
    
    plot_coefficients(tfinals = [0.3, 3.0, 30.0],  M=10, source_name = 'gaussian_s',  N_spaces = [128], 
    problem_name = 'su_olson_thick', rad_or_transport ='transfer', x0_or_sigma = 0.375,
    c = 0.0, cv0=0.0,mat_or_rad = 'rad', uncollided = False, s2 = False, moving = False, line = '-',
    legend = True, fign = 1)


    plt.close()
    plt.close()
    plt.close()
    plt.close()

    plot_coefficients(tfinals = [0.3, 3.0, 30.0],  M=10, source_name = 'gaussian_s',  N_spaces = [128], 
    problem_name = 'su_olson_thick_s2', rad_or_transport ='transfer', x0_or_sigma = 0.375,
    c = 0.0, cv0=0.0,mat_or_rad = 'rad', uncollided = False, s2 = True, moving = False, line = '-',
    legend = True, fign = 1)


    plt.close()
    plt.close()
    plt.close()
    plt.close()

    plot_coefficients(tfinals = [0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228, 100.0],  M=8, source_name = 'square_s',  N_spaces = [32], 
    problem_name = 'su_olson', rad_or_transport ='transfer', x0_or_sigma = 0.5,
    c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
    legend = True, fign = 1)

    plt.close()
    plt.close()
    plt.close()
    plt.close()
    
    plot_coefficients(tfinals = [0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228, 100.0],  M=8, source_name = 'square_s',  N_spaces = [32], 
    problem_name = 'su_olson_s2', rad_or_transport ='transfer', x0_or_sigma = 0.5,
    c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
    legend = True, fign = 1)

    plt.close()
    plt.close()
    plt.close()
    plt.close()


    plot_coefficients(tfinals = [0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228, 100.0],  M=8, source_name = 'gaussian_s',  N_spaces = [16], 
    problem_name = 'su_olson', rad_or_transport ='transfer', x0_or_sigma = 0.5,
    c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
    legend = True, fign = 1)

    plt.close()
    plt.close()
    plt.close()
    plt.close()
    
    plot_coefficients(tfinals = [0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228, 100.0],  M=8, source_name = 'gaussian_s',  N_spaces = [16], 
    problem_name = 'su_olson_s2', rad_or_transport ='transfer', x0_or_sigma = 0.5,
    c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
    legend = True, fign = 1)

    plt.close()
    plt.close()
    plt.close()
    plt.close()

    plot_coefficients(tfinals = [0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228, 100.0],  M=8, source_name = 'square_s',  N_spaces = [32], 
    problem_name = 'transfer_const_cv=0.03', rad_or_transport ='transfer', x0_or_sigma = 0.5,
    c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
    legend = True, fign = 1)

    plt.close()
    plt.close()
    plt.close()
    plt.close()
    
    plot_coefficients(tfinals = [0.1, 0.31623, 1.0, 3.16228, 10.0, 31.6228, 100.0],  M=8, source_name = 'square_s',  N_spaces = [32], 
    problem_name = 'transfer_const_cv=0.03', rad_or_transport ='transfer', x0_or_sigma = 0.5,
    c = 0.0, cv0=0.03,mat_or_rad = 'rad', uncollided = True, s2 = False, moving = True, line = '-',
    legend = True, fign = 1)

    plt.close()
    plt.close()
    plt.close()
    plt.close()

  
 
 

