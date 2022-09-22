import numpy as np
import matplotlib.pyplot as plt

from .loading_and_saving.load_solution import load_sol
from .plots.plot_functions.show import show
from .solver_classes import make_phi, build_problem, uncollided_solutions
from pathlib import Path
import quadpy
import csv
class plot:
    def __init__(self,tfinal, M,  N_space, problem_name, source_name, rad_or_transport, c, s2,
                cv0, x0_or_sigma , mat_or_rad, uncollided, moving, fign, name):
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

    def plot(self):
        data = load_sol(self.problem_name, self.source_name, self.rad_or_transport, self.c, self.s2, self.cv0)
        data.call_sol(self.tfinal, self.M, self.x0_or_sigma, self.N_space, self.mat_or_rad, self.uncollided, self.moving)
        self.N_ang = np.shape(data.coeff_mat)[0] -1
        self.coeff_mat = data.coeff_mat
        self.ws = data.ws
        self.edges = data.edges

        print('loading solutions for', self.N_space, 'spaces',self.N_ang, 'angles', self.M, 'M')
        
        plt.ion()
        plt.figure(self.fign)
        plt.plot(data.xs, data.phi, 'k-')
        plt.plot(data.xs, data.e, 'k--')

        show(self.name)
        plt.show()
        plt.close()

        return data.xs, data.phi

# plot(1,4,16,'transport', 'square_s', 'transport', 1.0, False, 0.0, 0.5, 'rad', True, True )

def plot_thin_nonlinear_problems(M=0, N_space = 16, problem_name = 'rad_transfer_const_cv', rad_or_transport = 'rad', 
                                c = 0.0, s2 = False, cv0=0.03, x0_or_sigma = 0.5, mat_or_rad ='rad', uncollided = True, moving = True):
    tfinal_list = [1.0,10.0,31.6228]
    source_name_list = ['gaussian_s', 'square_s']
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

def plot_su_olson(M=6, N_space = 16, problem_name = 'su_olson', rad_or_transport = 'rad', 
                                c = 0.0, s2 = False, cv0=0.0, x0_or_sigma = 0.5, mat_or_rad ='rad', uncollided = True, moving = True):
    tfinal_list = [0.1, 1.0,5.0,10.0,31.6228]
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


def make_tables_su_olson(M=0, N_space = 8, problem_name = 'su_olson', rad_or_transport = 'rad', 
                                c = 0.0, s2 = False, cv0=0.0, x0_or_sigma = 0.5, mat_or_rad ='rad', uncollided = True, moving = True):
    tfinal_list = np.array([0.1, 1.0, 10.0, 31.6228])
    xs_list = np.array([0.01, 0.1, 0.17783, 0.31623, 0.45, 0.5, 0.56234, 0.75, 1.0, 1.33352, 1.77828, 3.16228, 5.62341, 10.0, 17.78279])
    source_name_list = ['square_s']
    fign = 1
    delim = '_'
    path = Path("moving_mesh_transport")
    csv_file = 'su_olson.csv'
    data_phi = np.zeros((xs_list.size, tfinal_list.size))
    data_e = np.zeros((xs_list.size, tfinal_list.size))

    for count, tfinal in enumerate(tfinal_list):
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

            quick_build = build_problem.build(plotter.N_ang, N_space, M, tfinal, 0.5, 10, 1.0, np.array([0.0]), plotter.ws, xs_quad, ws_quad, np.array([1.0]), np.array([0.0]), 
            np.array([0,0,1,0,0,0,0,0,0,0,0,0,0,0,0]), uncollided, moving, np.array([1]), t_quad, t_ws, 1.0, np.array([1,0]), 0.0, 0, 1.0, 
            1, 0, False, np.zeros((1,1,1)), 1.0, 1.0)
            uncollided_class = uncollided_solutions.uncollided_solution(quick_build)

            output_maker = make_phi.make_output(tfinal, plotter.N_ang, plotter.ws, xs_list, plotter.coeff_mat, M, plotter.edges, uncollided_class)

            phi_new = output_maker.make_phi(uncollided)
            e_new = output_maker.make_e()
            
            data_phi[:, count] = phi_new
            data_e[:, count] = e_new
            plt.figure(fign)
            plt.plot(xs_list, phi_new)
            plt.show()
            fign+=1
    
    with open('su_olson_phi.csv', 'w') as file:
        writer = csv.writer(file)
        writer.writerows(data_phi)

    with open('su_olson_phi.e', 'w') as file:
        writer = csv.writer(file)
        writer.writerow(data_e)


 

