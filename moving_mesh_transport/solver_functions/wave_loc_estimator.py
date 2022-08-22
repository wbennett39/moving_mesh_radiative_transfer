import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline
from ..solver_classes.make_phi import make_output
from ..solver_classes.functions import find_nodes



class find_wave:
    """
    This class takes solutions at an array of times, creates interpolated solutions and derivatives, 
    and estimates the wave location at those times
    """
    def __init__(self, N_ang, N_space, ws, M, uncollided, mesh, uncollided_sol, thermal_couple, tfinal, x0, times):
        self.N_ang = N_ang
        self.N_space = N_space
        self.ws = ws
        self.M = M
        self.uncollided = uncollided
        self.mesh = mesh
        self.uncollided_sol = uncollided_sol
        self.thermal_couple = thermal_couple
        self.tfinal = tfinal
        self.x0 = x0
        self.times = times
        self.dx = 1e-3 # step for searching for the wave 

    def find_wave(self, sol):
        self.make_sol(sol)
        left_edge_list = np.zeros(sol.t.size)
        right_edge_list = np.zeros(sol.t.size)
        for it in range(sol.t.size):
            interpolated_sol = UnivariateSpline(self.xs_list[it], self.solutions[it, :], k=4)
            self.first_deriv = interpolated_sol.derivative(1)
            self.second_deriv = interpolated_sol.derivative(2)
            xs_range = [0, self.xs_list[it,-1]]
            x_left, x_right = self.find_wave_bounds(xs_range)
            left_edge_list[it] = x_left
            right_edge_list[it] = x_right
            xs_test = np.linspace(self.xs_list[it,0], self.xs_list[it,-1], 1000)
            plt.figure(22)
            plt.plot(xs_test, self.first_deriv(xs_test), '-')
            plt.plot(xs_test, self.second_deriv(xs_test), '--')
        plt.show()

        return left_edge_list, right_edge_list


    def make_sol(self, sol):
        self.mesh.move(sol.t[-1])
        self.edges = self.mesh.edges
        xs = find_nodes(self.edges, self.M)
        self.solutions = np.zeros((sol.t.size, xs.size))
        self.xs_list = np.zeros((sol.t.size, xs.size))
        for it in range(sol.t.size):
            self.mesh.move(sol.t[it])
            self.edges = self.mesh.edges
            xs = find_nodes(self.edges, self.M)
            self.xs_list[it,:] = xs
            t = sol.t[it]
            self.mesh.move(t)
            if self.thermal_couple == 0:
                sol_reshape = sol.y[:,it].reshape((self.N_ang,self.N_space,self.M+1))
            elif self.thermal_couple == 1:
                sol_reshape = sol.y[:,it].reshape((self.N_ang+1,self.N_space,self.M+1))

            output = make_output(t, self.N_ang, self.ws, xs, sol_reshape, self.M, self.edges, self.uncollided)
            phi = output.make_phi(self.uncollided_sol)

            self.solutions[it, :] = phi

    def find_wave_bounds(self, xs_range):
        x_left = self.x0
        x_right = self.x0
        edge = xs_range[-1]
        left_found = False
        right_found = False
        inflection_found = False
        tol = 1e-13

        while left_found == False:
            x_left -= self.dx
            print(x_left)
            if x_left <= xs_range[0]:
                left_found = True
            elif abs(self.first_deriv(x_left)) <= tol:
                left_found = True

        
        while right_found == False:
            x_right += self.dx
            print(x_right)
            if x_right >= xs_range[-1]:
                right_found = True
            elif abs(self.first_deriv(x_right)) <= tol:
                right_found = True

        while inflection_found == False:
            second_old = self.second_deriv(edge)
            edge -= dx
            if np.sign(second_old) != np.sign(self.second_deriv(edge)):
                inflection_found = True
                print(edge)
            elif edge <= 0:
                inflection_found = True


        return x_left, x_right


