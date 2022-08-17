"""
This notebook estimates the diffusive wavespeed of the 
scalar flux solution
"""
import numpy as np
import matplotlib.pyplot as plt
import math

from ..solver_classes.make_phi import make_output
from ..solver_classes.functions import find_nodes
from ..solver_classes.make_phi import make_output


def wavespeed_estimator(sol, N_ang, N_space, ws, M, uncollided, mesh, uncollided_sol, thermal_couple, tfinal, x0):
    t_points = sol.t
    mesh.move(sol.t[-1])
    edges = mesh.edges
    # xs = find_nodes(edges, M)
    xs = np.linspace(0, x0 + tfinal/(math.sqrt(3)) , 10000)
    solutions = np.zeros((sol.t.size, xs.size))
    wavespeeds = sol.t*0

    timesteps = sol.t[1:] - sol.t[:-1]

    # make solution at each time step
    for it in range(sol.t.size):
        
        t = sol.t[it]
        mesh.move(t)
        edges = mesh.edges
        # xs = find_nodes(edges, M)
        # xs = np.linspace(edges[0],edges[-1],1000)
        if thermal_couple == 0:
            sol_reshape = sol.y[:,it].reshape((N_ang,N_space,M+1))
        elif thermal_couple == 1:
            sol_reshape = sol.y[:,it].reshape((N_ang+1,N_space,M+1))

        output = make_output(t, N_ang, ws, xs, sol_reshape, M, edges, uncollided)
        phi = output.make_phi(uncollided_sol)

        solutions[it, :] = phi
    #     plt.figure(3)
    #     plt.plot(xs, phi, '-')
    #     plt.xlim(350,423)
    # plt.show()

    
    dx = np.zeros(sol.t.size)
    dt = np.zeros(sol.t.size)
    delta_t = sol.t[1] - sol.t[0]
    delta_x = abs(xs[1] - xs[0])

    left_edge_list = []
    right_edge_list = []

    left_edge_list.append(x0)
    right_edge_list.append(x0)
    ir_old = phi.size-1
    il_old = 0
    for it in range(1, sol.t.size):
        t = sol.t[it]
        spatial_deriv = derivative_estimator(xs,solutions[it,:],delta_x)
        il, ir = find_edge(t,xs, x0, spatial_deriv, il_old, ir_old)
        il_old = il
        ir_old = ir
        left_edge_list.append(xs[il])
        right_edge_list.append(xs[ir])
        for ix in range(phi.size-1):
            dt[it] += abs((solutions[it, ix] - solutions[it-1, ix])/delta_t)
            dx[it] += abs((solutions[it,ix+1]-solutions[it,ix-1])/(delta_x))
    # if right_edge_list != sol.t.size:
    #     print("mismatch")



    return dt/(0.5*dx), left_edge_list, right_edge_list




def derivative_estimator(xs,phi,delta_x):
    dx = phi*0
    dxx = phi*0
    for ix in range(1,phi.size-1):
        dx[ix] = (phi[ix] - phi[ix-1])/delta_x
    return dx
        
def find_edge(t,xs,x0, spatial_deriv, il_old, ir_old):
    x0_loc = np.argmin(np.abs(xs-x0))
    left_max_loc = np.argmin(np.abs(xs-(x0-t)))
    right_max_loc = np.argmin(np.abs(xs-(x0+t)))
    search_bounds = [left_max_loc, right_max_loc] 
    left = True
    right = True
    il = x0_loc
    ir = x0_loc
    mx = max(np.abs(spatial_deriv))
    mn = min(np.abs(spatial_deriv))
    tol = 1e-8



    while left == True:
        if abs(spatial_deriv[il]/mx) <= tol:
            left = False
        # elif il == search_bounds[0]:
        #     left = False
        else:
            il -= 1
    while right == True:
        if abs(spatial_deriv[ir]/mx) <= tol:
            right = False
        # elif ir == search_bounds[1]:
        #     right = False
        else:
            ir += 1
    return il, ir

    

        






