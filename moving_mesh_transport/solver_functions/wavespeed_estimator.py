"""
This notebook estimates the diffusive wavespeed of the 
scalar flux solution
"""
import numpy as np
import matplotlib.pyplot as plt

from ..solver_classes.make_phi import make_output
from ..solver_classes.functions import find_nodes
from ..solver_classes.make_phi import make_output


def wavespeed_estimator(sol, N_ang, N_space, ws, M, uncollided, mesh, uncollided_sol, thermal_couple, tfinal, x0):
    t_points = sol.t
    mesh.move(sol.t[-1])
    edges = mesh.edges
    # xs = find_nodes(edges, M)
    xs = np.linspace(-tfinal-x0,tfinal+x0,25 )
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


    for it in range(1, sol.t.size):
        for ix in range(phi.size-1):
            dt[it] += abs((solutions[it, ix] - solutions[it-1, ix])/delta_t)
            dx[it] += abs((solutions[it,ix+1]-solutions[it,ix-1])/(delta_x))
       
    
    return dt/(0.5*dx)





    

        






