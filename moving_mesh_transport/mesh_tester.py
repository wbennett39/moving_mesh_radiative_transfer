import numpy as np
import matplotlib.pyplot as plt
import quadpy

from .solver_classes import build_problem, mesh

def test_square_mesh(tfinal = 10.0, N_space = 8, x0 = 0.5, moving = True, move_type = np.array([1,0,0]), source_type = np.array([0, 0, 1]), edge_v = 1.0 / np.sqrt(3), thick = False, move_factor = 1.0, wave_loc_array = np.zeros((1,1,1)), pad = 25.0, leader_pad = 0.0 ):
    quad_thick_source = quadpy.c1.gauss_lobatto(int(N_space/2+1)).points
    quad_thick_edge = quadpy.c1.gauss_lobatto(int(N_space/4+1)).points

    mesh_ob =  mesh.mesh_class(N_space, x0, tfinal, moving, move_type, source_type, edge_v, thick, move_factor, wave_loc_array, pad, leader_pad, quad_thick_source, quad_thick_edge)

    time_points = 15 
    eval_times = np.linspace(0, tfinal, time_points)
    edges_array = np.zeros((time_points, N_space + 1))

    for count, t in enumerate(eval_times):
        mesh_ob.move(t)
        edges_array[count, :] = mesh_ob.edges

    print(np.linspace(-pad, pad, N_space + 1) - mesh_ob.edges, 'error')

    plt.ion()    
    plt.figure(1)
    for iedge in range(N_space + 1):
        mkr = 'bo--'
        if ( N_space / 4 <= iedge  <= 3 * N_space / 4):
            mkr = 'b^--'
        plt.plot(edges_array[:, iedge],eval_times, mkr , mfc = 'none')
        plt.plot(np.linspace(-pad, pad, 50), np.ones(50) * 10.0, 'k--')
    plt.scatter(np.linspace(-pad, pad, N_space+1), eval_times[-1] * np.ones(N_space+1), marker = '|', c= 'k', s= 64 )
    plt.xlabel('t')
    plt.ylabel('edge_pos')
    plt.show()


