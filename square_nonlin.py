# imports functions to run package from terminal 

from moving_mesh_transport.benchmarks import integrate_greens as intg
from moving_mesh_transport.plots import plotting_script as plotter
from moving_mesh_transport import solver
import matplotlib.pyplot as plt

from moving_mesh_transport.plots.plot_square_s_times import main as plot_square_s_times
from moving_mesh_transport.solution_plotter import plot as solplot
from moving_mesh_transport.solver_classes.functions import test_s2_sol



from moving_mesh_transport.solver_functions.run_functions import run

run = run()
# t_list = [0.1, 0.31623, 1, 3.16228, 10.0, 31.6228, 100.0]÷
t_list = [3.0]
# run.load('su_olson_s2')
# for t in t_list:
#     run.parameters['all']['tfinal'] = t
#     # run.mesh_parameters['choose_xs'] = True
#     run.square_source()
#     plt.close()
#     plt.close()
#     plt.close()
#     plt.close()


run.load('rad_transfer_const_cv_thick')
run.square_source(False, False)
run.mesh_parameters['estimate_wavespeed'] = False
run.mesh_parameters['find_wave_loc'] = False
for t in t_list:
    run.parameters['all']['tfinal'] = t
    run.parameters['square_source']['move_type'] = [0,1,0]

    # run.mesh_parameters['choose_xs'] = True
    run.square_source(False, True)
    plt.close()
    plt.close()
    plt.close()
    plt.close()
    
run.load('rad_transfer_const_cv_thick')
run.square_source(False, True)
for t in t_list:
    run.parameters['all']['tfinal'] = t
    run.parameters['square_source']['move_type'] = [0,1,0]
    # run.mesh_parameters['choose_xs'] = True
    run.mesh_parameters['estimate_wavespeed'] = False
    run.mesh_parameters['find_wave_loc'] = False
    run.square_source(False, True)
    plt.close()
    plt.close()
    plt.close()
    plt.close()





