
from moving_mesh_transport.benchmarks import integrate_greens as intg
from moving_mesh_transport.plots import plotting_script as plotter
from moving_mesh_transport import solver
import matplotlib.pyplot as plt

# from moving_mesh_transport.plots.plot_square_s_times import main as plot_square_s_times
from moving_mesh_transport.solution_plotter import plot_thin_nonlinear_problems as plot_thin
from moving_mesh_transport.solution_plotter import plot_thick_nonlinear_problems as plot_thick
from moving_mesh_transport.solution_plotter import plot_thick_suolson_problems as plot_sut
from moving_mesh_transport.solution_plotter import plot_su_olson as plot_su

from moving_mesh_transport.solution_plotter import make_tables_su_olson as tab_su

# from moving_mesh_transport.solver_classes.functions import test_s2_sol
from moving_mesh_transport.loading_and_saving.load_solution import load_sol as load



from moving_mesh_transport.solver_functions.run_functions import run



run = run()
run.load()

loader = load()

run.square_source()
plt.close()
plt.close()
plt.close()

# run.load('su_olson_thick')

t_list = [10.0, 31.6228, 100.0]
# t_list = [0.1, 0.31623, 1, 3.16228, 10.0]
# t_list = [0.3, 3.0, 30.0]
# factor_list = [2.5, 3.0, 6.0] # for thick gaussians nonlinear
# factor_list = [0.8, 0.9, 1.2]
# factor_list = [0.6, 0.6, 1.0, 1.0, 1.6, 1.6, 1.5] # thin gaussian ninlinear
factor_list = [2.4, 2.5, 1.7] # thin gaussian ninlinear

run.load('rad_transfer_const_cv')

for count, t in enumerate(t_list):
    run.parameters['all']['tfinal'] = t
    run.mesh_parameters['sqrt_t_move_factor'] = factor_list[count]
    run.parameters['all']['N_spaces'][0] = 64
    run.parameters['all']['Ms'][0] = 12
    run.gaussian_source(True, False)

    plt.close()
    plt.close()
    plt.close()

run.load('rad_transfer_const_cv_s2')

for count, t in enumerate(t_list):
    run.parameters['all']['tfinal'] = t
    run.mesh_parameters['sqrt_t_move_factor'] = factor_list[count]
    run.parameters['all']['N_spaces'][0] = 64
    run.parameters['all']['Ms'][0] = 12
    run.gaussian_source(True, False)

    plt.close()
    plt.close()
    plt.close()







