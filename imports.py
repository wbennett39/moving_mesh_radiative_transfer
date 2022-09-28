# imports functions to run package from terminal 

from moving_mesh_transport.benchmarks import integrate_greens as intg
from moving_mesh_transport.plots import plotting_script as plotter
from moving_mesh_transport import solver
import matplotlib.pyplot as plt

# from moving_mesh_transport.plots.plot_square_s_times import main as plot_square_s_times
from moving_mesh_transport.solution_plotter import plot_thin_nonlinear_problems as plot_thin
from moving_mesh_transport.solution_plotter import plot_thick_nonlinear_problems as plot_thick
from moving_mesh_transport.solution_plotter import plot_su_olson as plot_su

from moving_mesh_transport.solution_plotter import make_tables_su_olson as tab_su

# from moving_mesh_transport.solver_classes.functions import test_s2_sol
from moving_mesh_transport.loading_and_saving.load_solution import load_sol as load




from moving_mesh_transport.solver_functions.run_functions import run

run = run()
run.load()

loader = load()



