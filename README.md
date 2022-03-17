# moving_mesh_radiative_transfer
## An accurate and fast moving mesh Discontinuous Galerkin package for solving the 1D isotropic transport equation for the purpose of coupling to rad-transfer problems
### Quick start guide
### Installation 
Download the file ``moving_mesh_radiative_transfer`` and open the file. Invoke python via ``python3``.

### Testing
before invoking python, run ``pytest`` in the top level folder.

### Solver

To solve the transport equation for a specific source, 

``
import moving_mesh_transport.solver
``

Running 

``
solver.run_plane_IC(uncollided = True, Moving = True)
``

Will read in parameters from moving_mesh_transport/congfig.yaml and run an infininte plane pulse source with a moving mesh and using the uncollided solution. Setting ``uncollided = False`` does not use the uncollided solution and ``moving = False`` solves the equations with a static mesh. (note: the plane pulse takes much longer to run compared to the other sources due to the higher number of discrete angles required to converge).

The commands,

``
solver.run_square_IC(uncollided = True, Moving = True)
``

``
solver.run_square_s(uncollided = True, Moving = True)
``

``
solver.run_gaussian_IC(uncollided = True, Moving = True)
``

``
solver.run_gaussian_s(uncollided = True, Moving = True)
``

``
solver.run_MMS(uncollided = False, Moving = True)
``

run a square pulse, a square source, a Gaussian pulse, a Gaussian source, and a MMS (Method of Manufactured Solutions) problem.
(note: there is only one case for the MMS source)

The command 
``
solver.run_all()
``
runs all cases for every source.
The terminal will print a `RMSE` vallue (root mean square error) compared with the benchmark solution. The order of convergence is displayed as `Order`. The run data (RMSE, computation time, number of spaces, number of angles) will be saved to run_data_RMS.h5, overwriting previous runs that had the same parameters. A plot will be created of the benchmark solution and the solutions returned by the solver. 

To interact with the plots produced by the solver,

``
import matplotlib.pyplot as plt
``

run

``
plt.close()
``

to close the current plot.

### RMS Plotter

To visualize the accuracy of the results

``
from moving_mesh_transport.plots import make_plots
``

run 

``
make_plots.plot_all_rms_cells(tfinal, M)
``

To plot the results from running ``solver.run_all()`` where ``tfinal`` is the evaluation time and ``M`` is the number of basis functions. For the example case, choose ``tfinal = 1`` and `` M = 6``.

To plot a benchmark solutions, 

``
make_plots.plot_bench(tfinal, source_name, fign)
``

Where `source_name` can be `plane_IC`, `square_IC`, `square_source`, `gaussian_IC`, `gaussian_source`, or `MMS`.

### Benchmark maker

If you are interested in re-running the benchmark module (which takes quite a while)

``
from moving_mesh_transport.benchmarks import make_benchmarks
``

``
make_benchmarks.make_all()
``

will integrate the Greens function solution for the plane pulse and produce benchmark results for all of the sources run by the ``solver`` except for the MMS source at times 1, 5, and 10. For this reason, running the ``solver`` at other final times will return results but the benchmark solution will be 0 for all x. 

### Final notes
In ``config.yaml`` certain parameters may be easily changed, such as the final evaluation time, number of basis functions, number of cells in the mesh, and number of discrete angles. 
