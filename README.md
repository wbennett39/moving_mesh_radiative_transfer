# moving_mesh_radiative_transfer
## An accurate and fast moving mesh Discontinuous Galerkin package for solving the 1D isotropic transport equation for the purpose of coupling to rad-transfer problems
### Quick start guide
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
The terminal will print a RMSE vallue (root mean square error) compared with the benchmark solution. The run data (RMSE, computation time, number of spaces, number of angles) will be saved to run_data_RMS.h5, overwriting previous runs that had the same parameters.

To interact with the plots produced by the solver,

``
import matplotlib.pyplot as plt
``

run

``
plt.close()
``

to close the current plot.

### Plotter

