# moving_mesh_radiative_transfer
## An accurate and fast moving mesh Discontinuous Galerkin package for solving the 1D isotropic transport equation with or without material energy coupling

### Quick start guide
### Installation 
Download the file ``moving_mesh_radiative_transfer``.
Navigate to file location. 

Invoke python via ``python3 -i imports.py``.

The program will automatically load the input script ``transport.yaml``. More on input scripts in the next section. 

`loader` is initiated with:
``load_sol(self, problem_name = 'transport', source_name = 'square_s', rad_or_transfer = 'rad', c = 1.0, s2 = False, cv0 = 0.0)``


To simply load precomputed transport results, 
``loader.call_sol(tfinal, M, x0_or_sigma, N_space, mat_or_rad, uncollided, moving)``

`tfinal`: Evaluation time

`M`: Number of basis functions -1

`x0_or_sigma`: Relevant for the square and Gaussian sources

`N_space`: Number of cells

'mat_or_rad': Either `'mat'` for material energy density or `'rad'` for radiation energy density

`uncollided`: `True` or `False`

`moving`: `True` or  `False`

Now, `loader.xs` returns the spatial points, `loader.phi` gives the scalar flux, and `loader.e` gives the material energy density. It is also possible to load the coefficients of the polynomial expansion and the quadrature weights for calculating phi, `loader.coeff_mat` and `loader.ws`.


To load results from other radiative transfer problems,
``loader = load(problem_name, source_name, rad_or_transfer, c, s2, cv0)``
where

`problem_name`: is one of `['transport','su_olson', 'su_olson_s2', 'su_olson_s2_thick', 'rad_transfer_const_cv', 'rad_transfer_const_cv_thick', 'rad_transfer_const_cv_s2', 'rad_transfer_const_cv_thick_s2']`

`source_name`: is one of `['plane_IC', 'square_IC', 'gaussian_IC', 'square_source', 'gaussian_source']`

`rad_or_transfer`: is 'rad' for all transport problems and 'transfer' for radiative transfer problems.

`c`: is the scattering ratio

`s2`: `True` or `False`

`cv0`: only applies to constant cv cases. 


### Solver

The solver is automatically imported as `run` after running `imports.py`. `transport` is the default problem, but it can be changed by

``run.load(problem_name)``

where `problem name` is one of `['transport','su_olson', 'su_olson_s2', 'su_olson_s2_thick', 'rad_transfer_const_cv', 'rad_transfer_const_cv_thick', 'rad_transfer_const_cv_s2', 'rad_transfer_const_cv_thick_s2']'

To actually solve the problem with the parameters given in the chosen input script, run one of:

``
solver.run_square_IC(uncollided = True, Moving = True)
``

``
solver.run_square_source(uncollided = True, Moving = True)
``

``
solver.run_gaussian_IC(uncollided = True, Moving = True)
``

``
solver.run_gaussian_source(uncollided = True, Moving = True)
``

``
solver.run_MMS(uncollided = False, Moving = True)

``
solver.run_plane_IC(uncollided = True, Moving = True)
``

Setting ``uncollided = False`` does not use the uncollided solution and ``moving = False`` solves the equations with a static mesh. (note: the plane pulse takes much longer to run compared to the other sources due to the higher number of discrete angles required to converge).


(note: there is only one case for the MMS source: `uncollided=False, moving=True' and only for transport porlbems)


runs all cases for every source.
The terminal will print a `RMSE` vallue (root mean square error) compared with the benchmark solution. The order of convergence is displayed as `Order`. The run data (RMSE, computation time, number of spaces, number of angles) will be saved to run_data_RMS.h5, overwriting previous runs that had the same parameters. A plot will be created of the benchmark solution and the solutions returned by the solver. 

To close the plots produced by the solver,

``plt.close()``


### Input scripts
The `YAML` input scripts, found in the folder, `moving_mesh_transport/input_scripts`, allow the user to easily modify the problem parameters without quitting python (Or re-compiling all of the functions that `numba` compiles on the first run.)  Simply change the parameters and `run.load('problem_name')` to load them. 

`mesh_parameters.yaml` contains parameters that apply to all problems. Most of the parameters have to do with the optically thick source mesh.  


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



### Final notes
In ``config.yaml`` certain parameters may be easily changed, such as the final evaluation time, number of basis functions, number of cells in the mesh, and number of discrete angles. Please ensure that the lengths of the  `N_spaces`, `Ms`, and `N_angles` lists are all the same, or the program will complain. 


### Testing
before invoking python, run ``pytest`` in the top level folder.
