# moving_mesh_radiative_transfer
### An accurate and fast moving mesh Discontinuous Galerkin package for solving the 1D isotropic transport equation for the purpose of coupling to rad-transfer problems
**Quick start guide**

To solve the transport equation for a specific source, 

``
import moving_mesh_transport.solver
``

Running 

``
solver.run_plane_IC(uncollided = True, Moving = True)
``

Will read in parameters from moving_mesh_transport/congfig.yaml and run an infininte plane pulse source with a moving mesh and using the uncollided solution. Setting ``uncollided = False`` does not use the uncollided solution and ``moving = False`` solves the equations with a static mesh.
