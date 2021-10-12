Single simulation with SimContainer
====================================

Setup and running
-----------------------------

The following is a minimal example using the box_grid.py scenario.
It will run a single static simulation and save the cytokine field to file.


.. literalinclude:: ../../../../thesis/example_scripts/minimal/sim_container/run.py
    :linenos:

Line 6
    file path for simulation results
Line 9-14
    Call to setup function, which constructs the sim container, given some scenario specific parameters
Line 19
    Gets sim container from scenario for model_index = 0 (pde model)
Line 20/21
    Set file paths for pde model
Line 23
    Initially distributes cell types
Line 25
    Initializes sim objects
Line 26
    Runs the simulation from timeindex 0 to 1
Line 27
    Saves the result of timeindex 1 to file
Line 30-41
    Repeats same process for ode model

Parameters for the box_grid scenario
------------------------------------

The run file imports names from parameters.py, which are required for scenario construction.
For box_grid.py a basic parameter file looks like this

.. literalinclude:: ../../../../thesis/example_scripts/minimal/sim_container/parameters.py
    :linenos:


Results
------------

The resulting file structure will look like this

.. code-block::

    └── test_1
        ├── ode_model
        │   └── sol
        └── pde_model
            ├── mesh
            │   ├── boundary_markers.h5
            │   ├── mesh_il2.h5
            │   └── mesh_il2.xdmf
            ├── sol
            │   ├── distplot
            │   ├── field_il2_0.h5
            │   └── field_il2_0.xdmf
            └── solver_tmpil2


cache
    simulation mesh. In this case only one mesh is produced
entity_markers
    meshfunctions that can be loaded into paraview to visualize cell boundaries
sol
    solutions that can be visualized in paraview( .xmdf-files) or loaded again into fenics (.h5-files from distplot folder)
solver_tmpil2
    tmp directory for the solver processes. Sometimes useful for debugging




