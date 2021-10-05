Single simulation with SimContainer
====================================

Setup and running
-----------------------------

The following is a minimal example using the box_grid.py scenario.
It will run a single static simulation and save the cytokine field to file.


.. literalinclude:: ../../../../thesis/example_scripts/sim_container/run.py
    :linenos:

Line 6
    file path for simulation results
Line 7-14
    Call to setup function, which constructs the sim container, given some scenario specific parameters
Line 16
    Initially distributes cell types
Line 17
    Runs the simulation from timeindex 0 to 1
Line 18
    Saves the result of timeindex 1 to file
Line 19
    Saves the default markers (IL-2_surf_c and type_name) of timeindex 1 to file

Parameters for the box_grid scenario
------------------------------------

The run file imports names from parameters.py, which are required for scenario construction.
For box_grid.py a basic parameter file looks like this

.. literalinclude:: ../../../../thesis/example_scripts/sim_container/parameters.py
    :linenos:


Results
------------

The resulting file structure will look like this

.. code-block::

    ├── test_1
    │   ├── cache
    │   ├── entity_markers
    │   │   ├── IL-2_surf_c
    │   │   └── type_name
    │   ├── sol
    │   │   ├── distplot
    │   │   ├── field_IL-2_1.h5
    │   │   └── field_IL-2_1.xdmf
    │   └── solver_tmpil2

cache
    simulation mesh. In this case only one mesh is produced
entity_markers
    meshfunctions that can be loaded into paraview to visualize cell boundaries
sol
    solutions that can be visualized in paraview( .xmdf-files) or loaded again into fenics (.h5-files from distplot folder)
solver_tmpil2
    tmp directory for the solver processes. Sometimes useful for debugging




