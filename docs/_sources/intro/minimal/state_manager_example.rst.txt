Automation with StateManager and post processing
================================================

:class:`thesis.main.StateManager.StateManager` main usage is to automate parameter scans.
It is also necessary to use the post processing functions, even if no parameter scan is run.

Setup and running
-----------------------------

The previous script changes little if StateManager is used.
The SimContainer instance is now handed to StatManager instance and the run method is called on StateManager rather than SimContainer.
Additionally the function to initially distribute cell types has to be set to the pre_scan hook on StateManager and will be executed behind the scenes at an appropriate time.

.. literalinclude:: ../../../../thesis/example_scripts/state_manager/run.py
    :linenos:

Line 16
    instantiate StateManager
Line 17
    hand SimContainer instance to StateManager
Line 18
    set time range
Line 20-23
    define and set pre_scan hook to randomly distribute cells at t = 0
Line 24
    run simulation


Results
-------

The result directory has changed a little.
Even though no parameter scan was explicitly defined StateManger introduces an additional folder level.

.. code-block::

    └── test_1
        ├── cache
        ├── log.scan
        ├── records
        └── scan_0
            ├── entity_markers
            │   ├── IL-2_surf_c
            │   └── type_name
            ├── sol
            │   └── distplot
            ├── solver_tmpil2
            │   └── pickle
            └── timestep_logs

log.scan
    xml file that contains information for post processing
records
    records on resource usage and execution timing collected during parameter scan
scan_0
    Output produced by SimContainer for one scan sample plus the timestep_logs folder

PostProcessing
----------------

Withe log.scan file we can use :class:`thesis.main.PostProcess` to get some useful output from the simulation.

.. literalinclude:: ../../../../thesis/example_scripts/state_manager/post_process.py
    :linenos:

Line 5
    this needs to be set for correct unit conversion. The value -6 is hard coded in box_grid.py and simply needs to be repeat is.
    A better solution will come in the future.
Line 6
    Runs the post processing with max. 4 subprocesses.
    If Paraview is configured as described in the readme, some standard visualizations will be rendered.

Processed Results
------------------

The post processing creates the following files/directories

cell_df.h5
    Pandas Dataframe of agent internal variables and positions

global_df.h5
    Pandas Dataframe of post processing computation results together with global parameter values.
    By default these are mean concentration, standard deviation and coefficient of variation over mesh vertices.

PostPorcess.xml
    Less user friendly but more structured representation of post processing results.
    An intermediate step in creating global_df.h5 and contains the same information.

timing_df.h5
    More usable version of records/dump.json which contains information about execution timing during simulation

images
    contains render result from Paraview

The dataframes can be loaded into a custom plotting script for downstream analysis and visualization of results.
