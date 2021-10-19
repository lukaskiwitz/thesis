Automation with :class:`StateManager` and post processing
===========================================================

:class:`StateManager` main usage is to automate parameter scans.
It is also necessary to use the post processing functions, even if no parameter scan is run.

Setup and running
-----------------------------

The previous script changes little if :class:`StateManager` is used.
The :class:`SimContainer` instance is now handed to :class:`StatManager` instance and the run method is called
on :class:`StateManager` rather than :class:`SimContainer`.
Additionally the function to initially distribute cell types has to be set to the pre_scan
hook on :class:`StateManager` and will be executed behind the scenes at an appropriate time.

.. literalinclude:: ../../../../thesis/example_scripts/minimal/state_manager/run.py
    :linenos:

Line 16
    instantiate :class:`StateManager`
Line 17
    hand :class:`SimContainer` instance to :class:`StateManager`
Line 18
    set time range
Line 20-23
    define and set pre_scan hook to randomly distribute cells at t = 0
Line 24
    run simulation


Results
-------

The result directory has changed a little.
Even though no parameter scan was explicitly defined :class:`StateManger` introduces an additional folder level.

.. code-block::

    └── test_1
    ├── log.xml
    ├── mesh
    │   ├── boundary_markers.h5
    │   ├── mesh_il2.h5
    │   └── mesh_il2.xdmf
    ├── records
    │   ├── dump.json
    │   ├── eta.npy
    │   └── ruse.h5
    ├── scan_0
    │   ├── m_0
    │   │   ├── replicat_0
    │   │   │   ├── entity_markers
    │   │   │   │   ├── IL-2_surf_c
    │   │   │   │   │   ├── marker_0.h5
    │   │   │   │   │   └── marker_0.xdmf
    │   │   │   │   └── type_name
    │   │   │   │       ├── marker_0.h5
    │   │   │   │       └── marker_0.xdmf
    │   │   │   └── sol
    │   │   └── solver_tmpil2
    │   ├── m_1
    │   │   └── replicat_0
    │   │       └── sol
    │   │           └── field_il2.npy
    │   └── timestep_logs
    │       ├── step_0_0_0_0.xml
    │       └── step_1_0_0_0.xml


log.xml
    xml file that contains information for post processing
mesh
    stores mesh and boundary markers if no external cache is used
records
    records on resource usage and execution timing collected during parameter scan
scan_0
    Output produced by :class:`SimContainer` for one scan sample plus the timestep_logs folder
m_0
    Output from model with index = 0 (pde model in this case)
replicat_0
    Result from replicat with index = 0
entity_marker
    MeshFunctions which store properties on boundary pieces (cell type, concentration, ...)


PostProcessing
----------------

Withe log.scan file we can use :class:`PostProcess` to get some useful output from the simulation.

.. literalinclude:: ../../../../thesis/example_scripts/minimal/state_manager/post_process.py
    :linenos:

Line 5
    Sets the `PARAVIEW_PATH` environment variable to point to the Paraview directory. It needs to be a headless client (omesa).
    The rendering script was written for Paraview version 5.9.0.
Line 8
    this needs to be set for correct unit conversion. The value -6 is hard coded in box_grid.py and simply needs to be repeat is.
    A better solution will come in the future.
Line 10
    Adds the `ParaviewRender` post processing option to the computation, which should be run during post processing.
    If Paraview is configured correctly, some standard visualizations will be rendered.
Line 11
    Runs the post processing with max. 4 subprocesses.

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
