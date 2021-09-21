Runnig a full example
======================

1.Create a copy of one of the folders in "example_scripts/" in a location of your choice.

2.Open the parameters.py file and edit the path variable at the bottom to point to the desired top level directoy for your simulation results.

3.Run

.. code-block::

    python run.py; python post_process.py number_of_threads; python plot.py

to start the simulation and plot the results.

The first time any model is run it will need to

- create a mesh
- mark the boundary parts
- run the FFC compiler

this may take some time. For large meshes this can be considerable (1h or so), but
subsequent runs should be much faster.
