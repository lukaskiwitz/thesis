Overview
==================

.. note::
    The code structure described in this section is subject to change.
    Many of the classes will undergo refactoring and some generalization/expansion, but the general design principles will be maintained,
    as well as most of the external interface. The most notable changes will occur in scenario construction and scan handling.


The major components of this library are

:class:`thesis.main.SimContainer.SimContainer`
    Handles simulation of a single timeseries
:class:`thesis.main.StateManager.StateManager`
    Handles parameters scans and xml-file creation
:class:`thesis.main.PostProcess.PostProcessor`
    Handles post processing of simulation results into more usable form
:class:`thesis.main.MyPlotter.Plotter`
    Utility class for plotting the results from post processing





