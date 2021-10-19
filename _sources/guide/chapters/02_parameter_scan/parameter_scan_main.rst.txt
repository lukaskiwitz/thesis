StateManager
============

The main use of :class:`StateManager` is to run parameter scans in an organised manner and produces the xml
files necessary for the use of :class:`PostProcessor`.
When SimContainer is used in conjunction with StateManger to run parameter scans and produce xml output for post processing,
the simulation loop is invoke from with StateManger.run()
(:numref:`basic_scan_loop` and :numref:`state_manager_activity`).

.. uml:: state_manager_activity.puml
    :name: state_manager_activity
    :caption: Activity diagramm depicting the scan loop.

.. code-block::
    :caption: Basic parameter scan loop in pseudo code. Actual implementation in :py:func:`StateManager.run()`
    :name: basic_scan_loop
    :linenos:

    def run(self, ext_cache = "", model_names = None, number_of_replicats=1):

        for model_index in model_indicies:
            for scan_index in range(n_samples):
                self.sim_container = self.scenario.get_sim_container(
                    self.get_scan_sample(scan_index).p,
                    model_index=model_index
                )
                self.sim_container.top_path = self.path
                self.sim_container.initialize()
                self.update_sim_container(self.sim_container, scan_index, model_index)
                self.pre_scan(self, scan_index)
                self.sim_container.run(T, number_of_replicats=number_of_replicats)
                self.post_scan(self, scan_index)
                self.scan_tree.write_element_tree()
                self.save_records()


For each combination of model and scan_sample the :py:func:`SimContainer.run()` method is invoked once with the following steps

    #. (Lines 5-8) construct :class:`SimContainer` instance from scenario
    #. (Line 9) Set file paths
    #. (Line 10) Initialize :class:`SimContainer` instance
    #. (Line 11) update :class:`SimContainer` instance with :class:`ScanSample` instance
    #. (Line 12) run pre_scan hook
    #. (Line 13) call :py:func:`SimContainer.run()`
    #. (Line 14) run post_scan hook
    #. (Line 15) write xml output
    #. (Line 16) Save timing records


The call to step-methods propagates until it reaches the :py:func:`MySolver.step()` and :py:func:`InternalSolver.step()` methods respectively.
The actual solver for :class:`GlobalProblems` is implemented in the subclasses of :class:`MySolver` and the solver for :class:`Entity` behavior in
subclasses of :class:`InternalSolver`.


Parameter scan setup
------------------------------------------------------------------

.. uml:: scan_def_uml.puml
    :caption: The :class:`ScanContainer` class holds a list of :class:`ScanSample`, each representing a
        single point in the parameter scan. :class:`ScanSamples` objects are self contained
        and can be applied to :class:`SimContainer`, to change all (or just some)
        simulation parameters.
        :class:`ScanDefinition` instances can be passed to :class:`ScanContainer` as a more convenient
        way of adding scan points.

A parameter scan is produced by adding :class:`ScanSample` instances to :class:`ScanContainer`
and then handing the :class:`ScanContainer` instance to :class:`StateManager`.
This can be accomplished by either manually instantiating and population :class:`ScanSample`,
or through the use of :class:`ScanDefinition` and :py:func:`ScanContainer.add_single_parameter_scan`.

.. code-block::
    :caption: Code snippet that uses  :class:`ScanContainer` to setup a parameters scan over the
        Diffusion constant with logarithmically spaced samples ranging
        from parameter fold change 0.1 to 10.
    :name: D_scan_snippet

    template_D = pool.get_template("D")
    D = ScannableParameter(template_D(10), lambda x, v: x * v)
    scan_space = np.logspace(-1,1,10)
    D_def = ScanDefintion(D, "IL-2", scan_space, ScanType.GLOBAL, field_quantity="il2")
    scan_container.add_single_parameter_scan([D_def], scan_name="D")
    stMan.scan_container = scan_container


Instantiating :class:`ScanSamples` needs more code, but is more flexible.





