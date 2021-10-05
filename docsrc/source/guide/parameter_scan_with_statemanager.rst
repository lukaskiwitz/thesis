StateManager
============

When SimContainer is used in conjunction with StateManger to run parameter scans and produce xml output for post processing,
the simulation loop is invoke from with StateManger.run() (:ref:`basic_scan_loop`).

.. code-block::
    :caption: Basic parameter scan loop in pseudo code. Actual implementation in StateManager.run()
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


For each combination of model and scan_sample the SimContainer.run() method is invoked once with the following steps

    #. (Lines 5-8) construct SimContainer instance from scenario
    #. (Line 9) Set file paths
    #. (Line 10) Initialize SimContainer instance
    #. (Line 11) update SimContainer instance with ScanSample instance
    #. (Line 12) run pre_scan hook
    #. (Line 13) call SimContainer.run()
    #. (Line 14) run post_scan hook
    #. (Line 15) write xml output
    #. (Line 16) Save timing records


The call to step-methods propagates until it reaches the MySolver.step() and InternalSolver.step() methods respectively.
The actual solver for GlobalProblems is implemented in the subclasses of MySolver and the solver for Entity behavior in
subclasses of InternalSolver.

.. uml:: state_manager_activity.puml

Solvers for global problems
---------------------------