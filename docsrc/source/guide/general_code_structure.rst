General Code Structure
======================

.. note ::
The following code example are mostly pseudo code extracted from the actual implementation and may be incomplete.

Simulation
-----------------
.. uml:: simulation.puml
    :caption: Classes involved in timeseries and parameter scan


The simulation of a single timeseries (with replicats) is handled by SimContainer. When its run method is invoked,
the basic simulation loop starts ( :ref:`basic_sim_loop`). It executes the pre/post hooks and but primarily calls the step
method for each timestep. SimContainer.step() iterates over all global problems (a.k.a fields),in turn calling their step methods,
and then over all entities.

.. code-block::
    :name: basic_sim_loop
    :caption: Basic simulation loop in pseudo code. Actual implementation in SimContainer.run() and step()

    def run(self, T, number_of_replicats = 1):

        for replicat_index in range(number_of_replicats):

            for field in self.global_problems:
                field.initialize_run(self.p, self.top_path, self.path, self.get_tmp_path())

            self._pre_replicat(self, 0 + 1, replicat_index, T[0 + 1], T)  # internal use
            self.pre_replicat(self, 0 + 1, replicat_index, T[0 + 1], T)  # user defined
            self.t = T[0]

            for time_index, t in enumerate(T[1:]):

                self._pre_step(self, time_index + 1, replicat_index, T[time_index + 1], T)
                self.pre_step(self, time_index + 1, replicat_index, T[time_index + 1], T)
                dt = t - T[time_index]
                self.step(dt, time_index, replicat_index)
                self._post_step(self, time_index, replicat_index, t, T)
                self.post_step(self, time_index, replicat_index, t, T)

            self._post_replicat(self, T[-1], replicat_index, t, T)
            self.post_replicat(self, T[-1], replicat_index, t, T)

    def step(self, dt, time_index , replicat_index):

        self.apply_type_changes(replicat_index)

        for field in self.global_problems:
            tmp_path = self.get_tmp_path()
            field.update_step(self.p, self.path, time_index, tmp_path)
            field.step(self.t, dt, time_index, tmp_path)

        for i, entity in enumerate(self.entity_list):
            entity.step(self.t, dt)

        self.t = self.t + dt

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

Solvers for global problems
---------------------------




Scenario to SimContainer
-------------------------

.. uml:: scenario_to_simcontainer.puml
      :caption: Factory classes are bundled in the Scenarion class



The Scenario class is and aggregate of three abstract factory classes (MyGlobalModel, MyEntityLocator, ParameterTemplate),
the parameter pool, entity types (instances of EntityType)
and internal solver types (as types subclassed from InternalSovler).

When the actual SimConatiner instance ins needed, in custom code or from within StateManager, the
get_sim_container() method gets invoked on the Scenario instance (:ref:`get_sim_container`).

.. code-block::
    :caption: get_sim_container method from Scenario
    :name: get_sim_container

    def get_sim_container(self, p, model_index):

        global_model = self.global_models[model_index]

        parameter_set = deepcopy(self.global_parameters)
        if p is not None:
            parameter_set.update(p, overwrite=True)

        for locator in self.entity_locators:
            cell_list = locator.get_entity_list(self.entity_types[0], parameter_set)

        parameter_set.update(global_model.build_parameter_set(self.parameter_pool))

        sc = SimContainer(parameter_set)

        for p in global_model.get_problem_list(parameter_set):
            sc.add_problem(p)

        for i in self.internal_solvers:
            sc.add_internal_solver(i)

        for c in cell_list:
            sc.add_entity(c)

        for entity_type in self.entity_types:
            sc.add_entity_type(entity_type)

        default = deepcopy(ScanSample(parameter_set.collections, self.entity_types, {}))
        sc.default_sample = default

        return sc

