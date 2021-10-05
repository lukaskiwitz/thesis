Time series simulation with SimContainer
========================================

The simulation of a single timeseries (with replicats) is handled by SimContainer.
At its core it is an aggregate of the Entity, EntityType, ParameterSet, GlobalProblem and InternalSolver classes,
to bundle simulation objects and provide an abstract interface for higher level code.
The SimContainer object can be build manually or, usually, retrieved from a scenario object.

The UML diagram (:ref:`test_uml`) illustrates the code structure around SimConatiner and its relation to StateManager.
The class hierarchy behind the abstract base classes has been omitted for clarity,
but that is where much of the actual implementation can be found.

.. uml:: simulation_classdiagram.puml
    :caption: Classes involved in timeseries and parameter scan
    :name: test_uml


When the SimContainer.run() method is invoked,
the basic simulation loop starts. It executes the pre/post hooks but primarily calls the step
method for each timestep. SimContainer.step() iterates over all global problems (a.k.a fields),in turn calling their step methods,
and then over all entities.
The program flow  is illustrated in the activity diagram (:ref:`basic_sim_loop_activity`) and pseudo code snippet (:ref:`basic_sim_loop_class`) below.

.. uml:: simcontainer_activity.puml
    :caption: activity diagram illustrating the basic simulation loop
    :name: basic_sim_loop_activity


.. code-block::
    :name: basic_sim_loop_class
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

The step methods of GlobalProblem (rather its subclasses) and Entity (same) actually advance
the simulation by one timestep for each global problem (pde-field) and each entity (cell model) respectively.

Global-Entity coupling
-------------------------

The coupling between entities and fields is iterative and occurs through coupling properties.
After a global problem is solved the coupling propertie(s) are calculated for each entity, that interacts with this problem
and stored in that entities parameter set, where it can be accessed by the entities internal solver.
However, there is now strong formalism associated with coupling properties,
the are simply entries in entity.p like any other parameter and are referenced by name.