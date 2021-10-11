Scenario usage and construction
=================================

The :class:`Scenario` class is an object factory for :class:`SimContainer`.
It is and aggregate of two abstract factory classes (:class:`MyGlobalModel`, :class:`MyEntityLocator`),
one template object (:class:`ParameterTemplate`),
the parameter pool, entity types (instances of :class:`EntityType`)
and internal solver types (as types subclassed from :class:`InternalSovler`).

.. uml:: scenario_to_simcontainer.puml
      :caption: Factory classes are bundled in the :class:`Scenario` class

:class:`SimContainer` construction
------------------------------------

The scenario object contains information about the simulation scenario (in the form of template and factory objects),
but doesn't instantiate the costly simulation objects at creation.

Instead, when the actual :class:`SimConatiner` instance ins needed, in custom code or from within :class:`StateManager`,
the :py:func:`get_sim_container()` method gets invoked on the :class:`Scenario` instance (:numref:`get_sim_container`).
It returns the constructed :class:`SimContainer` instance,
which can run the simulation by invoking :py:func:`SimContainer.initialize()` and :py:func:`SimContainer.run()` methods.

This allows all simulation object to be hidden from the user,
though access/inspection is still possible entirely possible through a debugger.

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

:class:`Scenario` construction
-----------------------------------

To define a new model setup, a instance of the :class:`Scenario` class must
be populated with the appropriate template and factory objects.

The creation algorithm for the scenario, which is used in the examples,
is implemented in the setup function from box_grid.py.

In this case :py:func:`setup()` expects parameters as defined in the parameters.py
files of many models that use the box grid scenario,
but this is up to the programmer, the contents of box_grid.py can be anything that results in a valid scenario instance
and could be included in the main script.

In general the setup-routine:

    #. Defines several parameter templates for physical parameters used in the simulation
    #. parses a parameter file or recieves the customizable parameter values as arguments.
    #. Creates :class:`ParameterSet` instances for model-and entity-templates
    #. Creates and adds :class:`EntityType` instances
    #. Creates and adds :class:`FieldTemplate` and :class:`GlobalModel` instances
    #. Creates and adds one or more (non conflicting, please) :class:`EntityLocators`
    #. Returns the populated :class:`Scenario` instance


Template/factory objects contained in :class:`Scenario`
-------------------------------------------------------

:class:`ParameterTemplate`:

    A template object to derive independent :class:`Parameter` instances
    without redefining name an conversion factor/function

:class:`ParameterSet`:

    Holds the physical (with unit conversion) and misc (without conversion)
    :class:`Parameter` instances inside of :class:`ParameterCollection` objects.
    Typically each :class:`ParameterCollection` corresponds to a single
    field_quantity (il2, il6...) or a logical grouping of parameters (numeric ,geometry...).

:class:`EntityType`:

    Defines a "Type" that can be applied to a simulation entity,
    typically a cell type like Treg or Th1 to be applied to :class:`Cell` objects.
    It aggregates a parameter set and internal solver name to overwrite the
    values on :class:`Entity` object when applied.

:class:`FieldTemplate`:

    Object factory for :class:`GlobalProblem` objects. It stores name,
    field_quantity and a parameter set, specific for a given field.
    The choice of :class:`FieldTemplate` subclass decides what mathematical model is used
    (i.e. which subclass of :class:`GlobalProblem` will be instantiated)
    and what :class:`GlobalModel` objects it can be assigned to.

:class:`GlobalModel`:

    Aggregate class for :class:`FieldTemplate` objects, which together correspond to a given global model.

:class:`EntityLocator`:

    Object factory for :class:`Entity` objects. Takes some parameters and decides how and where to create simulation entities.
    Multiple locators can be used at once, but there is no checking
    wether the placement makes sense, i.e. geometries are non intersecting.



Entity-interactions (Boundary conditions)
-----------------------------------------




