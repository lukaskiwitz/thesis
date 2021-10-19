Cell behavior with :class:`InternalSolver`
============================================
Cell behavior can be implemented as arbirtray python code by subclassing from :class:`InternalSolver` and
adding the `type` to the simulation.


The :term:`ABC` :class:`InternalSolver` defines two :term:`abstract methods` :py:func:`on_type_change`
and :py:func:`step`, which must be implemented in every subclass.
The custom behavior implemented inside the :py:func:`step` method.
It has access to the timepoints between which to solve (`t1` and `t2`), the :class:`ParameterSet` instance attached to
the cell (entity) and a reference to the :class:`Entity` instance.

.. note:: The later is a bit of a hold over and may be removed at some point.

.. code-block::
    :linenos:

    class SimpleThresholdSolver(InternalSolver):
    name = "SimpleThresholdSolver"

    def on_type_change(self, p, replicat_index, entity=None):
        pass

    def step(self, t1, t2, dt, p, entity=None, **kwargs):
        ...
        return p

Parameters can be retrived from the parameter set. This includes :

    - The coupling properties defined by :class:`GlobalProblem` instance in the simulation (`il2_surf_c`,...)
    - Parameters, which were set during a previous call to :py:func:`step`
    - The  parameters derived from the scenario.


.. code-block::
    :linenos:

    def step(self, t1, t2, dt, p, entity=None, **kwargs):

        """Retrieves the number of IL-2R molecules on this cell"""
        il2_threshold = p.get_physical_parameter("R", "IL-2").get_in_post_unit()

        """Retrieves the IL-2 surface concentration (coupling property) from last iteration"""
        il2 = p.get_physical_parameter("surf_c", "IL-2").get_in_post_unit()

        return p

A complete example is given below. It uses the `box_gird.py` scenario
so the model parameters are defined in `parameters.py`.
The starting cells can become either responder cell (`abs`)
or secreting cells (`sec`) depending on the IL-2 concentration on their surface.
The overall differentiation is modelled as a Poisson process.

As the example takes a few minutes run you might want to run the plotting code separately.

.. literalinclude:: ../../../../thesis/example_scripts/minimal/cell_behavior/run.py
    :linenos:

The result should look like this

.. image:: cell_behavior_plot.png
