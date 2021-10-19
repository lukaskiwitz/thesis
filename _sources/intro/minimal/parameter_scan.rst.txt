Parameter Scans with :class:`StateManager`
==========================================

To run the same simulation over a parameter range we need to extend the previous example
by adding and populating a :class:`ScanContainer` instance.

The steps to define a new scan are as follows:

    #. Retrieve the :class:`MyParameterPool` instance from the :class:`Scenario` instance
    #. Search for the desired a template by name
    #. Create a :class:`ScannableParameter` instance with a default parameter value. In this case derived from the parameter template
    #. Create a :class:`ScanDefintion` from the scannable parameter
    #. Add a parameter scan to the :class:`ScanContainer` instance

.. code-block::
    :linenos:

    parameter_pool = scenario.parameter_pool
    t_D = parameter_pool.get_template("D")
    D = ScannableParameter(t_D(10), lambda x, v: x * v)
    D_def = ScanDefintion(D, "IL-2", scan_space, ScanType.GLOBAL, field_quantity="il2")
    scan_container.add_single_parameter_scan([D_def], scan_name="D")


The classes involved corresponds to different function, that are split up to increase modularity.

:class:`MyParameterPool`

    Stores a list of named :class:`ParameterTemplates`. Some of the will have been created during scenario construction,
    but custom templates can be added.

:class:`ScannableParameter`

    Stores a reference to a :class:`Parameter` that serves as the default value for scanning
    and a function with the signature `lambda x,v: ...`, that defines how the parameter base value (`x`)
    is point in the scan space (`v`).

:class:`ScanDefintion`

    Stores the :class:`ScannableParameter` instance, the collection name and `field_quantity` the scan should be applied to,
    the range over which to scan (`scan_space`) and the type of sim object to apply the parameter to (:class:`ScanType`).

The point of this separation is flexibility in setting up parameter scans.
The sample scan definition can be used to create multiple scan.

.. code-block::
    :linenos:

    parameter_pool = scenario.parameter_pool
    t_D = parameter_pool.get_template("D")
    D = ScannableParameter(t_D(10), lambda x, v: x * v)
    D_def = ScanDefintion(D, "IL-2", scan_space, ScanType.GLOBAL, field_quantity="il2")

    t_q = parameter_pool.get_template("q")
    q = ScannableParameter(t_q(10), lambda x, v: x * v)
    q_def = ScanDefintion(q, "IL-2", scan_space, ScanType.GLOBAL, field_quantity="il2")

    scan_container.add_single_parameter_scan([D_def], scan_name="D")
    scan_container.add_single_parameter_scan([D_def, q_def], scan_name="D_q")

Here the `D_def` scan definition is used to set up two scans. The first, which scans only over `D` and a second one,
which scans over `D` and the secretion rate a the same time.


A running example is given below. It demonstrate how to scan over entity properties aswell.

.. literalinclude:: ../../../../thesis/example_scripts/minimal/parameter_scan/run.py
    :linenos:


