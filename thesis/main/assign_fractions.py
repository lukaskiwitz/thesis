import random

from thesis.main.ParameterSet import MiscParameter


def assign_fractions(sc, t):
    """sets cell types according to the values given in fractions.
        The pseudo random seed depends on t, so that cell placement is repeatable. """

    ran = random.Random()
    ran.seed(t)

    for i, e in enumerate(sc.entity_list):

        fractions = sc.p.get_collection("fractions")
        e.change_type = fractions.parameters[0].name

        draw = ran.random()
        s = 0
        for f in fractions.parameters[1:]:
            s = s + f.get_in_sim_unit()
            if draw < s:
                e.change_type = f.name
                break
        e.p.add_parameter_with_collection(MiscParameter("id", int(i)))
