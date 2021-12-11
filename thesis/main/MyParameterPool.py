import logging
from typing import List, Mapping, Union

from thesis.main.ParameterSet import ParameterTemplate, ParameterCollection
from thesis.main.SimComponent import SimComponent

module_logger = logging.getLogger(__name__)


class MyParameterPool(SimComponent):
    """Stores parameter templates to derive Parameters and Collections from

    :ivar parameter_templates: List of parameter templates
    :vartype parameter_templates: List[ParameterTemplate]
    """

    def __init__(self):
        super(MyParameterPool, self).__init__()

        self.parameter_templates: List[ParameterTemplate] = []

    def join(self, pool: 'MyParameterPool', overwrite: bool = True):
        """
        Joins with a second parameter pool inplace.

        :param pool: Pool to join
        :param overwrite: Overwrite templates in this pool
        """

        own_names = [i.name for i in self.parameter_templates]

        for t in pool.parameter_templates:
            if (t.name in own_names and overwrite):
                for own_t in self.parameter_templates:
                    if t.name == own_t.name:
                        del own_t
                        break
                self.parameter_templates.append(t)

            elif t.name not in own_names:
                self.parameter_templates.append(t)

    def get_template(self, template_name: str) -> Union[ParameterTemplate, None]:

        """
        Retrieves first template that matches name

        :param template_name: name of the parameter template
        """
        for t in self.parameter_templates:
            if t.name == template_name:
                return t
        return None

    def add_template(self, template: ParameterTemplate):
        """
        Adds a template to this pool. Name must be unique.

        """
        assert isinstance(template, ParameterTemplate)
        assert self.get_template(template) is None

        self.parameter_templates.append(template)

    def get_as_collection(self, parameter_name_dict: Mapping[str, float], name: str = "dummy",
                          field_quantity: str = "") -> ParameterCollection:

        """
        Creates a parameter collection from key value pairs.

        :param parameter_name_dict: Keys must be template names in this pool and values must be appropriate parameter values
        :param name: collection name
        :param field_quantity: collection field quantity
        """
        c = ParameterCollection(name, [], field_quantity=field_quantity)
        for k, v in parameter_name_dict.items():
            t = self.get_template(k)
            if t is not None:
                if v is None:
                    c.set_parameter(t(in_sim=False), overwrite=True)
                else:
                    c.set_parameter(t(value=v, in_sim=False))
        return c
