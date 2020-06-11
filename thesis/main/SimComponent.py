import lxml.etree as et

from thesis.main.ParameterSet import ParameterSet


class SimComponent:

    def to_xml(self) -> et.Element:

        pass

    def from_xml(self, e: et.Element) -> None:

        pass
    def get_parameter(self):

        pass

    def update_parameter_set(self,parameter_set: ParameterSet):

        if self.p == None:

            self.p = parameter_set

        else:

            self.p.update_down(parameter_set)

