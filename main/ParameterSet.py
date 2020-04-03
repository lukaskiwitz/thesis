import json
from copy import deepcopy
from typing import List

import dill
import lxml.etree as ET
import pandas as pd

from MyError import DuplicateParameterError, DuplicateCollectionError, CollectionNotFoundInParameterSet
from my_debug import debug


class ParameterSet:

    def __init__(self, name: str, collections) -> None:

        self.collections: List[ParameterCollection] = collections
        self.name = name
        self.parent: ParameterSet = None
        self.children: List[ParameterSet] = []

    def add_parameter_with_collection(self, parameter):

        dummy = ParameterSet("dummy", [ParameterCollection(parameter.name, [parameter])])
        self.update(dummy)

    def update(self, parameter_set, override = False) -> None:
        debug("updating {n1} with {n2}".format(n1=self.name, n2=parameter_set.name))
        assert isinstance(parameter_set, ParameterSet)
        for update_c in parameter_set:
            if self.get_collection(update_c.name):

                self.get_collection(update_c.name).update(update_c,override=override)

            else:
                self.add_collection(deepcopy(update_c))
                debug("appending collection {n} to {n2}".format(n=update_c.name, n2=self.name))

    def add_collection(self, collection):

        self.collections.append(collection)

    def get_collection(self, name: str):

        result = []

        for c in self:
            if c.name == name:
                result.append(c)

        if len(result) > 1:

            raise DuplicateCollectionError(name)

        elif len(result) == 0:

            return None

        else:
            return result[0]

    def get_collections_by_field_quantity(self, quantity: str):
        result = []

        for c in self:
            if c.field_quantity == quantity:
                result.append(c)
        return result

    def get_physical_parameter_by_field_quantity(self, parameter_name: str, quantity_name: str):

        result = []
        collections = self.get_collections_by_field_quantity(quantity_name)

        for collection in collections:
            for parameter in collection:
                if parameter.name == parameter_name:
                    result.append(parameter)

        if len(result) > 1:
            raise DuplicateParameterError(quantity_name)

        elif len(result) == 0:

            return None

        else:
            return result[0]

    def __iter__(self):

        return self.collections.__iter__()

    def get_physical_parameter(self, parameter_name: str, collection_name: str):

        c = self.get_collection(collection_name)
        if c:
            return c.get_physical_parameter(parameter_name)
        else:
            return None

    def get_misc_parameter(self, parameter_name: str, collection_name: str):

        return self.get_collection(collection_name).get_misc_parameter(parameter_name)

    def serialize_to_xml(self):

        root = ET.Element("ParameterSet")
        root.set("name", self.name)
        for c in self:
            root.append(c.serialize_to_xml())
        return root

    def deserialize_from_xml(self, element: ET.Element):

        self.name = element.get("name")
        for collection in element.findall("ParameterCollection"):
            c = ParameterCollection("", [])
            c.deserialize_from_xml(collection)
            self.add_collection(c)

    def get_as_dictionary(self):

        result = {}
        for collection in self:
            for parameter in collection:
                name = collection.name+"_"+parameter.name
                result[name] = parameter.get_in_post_unit()

        return result


class ParameterCollection:

    def __init__(self, name: str, physical_parameters=[], field_quantity="") -> None:

        self.parameters: List[PhysicalParameter] = physical_parameters
        self.name: str = name
        self.field_quantity = field_quantity
        self.is_global = False

    def __iter__(self):
        return self.parameters.__iter__()

    def update(self, update_collection, override = False):
        for physical in update_collection:
            self.set_parameter(physical, override=override)
            debug("setting  parameter {n} on collection {n2}".format(n=physical.name, n2=self.name))

    def set_parameter(self, parameter, override = False):

        if self.get_parameter(parameter.name):

            if parameter.is_global:
                self.parameters.remove(self.get_parameter(parameter.name))
                self.parameters.append(parameter)
            if override:
                self.parameters.remove(self.get_parameter(parameter.name))
                self.parameters.append(deepcopy(parameter))

        else:
            if parameter.is_global:
                self.parameters.append(parameter)
            else:
                self.parameters.append(deepcopy(parameter))

    def set_physical_parameter(self, physical_parameter):

        assert type(physical_parameter) == PhysicalParameter

        self.set_parameter(physical_parameter)

    def set_misc_parameter(self, misc_parameter):

        assert type(misc_parameter) == MiscParameter
        self.set_parameter(misc_parameter)

    def get_parameter(self, name: str):

        result = []

        for physical in self:
            if physical.name == name:
                result.append(physical)

        if len(result) > 1:
            raise DuplicateParameterError(name)

        elif len(result) == 0:

            return None

        else:
            return result[0]

    def get_physical_parameter(self, parameter_name: str):

        result = self.get_parameter(parameter_name)
        assert isinstance(result, PhysicalParameter)
        return result

    def get_misc_parameter(self, misc_name):

        result = self.get_parameter(misc_name)
        assert isinstance(result, MiscParameter)
        return result

    def serialize_to_xml(self):

        root = ET.Element("ParameterCollection")

        root.set("name", json.dumps(self.name))
        root.set("field_quantity", json.dumps(self.field_quantity))
        root.set("is_global", json.dumps(self.is_global))

        for p in self:
            root.append(p.serialize_to_xml())
        return root

    def deserialize_from_xml(self, element: ET.Element):

        self.name = json.loads(element.get("name"))
        self.field_quantity = json.loads(element.get("field_quantity"))
        self.is_global = json.loads(element.get("is_global"))

        for physical in element.findall("PhysicalParameter"):
            p = PhysicalParameter("", 0)
            p.deserialize_from_xml(physical)
            self.parameters.append(p)


class Parameter:

    def set_in_sim_unit(self, value):
        self.value = value

    def get_in_sim_unit(self):
        return self.value

    def set_in_post_unit(self, value):
        self.value = value

    def get_in_post_unit(self):
        return self.value


class PhysicalParameter(Parameter):

    def __init__(self, name: str, value: float, to_sim=1, to_post=1, is_global=False) -> None:

        self.to_sim = to_sim
        self.name = name
        self.to_post = to_post
        self.is_global = is_global
        self.factor_conversion = False if callable(to_sim) else True
        self.value = None
        self.set_in_post_unit(deepcopy(value))

    def get_in_sim_unit(self) -> float:

        return self.value

    def set_in_sim_unit(self, value: float) -> None:

        self.value = value

    def get_in_post_unit(self) -> float:
        if self.factor_conversion:

            return self.value / self.to_sim
        else:
            return self.to_post(self.value)

    def set_in_post_unit(self, value: float):

        if self.factor_conversion:

            self.value = value * self.to_sim

        else:
            self.value = self.to_sim(value)

    def serialize_to_xml(self):

        root = ET.Element("PhysicalParameter")
        root.set("name", json.dumps(self.name))
        root.set("factor_conversion", json.dumps(self.factor_conversion))
        root.set("is_global", json.dumps(self.is_global))
        root.set("in_sim_units", json.dumps(self.get_in_sim_unit()))
        root.set("in_post_units", json.dumps(self.get_in_post_unit()))
        root.set("to_sim", str(dill.dumps(self.to_sim).hex()))
        root.set("to_post", str(dill.dumps(self.to_post).hex()))

        return root

    def deserialize_from_xml(self, element: ET.Element):

        self.name = json.loads(element.get("name"))
        self.factor_conversion = json.loads(element.get("factor_conversion"))
        self.is_global = json.loads(element.get("is_global"))
        self.set_in_sim_unit(json.loads(element.get("in_sim_units")))
        self.to_sim = dill.loads(bytes.fromhex(element.get("to_sim")))
        self.to_post = dill.loads(bytes.fromhex(element.get("to_post")))


class MiscParameter(Parameter):

    def __init__(self, name: str, value, is_global=False) -> None:
        self.name = name
        self.is_global = is_global
        self.value = value
        self.field_quantity = ""

    def serialize_to_xml(self):
        root = ET.Element("MiscParameter")

        root.set("name", json.dumps(self.name))
        root.set("is_global", json.dumps(self.is_global))
        root.set("value", json.dumps(self.value))
        root.set("field_quantity", json.dumps(self.field_quantity))

        return root

    def deserialize_from_xml(self, element: ET.Element):
        self.name = json.loads(element.get("name"))
        self.is_global = json.loads(element.get("is_global"))
        self.value = json.loads(element.get("value"))
        self.field_quantity = json.loads(element.get("field_quantity"))


class ScannablePhysicalParameter:
    def __init__(self,p,f, in_sim = False):

        self.p = deepcopy(p)
        self.f = f
        self.in_sim = in_sim

    def __call__(self, v):

        p = deepcopy(self.p)

        if self.in_sim:
            p.set_in_sim_unit(self.f(p.get_in_sim_unit(), v))
        else:
            p.set_in_post_unit(self.f(p.get_in_post_unit(), v))
        return p


class PhysicalParameterTemplate():

    def __init__(self, parameter):
        self.p = parameter

    def __call__(self, value, in_sim=False):

        p = deepcopy(self.p)
        if in_sim:
            p.set_in_sim_unit(value)
        else:
            p.set_in_post_unit(value)

        return p
def make_collection(parameter_tupel):
    result = ParameterCollection(parameter_tupel[0], list(parameter_tupel[1]))
    if len(parameter_tupel) > 2:
        result.field_quantity = parameter_tupel[2]
    return result

def parse_parameter(collection_list,name):

    return ParameterSet(name, [make_collection(i) for i in collection_list])
