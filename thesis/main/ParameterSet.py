import hashlib
import json
from copy import deepcopy
from typing import List

import dill
import lxml.etree as ET

from thesis.main.MyError import DuplicateParameterError, DuplicateCollectionError
from thesis.main.my_debug import debug
from thesis.main.my_debug import message


class ParameterSet:

    def __init__(self, name: str, collections) -> None:

        for s in collections:
            assert isinstance(s, ParameterCollection)

        self.collections: List[ParameterCollection] = collections
        self.name = name
        self.parent: ParameterSet = None
        self.children: List[ParameterSet] = []

    def add_parameter_with_collection(self, parameter):

        dummy = ParameterSet("dummy", [ParameterCollection(parameter.name, [parameter])])
        self.update(dummy, override=True)

    def update(self, input, override=False) -> None:

        """
        input can now be a Parameter, ParameterCollection or ParameterSet.
        When a Parameter object is passed a dummy collection is created.

        """
        debug("updating {n1} with {n2}".format(n1=self.name, n2=input.name))

        if isinstance(input, Parameter):
            dummy = ParameterCollection(input.name, [input])
            debug("Collection {name} created with parameter".format(name=input.name))
            input = dummy

        if isinstance(input, ParameterCollection):
            dummy = ParameterSet("dummy", [input])
            input = dummy

        assert isinstance(input, ParameterSet)
        for update_c in input:
            if self.get_collection(update_c.name):

                self.get_collection(update_c.name).update(update_c, override=override)

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

    def serialize_to_xml(self, global_collections=None, global_parameters=None):

        root = ET.Element("ParameterSet")
        root.set("name", self.name)
        for c in self:
            if c.is_global and (not global_collections is None):
                root.append(global_collections.serialize_collection(c))
            else:
                root.append(c.serialize_to_xml(global_parameters=global_parameters))
        return root

    @staticmethod
    def deserialize_from_xml(element: ET.Element):

        parameter_set = ParameterSet("dummy", [])

        parameter_set.name = element.get("name")
        for collection in element.findall("ParameterCollection"):
            c = ParameterCollection("", [])
            c.deserialize_from_xml(collection)
            parameter_set.add_collection(c)

        return parameter_set

    def get_as_dictionary(self, in_sim=False, with_collection_name=True, field_quantity=""):

        result = {}
        for collection in self:
            if not field_quantity == "":
                if not collection.field_quantity == field_quantity:
                    continue

            for parameter in collection:
                if with_collection_name:
                    name = collection.name + "_" + parameter.name
                else:
                    name = parameter.name
                if in_sim:
                    result[name] = parameter.get_in_sim_unit()
                else:
                    result[name] = parameter.get_in_post_unit()

        return result


class ParameterCollection:

    def __init__(self, name: str, physical_parameters=[], field_quantity: str = "", is_global=False) -> None:

        for s in physical_parameters:
            assert isinstance(s, Parameter)
        assert isinstance(field_quantity, str)
        self.parameters: List[PhysicalParameter] = physical_parameters
        self.name: str = name
        self.field_quantity = field_quantity
        self.is_global = is_global

    def __iter__(self):
        return self.parameters.__iter__()

    def update(self, update_collection, override=False):

        for physical in update_collection:
            self.set_parameter(physical, override=override)
            debug("setting  parameter {n} on collection {n2}".format(n=physical.name, n2=self.name))

    def set_parameter(self, parameter, override=False):

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
        if result is None:
            return None

        assert isinstance(result, PhysicalParameter)
        return result

    def get_misc_parameter(self, misc_name):

        result = self.get_parameter(misc_name)
        if result is None:
            return None

        assert isinstance(result, MiscParameter)
        return result

    def serialize_to_xml(self, global_parameters=None):

        root = ET.Element("ParameterCollection")

        root.set("name", json.dumps(self.name))
        root.set("field_quantity", json.dumps(self.field_quantity))
        root.set("is_global", json.dumps(self.is_global))

        for p in self:
            if p.is_global and (not global_parameters is None):
                root.append(global_parameters.serialize_parameter(p))
            else:
                root.append(p.serialize_to_xml())
        return root

    def deserialize_from_xml(self, element: ET.Element):

        def replace_global(element, globals):

            def get_global_dict(e):

                collections = e.find(globals)
                if not collections is None:
                    return collections
                else:
                    return get_global_dict(e.getparent())

            if "global_key" in element.keys():
                global_collections = get_global_dict(element)
                element = global_collections.find("{tag}[@key='{key}']".format(
                    tag=element.tag,
                    key=element.get("global_key")
                ))
                return element
            return element

        element = replace_global(element, "GlobalCollections")

        self.name = json.loads(element.get("name"))
        self.field_quantity = json.loads(element.get("field_quantity"))
        self.is_global = json.loads(element.get("is_global"))

        for physical in element.findall("PhysicalParameter"):
            physical = replace_global(physical, "GlobalParameters")
            p = PhysicalParameter("", 0)
            p.deserialize_from_xml(physical)
            self.parameters.append(p)

        for misc in element.findall("MiscParameter"):
            misc = replace_global(misc, "GlobalParameters")
            p = MiscParameter("", 0)
            p.deserialize_from_xml(misc)
            self.parameters.append(p)


class Parameter:

    def __init__(self, name, value):
        self.name = name
        self.value = value

    def __call__(self, *args, in_sim=False):
        if len(args > 0):
            if in_sim:
                self.set_in_sim_unit(args[0])
            else:
                self.set_in_post_unit(args[0])
        else:
            if in_sim:
                return self.get_in_sim_unit()
            else:
                return self.get_in_post_unit()

    def set_in_sim_unit(self, value):
        self.value = value

    def get_in_sim_unit(self):
        return self.value

    def set_in_post_unit(self, value):
        self.value = value

    def get_in_post_unit(self):
        return self.value

    def _get_serializiable_value(self):

        value = "null"

        try:
            value = json.dumps(self.value)
        except TypeError as e:
            message(
                "Trying to serialize Parameter {n} as float because it was not json serializable".format(n=self.name))
            try:
                value = json.dumps(float(self.value))
            except Exception as e_2:
                message("Casting to float failed; cannot serialize Parameter {n}".format(n=self.name))

        return value


class PhysicalParameter(Parameter):

    def __init__(self, name: str, value: float, to_sim=1, to_post=1, is_global=False) -> None:

        self.to_sim = to_sim
        self.name = name
        self.to_post = to_post
        self.is_global = is_global
        self.factor_conversion = False if callable(to_sim) else True
        self.value = None
        self.set_in_post_unit(deepcopy(value))

    def my_cast(self, value):

        if not (isinstance(value, float) or isinstance(value, int)):
            message(""
                    "The value of physical parameter "
                    "{n} was type {t}. Casting to float".format(n=self.name, t=type(self.value))
                    )
            value = float(value)
        return value

    def get_in_sim_unit(self) -> float:

        return self.value

    def set_in_sim_unit(self, value: float) -> None:

        value = self.my_cast(value)

        self.value = value

    def get_in_post_unit(self) -> float:
        if self.factor_conversion:

            return self.value / self.to_sim
        else:
            return self.to_post(self.value)

    def set_in_post_unit(self, value: float):

        value = self.my_cast(value)

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
        self.to_sim = dill.loads(bytes.fromhex(element.get("to_sim")))
        self.to_post = dill.loads(bytes.fromhex(element.get("to_post")))
        self.set_in_sim_unit(json.loads(element.get("in_sim_units")))


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

        root.set("value", self._get_serializiable_value())
        root.set("field_quantity", json.dumps(self.field_quantity))

        return root

    def deserialize_from_xml(self, element: ET.Element):
        self.name = json.loads(element.get("name"))
        self.is_global = json.loads(element.get("is_global"))
        self.value = json.loads(element.get("value"))
        self.field_quantity = json.loads(element.get("field_quantity"))

    def get_in_sim_unit(self, type=str):
        value = super().get_in_sim_unit()
        return type(value)


class ScannablePhysicalParameter:
    def __init__(self, p, f, in_sim=False):

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


class ParameterTemplate:

    def __init__(self):
        pass


class PhysicalParameterTemplate(ParameterTemplate):

    def __init__(self, parameter):
        self.p = parameter
        self.name = parameter.name

    def __call__(self, value=None, in_sim=False):

        p = deepcopy(self.p)
        if value is not None:
            if in_sim:
                p.set_in_sim_unit(value)
            else:
                p.set_in_post_unit(value)

        return p


class GlobalParameters:

    def __init__(self):
        self.parameters = {}

    def add(self, p: Parameter):
        assert isinstance(p, Parameter)
        self.parameters[self.get_key(p)] = p

    def get_key(self, p: Parameter):
        return hashlib.md5(ET.tostring(p.serialize_to_xml())).hexdigest()

    def serialize_to_xml(self):
        root = ET.Element("GlobalParameters")
        for i, p in self.parameters.items():
            e = p.serialize_to_xml()
            e.set("key", self.get_key(p))
            root.append(e)
        return root

    def serialize_parameter(self, p: Parameter):

        if not self.get_key(p) in self.parameters.keys():
            self.add(p)

        root = ET.Element(p.serialize_to_xml().tag)
        root.set("global_key", self.get_key(p))

        return root


class GlobalCollections:

    def __init__(self):
        self.collections = {}

    def add(self, c: ParameterCollection):
        assert isinstance(c, ParameterCollection)
        self.collections[self.get_key(c)] = c

    def get_key(self, c: ParameterCollection):
        return hashlib.md5(ET.tostring(c.serialize_to_xml())).hexdigest()

    def serialize_to_xml(self):
        root = ET.Element("GlobalCollections")
        for i, c in self.collections.items():
            e = c.serialize_to_xml()
            e.set("key", self.get_key(c))
            root.append(e)
        return root

    def serialize_collection(self, c: ParameterCollection):

        if not self.get_key(c) in self.collections.keys():
            self.add(c)

        root = ET.Element(c.serialize_to_xml().tag)
        root.set("global_key", self.get_key(c))

        return root


def make_collection(parameter_tupel):
    result = ParameterCollection(parameter_tupel[0], list(parameter_tupel[1]))
    if len(parameter_tupel) > 2:
        result.field_quantity = parameter_tupel[2]
    return result


def parse_parameter(collection_list, name):
    return ParameterSet(name, [make_collection(i) for i in collection_list])
