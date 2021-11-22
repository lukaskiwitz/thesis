import hashlib
import json
from abc import abstractmethod, ABC
from copy import deepcopy
from typing import List, Union, Dict, Any, TypeVar, Callable

import dill
import lxml.etree as ET
import numpy as np

from thesis.main.MyError import DuplicateParameterError, DuplicateCollectionError
from thesis.main.my_debug import debug
from thesis.main.my_debug import message


class ParameterSet:
    """Container object that holds parameter collections

    :ivar collections: List of collection in this parameter set
    :ivar name: name of this parameter set

    :vartype collections: List[ParameterCollection]
    :vartype name: str
    """

    def __init__(self, name: str, collections: List['ParameterCollection']) -> None:

        for s in collections:
            assert isinstance(s, ParameterCollection)

        self.collections: List[ParameterCollection] = collections
        self.name: str = name

    def add_parameter_with_collection(self, parameter: 'Parameter'):
        """
        Adds a parameter with its own ParameterCollection

        :param parameter: parameter to add
        """
        dummy = ParameterSet("dummy", [ParameterCollection(parameter.name, [parameter])])
        self.update(dummy, overwrite=True)

    def update(self, input_parameter_object: Union['ParameterSet', 'ParameterCollection', 'Parameter'],
               overwrite: bool = False) -> None:

        """
        Updates this parameter set inplace. Values are not overwritten by naive.

        :param input_parameter_object: can be a Parameter, ParameterCollection or ParameterSet. When a Parameter object is passed a dummy collection is created.
        :param overwrite: overwrite values in this parameter set
        """

        if isinstance(input_parameter_object, Parameter):
            dummy = ParameterCollection(input_parameter_object.name, [input_parameter_object])
            debug("Collection {name} created with parameter".format(name=input_parameter_object.name))
            input_parameter_object = dummy

        if isinstance(input_parameter_object, ParameterCollection):
            dummy = ParameterSet("dummy", [input_parameter_object])
            input_parameter_object = dummy

        assert isinstance(input_parameter_object, ParameterSet)
        for update_c in input_parameter_object:
            if self.get_collection(update_c.name):

                self.get_collection(update_c.name).update(update_c, overwrite=overwrite)

            else:
                self.add_collection(deepcopy(update_c))
                debug("appending collection {n} to {n2}".format(n=update_c.name, n2=self.name))

    def add_collection(self, collection: 'ParameterCollection') -> None:
        """
        Adds a collection to this parameter set.

        :param collection: collection to add

        """

        if self.get_collection(collection.name) is None:
            self.collections.append(collection)
        else:
            self.get_collection(collection.name).update(collection, overwrite=True)

    def get_collection(self, name: str) -> Union['ParameterCollection', None]:

        """
        Retrieves as parameter collection by name.

        :param name: collection name to look for
        """

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

    def get_collections_by_field_quantity(self, field_quantity: str) -> List['ParameterCollection']:
        """
        Retrieves a list of parameter collections by associated field quantity.

        :param field_quantity: field quantity
        :return:
        """
        result = []
        for c in self:
            if c.field_quantity == field_quantity:
                result.append(c)

        return result

    def get_physical_parameter_by_field_quantity(self, parameter_name: str, field_quantity_name: str) -> Union[
        'PhysicalParameter', None]:

        """
        Retrieves a parameter by associated field quantity and name.

        :param parameter_name: parameter name
        :param field_quantity_name: field quantity
        """

        result = []
        collections = self.get_collections_by_field_quantity(field_quantity_name)

        for collection in collections:
            for parameter in collection:
                if parameter.name == parameter_name:
                    result.append(parameter)

        if len(result) > 1:
            raise DuplicateParameterError(field_quantity_name)

        elif len(result) == 0:

            return None

        else:
            return result[0]

    def __iter__(self):

        return self.collections.__iter__()

    def get_physical_parameter(self, parameter_name: str, collection_name: str) -> Union['PhysicalParameter', None]:
        """
        Retrieves a physical parameter by parameter and collection name.

        :param parameter_name:
        :param collection_name:
        """

        c = self.get_collection(collection_name)
        if c:
            return c.get_physical_parameter(parameter_name)
        else:
            return None

    def get_misc_parameter(self, parameter_name: str, collection_name: str) -> Union['MiscParameter', None]:

        """
       Retrieves a misc parameter by parameter and collection name.

       :param parameter_name:
       :param collection_name:
       """

        return self.get_collection(collection_name).get_misc_parameter(parameter_name)

    def serialize_to_xml(self, global_collections: 'GlobalCollections' = None,
                         global_parameters: 'GlobalParameters' = None) -> ET.Element:
        """
        Returns an xml representation of this parameter set.

        :param global_collections: global collection object (optional)
        :param global_parameters: global parameters object (optional)
        """

        root = ET.Element("ParameterSet")
        root.set("name", self.name)
        for c in self:
            if c.is_global and (not global_collections is None):
                root.append(global_collections._serialize_collection(c))
            else:
                root.append(c.serialize_to_xml(global_parameters=global_parameters))
        return root

    @staticmethod
    def deserialize_from_xml(element: ET.Element, parent_tree: ET.ElementTree = None) -> 'ParameterSet':

        """
        Creates a parameter set from xml representation

        :param element: the xml elment
        :param parent_tree: an xml tree to integrate global collections/parameters into (passed to ParameterCollection)
        """

        parameter_set = ParameterSet("dummy", [])

        parameter_set.name = element.get("name")
        for collection in element.findall("ParameterCollection"):
            c = ParameterCollection("", [])
            c.deserialize_from_xml(collection, parent_tree=parent_tree)
            parameter_set.add_collection(c)

        return parameter_set

    def get_as_dictionary(self, in_sim: bool = False, with_collection_name: bool = True, field_quantity: str = "") -> \
            Dict[str, Any]:

        """
        Returns dicitonary version of this parameter set.
        Collections are flattened into dict keys: {collectionname_parametername: value...}

        :param in_sim: return values in sim units
        :param with_collection_name: include collection name in dict keys
        :param field_quantity: only include collection with this field quantity
        """

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
                    try:
                        result[name] = parameter.get_in_sim_unit()
                    except:
                        pass
                else:
                    result[name] = parameter.get_in_post_unit()

        return result


class ParameterCollection:
    """
    Container object that holds parameter objects

    :ivar parameters: List of parameters in this collection
    :ivar name: collection name
    :ivar field_quantity: field quantity associated with this collection
    :ivar is_global: If True: parameter objects are handed down by reference. If False: parameter objects are handed down by copy

    :vartype: parameters: List[Parameter]
    :vartype: name: str
    :vartype: field_quantity: str
    :vartype: is_global: bool

    """

    def __init__(self, name: str, parameters: List['Parameter'] = None, field_quantity: str = "",
                 is_global: bool = False) -> None:

        assert isinstance(field_quantity, str)
        self.parameters: List[Parameter] = parameters
        self.name: str = name
        self.field_quantity = field_quantity
        self.is_global = is_global

    def __iter__(self):
        return self.parameters.__iter__()

    def get_as_dictionary(self, in_sim: bool = False, with_collection_name: bool = False) -> \
            Dict[str, Any]:

        """
        Returns dicitonary version of this collection.
        Collections are flattened into dict keys: {collectionname_parametername: value...}

        :param in_sim: return values in sim units
        :param with_collection_name: include collection name in dict keys
        :param field_quantity: only include collection with this field quantity
        """
        result = {}

        for parameter in self:
            if with_collection_name:
                name = self.name + "_" + parameter.name
            else:
                name = parameter.name
            if in_sim:
                try:
                    result[name] = parameter.get_in_sim_unit()
                except:
                    pass
            else:
                result[name] = parameter.get_in_post_unit()

        return result



    def update(self, update_collection: 'ParameterCollection', overwrite: bool = False):

        """
        Updates this parameter collection inplace. Values are not overwritten by naive.

        :param update_collection: collection with new parameters
        :param overwrite: overwrite values in this collection

        """
        for physical in update_collection:
            self.set_parameter(physical, overwrite=overwrite)
            debug("setting  parameter {n} on collection {n2}".format(n=physical.name, n2=self.name))

    def set_parameter(self, parameter: 'Parameter', overwrite: bool = False):

        """
        Sets a parameter in this collection

        :param parameter: parameter to set.
        :param overwrite: overwrite parameter with same name

        """

        if self.get_parameter(parameter.name):

            if parameter.is_global:
                self.parameters.remove(self.get_parameter(parameter.name))
                self.parameters.append(parameter)
            if overwrite:
                self.parameters.remove(self.get_parameter(parameter.name))
                self.parameters.append(deepcopy(parameter))

        else:
            if parameter.is_global:
                self.parameters.append(parameter)
            else:
                self.parameters.append(deepcopy(parameter))

    def set_physical_parameter(self, physical_parameter: 'PhysicalParameter', overwrite: bool = False):

        """
        Sets a physical parameter in this collection

        :param physical_parameter: physical parameter to set.
        :param overwrite: overwrite parameter with same name

        """

        assert type(physical_parameter) == PhysicalParameter

        self.set_parameter(physical_parameter, overwrite=overwrite)

    def set_misc_parameter(self, misc_parameter: 'MiscParameter', overwrite: bool = False):

        """
       Sets a mics parameter in this collection

       :param misc_parameter: misc parameter to set.
       :param overwrite: overwrite parameter with same name

       """

        assert type(misc_parameter) == MiscParameter
        self.set_parameter(misc_parameter, overwrite=overwrite)

    def get_parameter(self, parameter_name: str) -> Union['Parameter', None]:

        """
        Gets parameter by name

        :param parameter_name: parameter name
        """
        result = []

        for physical in self:
            if physical.name == parameter_name:
                result.append(physical)

        if len(result) > 1:
            raise DuplicateParameterError(parameter_name)

        elif len(result) == 0:

            return None

        else:
            return result[0]

    def get_physical_parameter(self, parameter_name: str) -> 'PhysicalParameter':

        """
        Gets physical parameter by name

        :param parameter_name: physical parameter name
        """

        result = self.get_parameter(parameter_name)
        if result is None:
            return None

        assert isinstance(result, PhysicalParameter)
        return result

    def get_misc_parameter(self, misc_name: str) -> 'MiscParameter':

        """
        Gets mics parameter by name

        :param misc_name: misc parameter name

        """

        result = self.get_parameter(misc_name)
        if result is None:
            return None

        assert isinstance(result, MiscParameter)
        return result

    def serialize_to_xml(self, global_parameters: 'GlobalParameters' = None) -> ET.Element:

        """Returns an xml representation of this parameter collection.

        :param global_parameters: global parameters object (optional)

        """
        root = ET.Element("ParameterCollection")

        root.set("name", json.dumps(self.name))
        root.set("field_quantity", json.dumps(self.field_quantity))
        root.set("is_global", json.dumps(self.is_global))

        for p in self:
            if p.is_global and (not global_parameters is None):
                root.append(global_parameters._serialize_parameter(p))
            else:
                root.append(p.serialize_to_xml())
        return root

    def deserialize_from_xml(self, element: ET.Element, parent_tree: ET.ElementTree = None) -> None:
        """
        Creates a parameter collection from xml representation

        :param element: the xml elment
        :param parent_tree: an xml tree to integrate global collections/parameters into

        """

        def replace_global(element, globals, parent_tree=None):

            def get_global_dict(e):

                collections = e.find(globals)
                if not collections is None:
                    return collections
                else:
                    return get_global_dict(e.getparent())

            if "global_key" in element.keys():

                if parent_tree is None:
                    global_collections = get_global_dict(element)
                else:
                    global_collections = get_global_dict(parent_tree)

                element = global_collections.find("{tag}[@key='{key}']".format(
                    tag=element.tag,
                    key=element.get("global_key")
                ))
                return element
            return element

        element = replace_global(element, "GlobalCollections", parent_tree=parent_tree)

        self.name = json.loads(element.get("name"))
        self.field_quantity = json.loads(element.get("field_quantity"))
        self.is_global = json.loads(element.get("is_global"))

        for physical in element.findall("PhysicalParameter"):
            physical = replace_global(physical, "GlobalParameters", parent_tree=parent_tree)
            p = PhysicalParameter("", 0)
            p.deserialize_from_xml(physical)
            self.parameters.append(p)

        for misc in element.findall("MiscParameter"):
            misc = replace_global(misc, "GlobalParameters", parent_tree=parent_tree)
            p = MiscParameter("", 0)
            p.deserialize_from_xml(misc)
            self.parameters.append(p)


T = TypeVar('T')
PhysicalValue = TypeVar('PhysicalValue', float, int, np.float, np.int, None)


class Parameter(ABC):
    """
    Abstract base class for container object that bundles parameter value, name and conversion function/factor.

    :ivar name: parameter name
    :ivar value: parameter value (must be at least json serializable)
    :ivar is_global: make this a global parameter. If True: parameter is handed down by reference. If False: parameter is handed down by copy

    :vartype name: str
    :vartype value: T
    :vartype is_global: bool
    """

    def __init__(self, name: str, value: T):
        self.name: str = name
        self.value: T = value
        self.is_global: bool = False

    def __call__(self, *args, in_sim: bool = False):

        """
        Gets parameter value

        :param args:
        :param in_sim:
        """
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

    def set_in_sim_unit(self, value: T):
        """
        Sets parameter value in sim units

        :param value: parameter value

        """
        self.value = value

    def get_in_sim_unit(self):
        """
        Gets parameter value in sim units
        """

        return self.value

    def set_in_post_unit(self, value: T):
        """
        Sets parameter value in post units

        :param value: parameter value

        """

        self.value = value

    def get_in_post_unit(self):
        """
        Gets parameter value in post units
        """

        return self.value

    def _get_serializiable_value(self) -> T:

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

    @abstractmethod
    def serialize_to_xml(self) -> ET.Element:
        """Returns an xml representation of this parameter"""
        pass

    @abstractmethod
    def deserialize_from_xml(self, element: ET.Element):
        """Creates a parameter from xml representation"""
        pass


class PhysicalParameter(Parameter):
    """
    Container class for numerical parameters that have physical units.

    :ivar to_sim: function of factor to convert from post to sim units
    :ivar to_post: function of factor to convert from sim to post units
    :ivar factor_conversion: wether conversion is by factor or callable

    :vartype to_sim: Union[float, Callable]
    :vartype to_post: Union[float, Callable]
    :vartype factor_conversion: bool

    """

    def __init__(self, name: str, value: PhysicalValue, to_sim: Union[float, Callable] = 1,
                 to_post: Union[float, Callable] = 1, is_global: bool = False) -> None:

        self.name: str = name
        self.to_sim: Union[float, Callable] = to_sim
        self.to_post: Union[float, Callable] = to_post
        self.is_global: bool = is_global
        self.factor_conversion: bool = False if callable(to_sim) else True
        self.value: PhysicalValue = None
        self.set_in_post_unit(deepcopy(value))

    def my_cast(self, value: Any) -> PhysicalValue:
        """Tries to cast value for serialization

        :param value: parameter value to cast
        """
        if not (isinstance(value, float) or isinstance(value, int)):
            message(""
                    "The value of physical parameter "
                    "{n} was type {t}. Casting to float".format(n=self.name, t=type(self.value))
                    )
            value = float(value)
        return value

    def get_in_sim_unit(self) -> PhysicalValue:

        return self.value

    def set_in_sim_unit(self, value: PhysicalValue) -> None:

        value = self.my_cast(value)

        self.value = value

    def get_in_post_unit(self) -> PhysicalValue:
        if self.factor_conversion:

            return self.value / self.to_sim
        else:
            return self.to_post(self.value)

    def set_in_post_unit(self, value: PhysicalValue):

        value = self.my_cast(value)

        if self.factor_conversion:

            self.value = value * self.to_sim

        else:
            self.value = self.to_sim(value)

    def serialize_to_xml(self) -> ET.Element:

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
    """
    Container class for miscellaneous parameters that do not have physical units.

    """

    def __init__(self, name: str, value: Any, is_global: bool = False) -> None:
        self.name: str = name
        self.is_global: bool = is_global
        self.value: Any = value

    def serialize_to_xml(self) -> ET.Element:
        root = ET.Element("MiscParameter")

        root.set("name", json.dumps(self.name))
        root.set("is_global", json.dumps(self.is_global))

        root.set("value", self._get_serializiable_value())

        return root

    def deserialize_from_xml(self, element: ET.Element):
        self.name = json.loads(element.get("name"))
        self.is_global = json.loads(element.get("is_global"))
        self.value = json.loads(element.get("value"))

    def get_in_sim_unit(self, type=str):
        value = super().get_in_sim_unit()
        return type(value)


class ScannableParameter:
    """
    Container class that combine parameter with a function for scanning over parameter ranges.

    :ivar p: base physical parameter
    :ivar f: scaling function; Must have signature f(x,v), were x is the base parameter value and v is the scaling parameter
    :ivar in_sim: wether to use sim or post units for scaling function

    :vartype p: Parameter
    :vartype f: Callable
    :vartype in_sim: bool
    """

    def __init__(self, p: Parameter, f: Callable, in_sim: bool = False):

        self.p = deepcopy(p)
        self.f = f
        self.in_sim = in_sim

    def __call__(self, v) -> Parameter:
        """
        Gets parameter object with value scaled according to scaling function f(x,v)

        :param v: scaling parameter
        """

        p = deepcopy(self.p)

        if self.in_sim:
            p.set_in_sim_unit(self.f(p.get_in_sim_unit(), v))
        else:
            p.set_in_post_unit(self.f(p.get_in_post_unit(), v))
        return p


class ParameterTemplate(ABC):
    """
    Abstract base class for parameter object factories

    :ivar name: template name

    :vartype name: str
    """

    @abstractmethod
    def __init__(self):
        self.name: Union[str, None] = None

    @abstractmethod
    def __call__(self, *args, **kwargs): pass


class PhysicalParameterTemplate(ParameterTemplate):
    """
    Object factory for physical parameters

    :ivar p: physical parameter prototype
    :ivar name: template name

    :vartype p: PhysicalParameter
    :vartype name: str
    """

    def __init__(self, parameter: PhysicalParameter):

        self.p: PhysicalParameter = parameter
        self.name: str = parameter.name

    def __call__(self, value: PhysicalValue = None, in_sim: bool = False) -> PhysicalParameter:
        """
        Returns physical parameter object

        :param value: parameter value
        :param in_sim: set value in sim units
        """
        p = deepcopy(self.p)
        if value is not None:
            if in_sim:
                p.set_in_sim_unit(value)
            else:
                p.set_in_post_unit(value)

        return p


class MiscParameterTemplate(ParameterTemplate):
    """
    Object factory for misc parameters

    :ivar p: misc parameter prototype
    :ivar name: template name

    :vartype p: MiscParameter
    :vartype name: str
    """

    def __init__(self, parameter: MiscParameter):
        self.p: MiscParameter = parameter
        self.name: str = parameter.name

    def __call__(self, value: Any = None) -> MiscParameter:
        """
        Returns physical parameter object

        :param value: parameter value
        """
        p = deepcopy(self.p)
        if value is not None:
            p.set_in_post_unit(value)

        return p


class GlobalParameters:
    """
    Container class ot store global parameters

    :ivar parameters: Dict containing stored parameters

    :vartype parameters: Dict[str, Parameter]
    """

    def __init__(self):
        self.parameters: Dict[str, Parameter] = {}

    def add(self, p: Parameter) -> None:

        """
        Adds parameter

        :param p: parameter Object
        """

        assert isinstance(p, Parameter)
        self.parameters[self.get_key(p)] = p

    def get_key(self, p: Parameter) -> str:
        """
        computes md5 hash on the string xml representation

        :param p: parameter to hash
        :return: parameter md5 hash
        """
        return hashlib.md5(ET.tostring(p.serialize_to_xml())).hexdigest()

    def serialize_to_xml(self) -> ET.Element:
        """
        :return: xml representation of this object
        """
        root = ET.Element("GlobalParameters")
        for i, p in self.parameters.items():
            e = p.serialize_to_xml()
            e.set("key", self.get_key(p))
            root.append(e)
        return root

    def _serialize_parameter(self, p: Parameter) -> ET.Element:
        """
        Serializes given parameter

        :param p: parameter to serialize
        :return: xml representation of the input parameter object
        """
        if not self.get_key(p) in self.parameters.keys():
            self.add(p)

        root = ET.Element(p.serialize_to_xml().tag)
        root.set("global_key", self.get_key(p))

        return root


class GlobalCollections:
    """
   Container class ot store global parameter collections

   :ivar collections: Dict containing stored parameter collections

   :vartype parameters: Dict[str, ParameterCollection]
   """

    def __init__(self):
        self.collections: Dict[str, ParameterCollection] = {}

    def add(self, c: ParameterCollection):
        """adds parameter collection"""

        assert isinstance(c, ParameterCollection)
        self.collections[self.get_key(c)] = c

    def get_key(self, c: ParameterCollection) -> str:

        """
        computes md5 hash on the string xml representation

        :param c: parameter collection to hash
        :return: collection md5 hash
        """
        return hashlib.md5(ET.tostring(c.serialize_to_xml())).hexdigest()

    def serialize_to_xml(self) -> ET.Element:
        """
        :return: xml representation of this object
        """
        root = ET.Element("GlobalCollections")
        for i, c in self.collections.items():
            e = c.serialize_to_xml()
            e.set("key", self.get_key(c))
            root.append(e)
        return root

    def _serialize_collection(self, c: ParameterCollection) -> ET.Element:

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
