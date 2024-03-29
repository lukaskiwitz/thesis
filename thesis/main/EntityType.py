import json
from abc import ABC
from copy import deepcopy
from typing import List

import lxml.etree as ET

from thesis.main.BC import BC
from thesis.main.ParameterSet import ParameterSet, ParameterCollection
from thesis.main.SimComponent import SimComponent


# module_logger = logging.getLogger(__name__)


class EntityType(ABC, SimComponent):

    def __init__(self, p: ParameterSet, name: str):

        super.__init__()

        self.p = deepcopy(p)
        self.name = name
        self.internal_solver_names: List[str] = []
        self.global_solver_names: List[str] = []
        self.interactions: List[BC] = []

    def get_updated(self, update):

        entity_type = deepcopy(self)
        if isinstance(update, ParameterCollection):
            update = [update]

        if isinstance(update, List):

            update_set = ParameterSet("dummy", update)

        elif isinstance(update, ParameterSet):

            update_set = update

        entity_type.p.update(update_set, overwrite=True)
        return entity_type


class CellType(EntityType):

    def __init__(self, p: ParameterSet, name: str, solver_name: str):
        super(EntityType, self).__init__()

        self.p = deepcopy(p)
        self.name = name
        self.internal_solver = solver_name
        self.interactions = []

    def serialize_to_xml(self):
        root = ET.Element("CellType")
        root.set("name", json.dumps(self.name))
        root.set("internal_solver_name", json.dumps(self.internal_solver))
        root.append(self.p.serialize_to_xml())

        return root

    def deserialize_from_xml(self, element: ET.Element):
        self.name = json.loads(element.get("name"))
        self.internal_solver = json.loads(element.get("internal_solver_name"))

        self.p = ParameterSet.deserialize_from_xml(element.find("ParameterSet"))
