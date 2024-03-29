import os
from copy import deepcopy
from typing import Union, List, Mapping

from thesis.main.MyParameterPool import MyParameterPool
from thesis.main.ParameterSet import ParameterSet
from thesis.main.ScanContainer import ScanSample
from thesis.main.SimComponent import SimComponent
from thesis.main.SimContainer import SimContainer


# module_logger = logging.getLogger(__name__)

class ModelIndexOutOfRangeError(Exception): pass


class ModelNameNotFoundError(Exception): pass


class MyScenario(SimComponent):

    def __init__(self, parameter_pool):

        assert isinstance(parameter_pool, MyParameterPool)
        super(MyScenario, self).__init__()

        self.global_models = []
        self.internal_solvers = []
        self.entity_types = []
        self.entity_locators = []
        self.parameter_pool = parameter_pool
        self.global_parameters: ParameterSet = ParameterSet("scenario_dummy", [])
        self.markers = []
        self.marker_lookup: Mapping[str, int] = {}

    def get_model_indicies(self) -> List[int]:
        """Gets a list of all available model indices"""

        return list(range(len(self.global_models)))

    def get_model_name(self, model_index: int) -> str:
        """
        Gets model index for given model name if avaiable

        :raises ModelIndexOutOfRangeError
        """
        if model_index < len(self.global_models):

            return self.global_models[model_index].name
        else:
            raise ModelIndexOutOfRangeError("Theres is now global model with this index in this scenario")

    def get_model_index(self, model_name: str) -> int:
        """
        Returns index of first global model that matches with get_model_name

        :raises ModelNameNotFoundError
        """

        for i, model in enumerate(self.global_models):
            if model.name == model_name:
                return i
        raise ModelNameNotFoundError("There is no model with this name in this scenario")

    def get_sim_container(self, path: str, scan_sample: Union[ScanSample, None], model_index, ext_cache=""):

        if scan_sample is not None:
            p = scan_sample.p
            overwrite_locator_cache = scan_sample.remesh_scan_sample or scan_sample.remesh_timestep
        else:
            p = ParameterSet("dummy", [])
            overwrite_locator_cache = False
        global_model = self.global_models[model_index]

        parameter_set = deepcopy(self.global_parameters)
        if p is not None:
            parameter_set.update(p, overwrite=True)

        sc = SimContainer(parameter_set)
        sc.path = path
        sc.top_path = path
        parameter_set.update(global_model.build_parameter_set(self.parameter_pool))

        for p in global_model.get_problem_list(parameter_set):
            sc.add_problem(p)

        for i in self.internal_solvers:
            sc.add_internal_solver(i)

        if not ext_cache == "":
            sc.set_ext_cache(ext_cache)

        if hasattr(sc.global_problems[0], "get_mesh_path"):
            mesh_folder_path = "/".join(sc.global_problems[0].get_mesh_path(sc.path, 0, abspath=True).split("/")[:-1])
        else:
            mesh_folder_path = None

        for locator in self.entity_locators:
            cell_list = locator.get_entity_list(self.entity_types[0], parameter_set, os.path.abspath(sc.top_path),
                                                mesh_folder_path,
                                                overwrite_cache=overwrite_locator_cache)

        for c in cell_list:
            sc.add_entity(c)

        for entity_type in self.entity_types:
            sc.add_entity_type(entity_type)

        default = deepcopy(ScanSample(parameter_set.collections, self.entity_types, {}))
        sc.default_sample = default

        sc.marker_lookup = self.marker_lookup
        sc.markers = list(set(sc.markers + self.markers))

        return sc

    def get_entity_type_by_name(self, name: str):
        for e in self.entity_types:
            if e.name == name:
                return e

        return None

    def serialize_to_xml(self):
        raise NotImplementedError

    def deserialize_from_xml(self):
        raise NotImplementedError
