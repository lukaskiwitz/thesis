import json
from typing import List, Dict

import lxml.etree as ET

from thesis.main.EntityType import EntityType, CellType
from thesis.main.ParameterSet import ParameterSet, ParameterCollection, MiscParameter

class ScanContainer:

    def __init__(self):
        self.scan_samples = []

    def __iter__(self):
        return self.scan_samples.__iter__()

    def add_sample(self, sample):

        self.scan_samples.append(sample)

    def add_single_sim_parameter_scan(self, scanable, collection_name, field_quantity, scan_space, scan_name=None):


        for v in scan_space:
            sim_parameters = [
                ParameterCollection(collection_name, [scanable(v)], field_quantity=field_quantity),
            ]

            sample = ScanSample(sim_parameters, [], {}, scan_name=scan_name)
            self.scan_samples.append(sample)

    def add_multiple_sim_parameter_scan(self, scanable_collections, field_quantity, scan_space, scan_name=None):

        """

        :param scanable_collections: [(collection_name,field_quantity, [scanbles]), ...]
        """
        for v in scan_space:
            for c_name, fq, scanable in scanable_collections:

                sim_parameters = [
                    ParameterCollection(c_name,
                    [s(v) for s in scanable]
                    , field_quantity=fq),
                ]

            sample = ScanSample(sim_parameters, [], {}, scan_name=scan_name)
            self.scan_samples.append(sample)

    def add_single_entity_scan(self, entities,scanable,collection_name,field_quantity,scan_space, scan_name = None):

        for v in scan_space:
            entity_types = []
            for e in entities:
                if not field_quantity is None:
                    e = e.get_updated([ParameterCollection(collection_name, [scanable(v)], field_quantity=field_quantity)])
                else:
                    e = e.get_updated([ParameterCollection(collection_name, [scanable(v)])])
                entity_types.append(e)

            sample = ScanSample([], entity_types, {}, scan_name=scan_name)
            self.scan_samples.append(sample)


    def serialize_to_xml(self, scan_folder_pattern: str):

        root = ET.Element("ScanContainer")
        for n, s in enumerate(self.scan_samples):

            sample = s.serialize_to_xml(scan_folder_pattern.format(n=n))
            sample.set("scan_index", json.dumps(n))

            root.append(sample)

        return root

    def deserialize_from_xml(self,root):

        for n,s in enumerate(root.findall("ScanSample")):

            sample = ScanSample([],[],{})
            sample.deserialize_from_xml(s)
            self.scan_samples.append(sample)

class ScanSample:

    def __init__(self,collection_list: List[ParameterCollection], entity_types: List[EntityType], outer_domain_dict: Dict, scan_name = None):
        self.p: ParameterSet = ParameterSet("dynamic",collection_list)
        if not scan_name is None:
            name = MiscParameter("scan_name",scan_name)
            self.p.add_parameter_with_collection(name)

        self.entity_types: List[EntityType] = entity_types
        self.outer_domain_parameter_dict = {}
        for k,v in outer_domain_dict.items():
            # assert isinstance(v,ScannablePhysicalParameter)
            self.outer_domain_parameter_dict[k]  = ParameterSet("dummy",v)


    def serialize_to_xml(self, sub_path: str):

        root = ET.Element("ScanSample")

        root.set("sub_path",json.dumps(sub_path))

        parameters: ET.Element = ET.SubElement(root, "Parameters")
        parameters.append(self.p.serialize_to_xml())

        outer_domain_parameters: ET.Element = ET.SubElement(root, "OuterDomainParameters")
        for k,v in self.outer_domain_parameter_dict.items():

            patch_element = ET.SubElement(outer_domain_parameters,"PatchParameters")
            patch_element.set("key",k)
            patch_element.append(v.serialize_to_xml())

        entity_types: ET.Element = ET.SubElement(root, "EntityTypes")

        for entity_type in self.entity_types:
            entity_types.append(entity_type.serialize_to_xml())

        return root

    def deserialize_from_xml(self,root):


        self.p = ParameterSet("dummy",[])
        self.p.deserialize_from_xml(root.find("./Parameters/ParameterSet[@name='dynamic']"))

        for patch_element in root.findall("./OuterDomainParameters/PatchParameters"):
            k = patch_element.get("key")

            dummy = ParameterSet("dummy",[])
            dummy.deserialize_from_xml(patch_element[0])

            self.outer_domain_parameter_dict[k] = dummy


        for entity_type in root.find("EntityTypes"):
            if entity_type.tag == "CellType":
                t = CellType(None,"","")
                t.deserialize_from_xml(entity_type)
                self.entity_types.append(t)



