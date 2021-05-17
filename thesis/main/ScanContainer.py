import json
from typing import List, Dict
from enum import Enum

import lxml.etree as ET

from thesis.main.EntityType import EntityType, CellType
from thesis.main.ParameterSet import ParameterSet, ParameterCollection, MiscParameter, PhysicalParameterTemplate, ScannablePhysicalParameter



class ScanContainer:

    def __init__(self):
        self.scan_samples = []

    def __iter__(self):
        return self.scan_samples.__iter__()

    def add_sample(self, sample):

        self.scan_samples.append(sample)

    def add_single_parameter_scan(self, scan_list, scan_name="scan", add_axis = True, dynamic_mesh = False):

        """

        :param scan_list: [ScanDefintion])
        """

        if isinstance(scan_list,ScanDefintion):
            scan_list = [scan_list]



        scan_length = max([len(i.scan_space) for i in scan_list])


        for i in range(scan_length):

            sim_parameters = []
            entity_types = {}
            boundary = {}


            for d in scan_list:
                if i < len(d.scan_space):
                    v = d.scan_space[i]
                else:
                    v = d.scan_space[-1]

                assert isinstance(d,ScanDefintion)

                if add_axis:
                    axis  = ScannablePhysicalParameter(
                        PhysicalParameterTemplate(MiscParameter("value", 0, is_global=True))(0),
                        lambda x, v: v
                    )
                    sim_parameters.append(
                        ParameterCollection("scan", [axis(v)])
                    )

                if d.scan_type == ScanType.GLOBAL:
                    sim_parameters.append(
                        ParameterCollection(d.collection_name, [d.scannable(v)], field_quantity=d.field_quantity)
                    )
                elif d.scan_type == ScanType.ENTITY:

                    if d.entity_type.name not in entity_types.keys():
                        entity_type = d.entity_type
                    else:
                        entity_type = entity_types[d.entity_type.name]

                    entity_types[d.entity_type.name] = entity_type.get_updated(
                        [ParameterCollection(d.collection_name, [d.scannable(v)], field_quantity=d.field_quantity)]
                    )

                elif d.scan_type == ScanType.BOUNDARY:
                    boundary[d.boundary_pieces_name] = ParameterCollection(d.collection_name, [d.scannable(v)], field_quantity=d.field_quantity)


            sample = ScanSample(sim_parameters, list(entity_types.values()), boundary, scan_name=scan_name)
            sample.dynamic_mesh = dynamic_mesh
            self.scan_samples.append(sample)

    def add_2d_parameter_scan(self, a1, a2, scan_name = "scan"):

        """

        :param scan_list: (([ScanDefintion],axis_name), ([ScanDefintion],axis_name))

        """

        if len(a1) == 3:
            a1_axis = ScanDefintion(ScannablePhysicalParameter(
                PhysicalParameterTemplate(MiscParameter(a1[1], 0, is_global=True))(0),
                lambda x, v: v
            ), "axis", a1[2], ScanType.GLOBAL)
        else:
            a1_axis = ScanDefintion(ScannablePhysicalParameter(
                PhysicalParameterTemplate(MiscParameter(a1[1], 0, is_global=True))(0),
                lambda x, v: v
            ), "axis", range(max([len(i.scan_space) for i in a1[0]])), ScanType.GLOBAL)

        if len(a2) == 3:
            a2_axis = ScanDefintion(ScannablePhysicalParameter(
                PhysicalParameterTemplate(MiscParameter(a2[1], 0, is_global=True))(0),
                lambda x,v: v
            ),"axis",a2[2],ScanType.GLOBAL)
        else:
            a2_axis = ScanDefintion(ScannablePhysicalParameter(
                PhysicalParameterTemplate(MiscParameter(a2[1], 0, is_global=True))(0),
                lambda x, v: v
            ), "axis", range(max([len(i.scan_space) for i in a2[0]])), ScanType.GLOBAL)


        a1 = a1[0]
        a1.append(a1_axis)

        a2 = a2[0]
        a2.append(a2_axis)

        from copy import deepcopy
        length = max([len(i.scan_space) for i in a1])
        for i in range(length):

            a = deepcopy(a2)

            for d1 in a1:
                d = deepcopy(d1)
                d.scan_space = [d.scan_space[i]]
                a.append(d)

            self.add_single_parameter_scan(a, add_axis=False, scan_name=scan_name)

    def _add_single_entity_scan(self, entities,scanable,collection_name,field_quantity,scan_space, scan_name = None):

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

    def __init__(self,collection_list: List[ParameterCollection], entity_types: List[EntityType], outer_domain_dict: Dict, scan_name = "UnnamedScanSample"):

        self.p: ParameterSet = ParameterSet("dynamic",collection_list)
        self.dynamic_mesh = False

        if not scan_name is None:
            name = MiscParameter("scan_name",scan_name)
            self.p.add_parameter_with_collection(name)

        for e in entity_types:
            assert isinstance(e,EntityType)

        self.entity_types: List[EntityType] = entity_types
        self.outer_domain_parameter_dict = {}
        for k,v in outer_domain_dict.items():
            # assert isinstance(v,ScannablePhysicalParameter)
            self.outer_domain_parameter_dict[k]  = ParameterSet("dummy",[v])


    def serialize_to_xml(self, sub_path: str):

        root = ET.Element("ScanSample")
        dynamic_mesh = ET.SubElement(root, "DynamicMesh")
        dynamic_mesh.text = str(self.dynamic_mesh)

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


        self.p = ParameterSet.deserialize_from_xml(root.find("./Parameters/ParameterSet[@name='dynamic']"))
        dynamic_mesh = root.find("./DynamicMesh").text
        self.dynamic_mesh = True if dynamic_mesh == "True" else False

        for patch_element in root.findall("./OuterDomainParameters/PatchParameters"):
            k = patch_element.get("key")

            dummy = ParameterSet.deserialize_from_xml(patch_element[0])
            self.outer_domain_parameter_dict[k] = dummy


        for entity_type in root.find("EntityTypes"):
            if entity_type.tag == "CellType":
                t = CellType(None,"","")
                t.deserialize_from_xml(entity_type)
                self.entity_types.append(t)

class ScanType(Enum):

    GLOBAL = 1
    ENTITY = 2
    BOUNDARY = 3

class ScanDefintion:

    def __init__(self,scannable,collection_name,scan_space, scan_type: ScanType,
                 entity_type = None, boundary_pieces_name = None, field_quantity = ""):

        assert isinstance(scannable, ScannablePhysicalParameter)

        self.scannable = scannable
        self.collection_name = collection_name
        self.field_quantity = field_quantity
        self.scan_space = scan_space
        self.scan_type = scan_type
        self.entity_type = entity_type
        self.boundary_pieces_name = boundary_pieces_name


        if scan_type == ScanType.ENTITY:
            assert isinstance(entity_type,EntityType)
        elif scan_type == ScanType.BOUNDARY:
            assert isinstance(boundary_pieces_name, str)





