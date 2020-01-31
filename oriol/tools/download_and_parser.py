import os
import numpy as np
from urllib import request
import xml.etree.ElementTree as ET


def splitting_downloaded_xml(database, list_of_xml):
    for file in list_of_xml:
        tree = ET.parse(file)
        root = tree.getroot()
        for child in root:
            id = child.find("GBSeq_locus").text
            if database == "ncbi":
                filename = f"oriol/data/ncbi_xml/single_files/{id}.xml"
            elif database == "genbank":
                filename = f"oriol/data/genbank_xml/single_files/{id}.xml"
            os.mknod(filename)
            with open(filename, "w") as file:
                file.write(ET.tostring(child).decode("ASCII"))


class GetXml:
    def __init__(self, baltimore, database):
        self.baltimore = baltimore
        self.database = database

    def download(self, identifier):
        splitted_array = np.array_split(list(self.baltimore.keys()), np.ceil((len(self.baltimore) / 400)))
        for i, _ in enumerate(splitted_array):
            element_ids = ",".join(splitted_array[i])
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={element_ids}&format=xml"
            if self.database == "ncbi":
                request.urlretrieve(url, f"oriol/data/ncbi_xml/{identifier}_{i}.xml")
            elif self.database == "genbank":
                request.urlretrieve(url, f"oriol/data/genbank_xml/{identifier}_{i}.xml")
            print(f"{i + 1}/{len(splitted_array)} XML-files created")


class ParsingXml:
    def __init__(self, xml_files):
        self.xml = xml_files

    def definition(self):
        id_definition = {}
        for file in self.xml:
            "repeating element"
            _, file_id = os.path.split(file)
            id = file_id[:-4]
            tree = ET.parse(file)
            root = tree.getroot()
            id_definition[id] = root.find("GBSeq_definition").text
        return id_definition

    def length(self):
        id_length = {}
        for file in self.xml:
            "repeating element"
            _, file_id = os.path.split(file)
            id = file_id[:-4]
            tree = ET.parse(file)
            root = tree.getroot()
            id_length[id] = int(root.find("GBSeq_length").text)
        return id_length

    def lineage(self):
        id_lineage = {}
        for file in self.xml:
            "repeating element"
            _, file_id = os.path.split(file)
            id = file_id[:-4]
            tree = ET.parse(file)
            root = tree.getroot()
            id_lineage[id] = root.find("GBSeq_taxonomy").text.split("; ")
        return id_lineage

    def mol_type(self):
        id_moltype = {}
        for file in self.xml:
            "repeating element"
            _, file_id = os.path.split(file)
            id = file_id[:-4]
            tree = ET.parse(file)
            root = tree.getroot()
            id_moltype[id] = root.find("GBSeq_moltype").text
        return id_moltype

    def cds(self):
        id_cds = {}
        for file in self.xml:
            "repeating element"
            _, file_id = os.path.split(file)
            id = file_id[:-4]
            tree = ET.parse(file)
            root = tree.getroot()
            cds = []
            for gbtable in root.findall('GBSeq_feature-table'):
                for gbfeature in gbtable.findall("GBFeature"):
                    for location, feature in zip(gbfeature.findall("GBFeature_location"), gbfeature.findall("GBFeature_key")):
                        if feature.text == "CDS":
                            cds.append(location.text)
                id_cds[id] = cds
        return id_cds

    def sequence(self):
        id_sequence = {}
        for file in self.xml:
            "repeating element"
            _, file_id = os.path.split(file)
            id = file_id[:-4]
            tree = ET.parse(file)
            root = tree.getroot()
            id_sequence[id] = root.find("GBSeq_sequence").text.upper()
        return id_sequence

    def builder(self, ids, baltimore, definition, length, lineage, mol_type, cds, sequence):
        json_dict = {}
        for id in ids:
            json_dict[id] = {}

        for id in baltimore:
            for key in json_dict:
                if id == key:
                    json_dict[id]["baltimore_group"] = baltimore[id]

        for id in definition:
            for key in json_dict:
                if id == key:
                    json_dict[id]["definition"] = definition[id]

        for id in length:
            for key in json_dict:
                if id == key:
                    json_dict[id]["length"] = length[id]

        for id in lineage:
            for key in json_dict:
                if id == key:
                    json_dict[id]["lineage"] = lineage[id]

        for id in mol_type:
            for key in json_dict:
                if id == key:
                    json_dict[id]["mol_type"] = mol_type[id]

        for id in cds:
            for key in json_dict:
                if id == key:
                    json_dict[id]["CDS"] = cds[id]

        for id in sequence:
            for key in json_dict:
                if id == key:
                    json_dict[id]["sequence"] = sequence[id]

        return json_dict
