import os
from urllib import request
import xml.etree.ElementTree as ET


class GetXml:
    def __init__(self, ncbi_id, person, database):
        self.ncbi_id = ncbi_id
        self.person = person
        self.database = database

    def download(self):
        for i, id in enumerate(self.ncbi_id.keys()):
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={id}&format=xml"
            if self.person == "katrin":
                if self.database == "ncbi":
                    os.mknod(f"oriol/data/katrin/ncbi_xml/{id}.xml")
                    created_file = f"oriol/data/katrin/ncbi_xml/{id}.xml"
                if self.database == "genbank":
                    os.mknod(f"oriol/data/katrin/genbank_xml/{id}.xml")
                    created_file = f"oriol/data/katrin/genbank_xml/{id}.xml"
            request.urlretrieve(url, created_file)
            print(f"{i+1}/{len(self.ncbi_id.keys())} XML-files created")


class ParsingXml:
    def __init__(self, xml_files):
        self.xml = xml_files

    def definition(self):
        id_definition = {}
        for file in self.xml:
            "repeating element"
            _, file_id = os.path.split(file)
            id = file_id[:9]
            tree = ET.parse(file)
            root = tree.getroot()
            for part in root.findall('GBSeq'):
                definition = part.find("GBSeq_definition").text
                id_definition[id] = definition
        return id_definition

    def length(self):
        id_length = {}
        for file in self.xml:
            "repeating element"
            _, file_id = os.path.split(file)
            id = file_id[:9]
            tree = ET.parse(file)
            root = tree.getroot()
            for part in root.findall('GBSeq'):
                length = part.find("GBSeq_length").text
                id_length[id] = int(length)
        return id_length

    def lineage(self):
        id_lineage = {}
        for file in self.xml:
            "repeating element"
            _, file_id = os.path.split(file)
            id = file_id[:9]
            tree = ET.parse(file)
            root = tree.getroot()
            for part in root.findall('GBSeq'):
                taxonomy = part.find("GBSeq_taxonomy").text
                id_lineage[id] = taxonomy.split("; ")
        return id_lineage

    def mol_type(self):
        id_moltype = {}
        for file in self.xml:
            "repeating element"
            _, file_id = os.path.split(file)
            id = file_id[:9]
            tree = ET.parse(file)
            root = tree.getroot()
            for part in root.findall('GBSeq'):
                moltype = part.find("GBSeq_moltype").text
                id_moltype[id] = moltype
        return id_moltype

    def cds(self):
        id_cds = {}
        for file in self.xml:
            "repeating element"
            _, file_id = os.path.split(file)
            id = file_id[:9]
            tree = ET.parse(file)
            root = tree.getroot()
            for gbseq in root.findall("GBSeq"):
                cds = []
                for gbtable in gbseq.findall('GBSeq_feature-table'):
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
            id = file_id[:9]
            tree = ET.parse(file)
            root = tree.getroot()
            for part in root.findall('GBSeq'):
                sequence = part.find("GBSeq_sequence").text
                id_sequence[id] = sequence.upper()
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


    #tree = ET.parse('country_data.xml')
    #root = tree.getroot()


#def download_xml(ncbi_id):
    #splitted_array = np.array_split(list(ncbi_id.keys()), int((len(ncbi_id)/400) + 1))
    #print(splitted_array)
    #all_ids = ",".join(list(ncbi_id.keys())[:400])
    #url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={all_ids}&format=xml"
    #urlreq.urlretrieve(url, f"data.xml")