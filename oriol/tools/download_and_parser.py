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

    def ids(self):
        splitted = [os.path.split(x) for x in self.xml]
        ids = [x[:9] for _, x in splitted]
        return ids

    def baltimore(self, ids, old_baltimore):
        baltimore = "need to continue here"
        return baltimore

    def builder(self, ids, baltimore, length, lineage, cds, mol_type):
        pass


    #tree = ET.parse('country_data.xml')
    #root = tree.getroot()


#def download_xml(ncbi_id):
    #splitted_array = np.array_split(list(ncbi_id.keys()), int((len(ncbi_id)/400) + 1))
    #print(splitted_array)
    #all_ids = ",".join(list(ncbi_id.keys())[:400])
    #url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={all_ids}&format=xml"
    #urlreq.urlretrieve(url, f"data.xml")
