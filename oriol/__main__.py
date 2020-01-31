import argparse
import glob
import json
import os
from Bio import Entrez
from oriol.tools.get_virus_names import ViralZone
from oriol.tools.download_and_parser import GetXml
from oriol.tools.download_and_parser import ParsingXml
from oriol.tools.download_and_parser import splitting_downloaded_xml
from urllib import request as urlreq

parser = argparse.ArgumentParser()
parser.add_argument("--id", required=True, help="NCBI IDs or Genbank IDs", choices=["ncbi", "genbank"])

a = parser.parse_args()

if __name__ == "__main__":
    Entrez.email = 'YourID@hs-mittweida.de'

    url_list = {"baltimore_4": "https://viralzone.expasy.org/245", "baltimore_5": "https://viralzone.expasy.org/240",
                "baltimore_6": "https://viralzone.expasy.org/246"}

    if not os.path.exists("oriol/data"):
        os.mkdir("oriol/data")

    """checking whether viral.txt exists or not --> if not download the file from ncbi"""
    if not os.path.exists("oriol/data/viral.txt"):
        print("Beginning file download...")
        url_ncbi = "https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=10239&cmd=download2"
        urlreq.urlretrieve(url_ncbi, "oriol/data/viral.txt")
        print("Finished file download!")

    viruses = ViralZone(a.id)
    virus_names = viruses.crawler(url_list)
    baltimore_4 = viruses.extracting_accession_numbers(virus_names[0])
    baltimore_5 = viruses.extracting_accession_numbers(virus_names[1])
    baltimore_6 = viruses.extracting_accession_numbers(virus_names[2])

    if not os.path.exists("oriol/data"):
        os.mkdir("oriol/data")
    if not os.path.exists("oriol/data/ncbi_xml"):
        os.mkdir("oriol/data/ncbi_xml")
    if not os.path.exists("oriol/data/genbank_xml"):
        os.mkdir("oriol/data/genbank_xml")

    xml_download_baltimore_4 = GetXml(baltimore_4, a.id)
    xml_download_baltimore_4.download("balt4")
    xml_download_baltimore_5 = GetXml(baltimore_5, a.id)
    xml_download_baltimore_5.download("balt5")
    xml_download_baltimore_6 = GetXml(baltimore_6, a.id)
    xml_download_baltimore_6.download("balt6")

    if a.id == "ncbi":
        list_of_xml = glob.glob(os.path.join("oriol/data/ncbi_xml/", "*.xml"))
    elif a.id == "genbank":
        list_of_xml = glob.glob(os.path.join("oriol/data/genbank_xml/", "*.xml"))

    list_of_xml.sort()

    if not os.path.exists("oriol/data/ncbi_xml/single_files"):
        os.mkdir("oriol/data/ncbi_xml/single_files")

    if not os.path.exists("oriol/data/genbank_xml/single_files"):
        os.mkdir("oriol/data/genbank_xml/single_files")

    splitting_downloaded_xml(a.id, list_of_xml)

    if a.id == "ncbi":
        xml_files = glob.glob(os.path.join("oriol/data/ncbi_xml/single_files", "*.xml"))
    elif a.id == "genbank":
        xml_files = glob.glob(os.path.join("oriol/data/genbank_xml/single_files", "*.xml"))

    splitted = [os.path.split(x) for x in xml_files]

    ids = [x[:-4] for _, x in splitted]

    baltimore = {}
    balti_4 = [key for key in baltimore_4]
    balti_5 = [key for key in baltimore_5]
    balti_6 = [key for key in baltimore_6]

    for element in balti_4:
        baltimore[element] = 4
    for element in balti_5:
        baltimore[element] = 5
    for element in balti_6:
        baltimore[element] = 6

    parser = ParsingXml(xml_files)
    definition = parser.definition()
    length = parser.length()
    lineage = parser.lineage()
    mol_type = parser.mol_type()
    cds = parser.cds()
    sequence = parser.sequence()

    json_dict = parser.builder(ids, baltimore, definition, length, lineage, mol_type, cds, sequence)

    if a.id == "ncbi":
        os.mknod("oriol/data/ncbi.json")
        with open('oriol/data/ncbi.json', 'w') as f:
            json.dump(json_dict, f, indent=3)
    elif a.id == "genbank":
        os.mknod("oriol/data/genbank.json")
        with open('oriol/data/genbank.json', 'w') as f:
            json.dump(json_dict, f, indent=3)




