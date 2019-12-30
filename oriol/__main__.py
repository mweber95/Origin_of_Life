import argparse
import glob
import json
import numpy as np
import os
from Bio import Entrez
from oriol.tools.get_virus_names import ViralZone
from oriol.tools.download_and_parser import GetXml
from oriol.tools.download_and_parser import ParsingXml
from urllib import request as urlreq

parser = argparse.ArgumentParser()
parser.add_argument("--person", required=True, help="Who's that pokemon?", choices=["katrin", "fritz"])
parser.add_argument("--id", required=True, help="NCBI IDs or Genbank IDs", choices=["ncbi", "genbank"])

a = parser.parse_args()

if __name__ == "__main__":
    Entrez.email = 'YourID@hs-mittweida.de'

    katrin = {"baltimore_4": "https://viralzone.expasy.org/245", "baltimore_5": "https://viralzone.expasy.org/240",
              "baltimore_6": "https://viralzone.expasy.org/246"}

    fritz = {"baltimore_1": "https://viralzone.expasy.org/238"}


    """checking whether viral.txt exists or not --> if not download the file from ncbi"""
    if not os.path.exists("oriol/data/viral.txt"):
        print("Beginning file download...")
        url_ncbi = "https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=10239&cmd=download2"
        urlreq.urlretrieve(url_ncbi, "oriol/data/viral.txt")
        print("Finished file download!")

    if a.person == "katrin":
        viruses = ViralZone(katrin, a.id)
        virus_names = viruses.crawler()
        baltimore_4 = viruses.extracting_accession_numbers(virus_names[0])
        baltimore_5 = viruses.extracting_accession_numbers(virus_names[1])
        baltimore_6 = viruses.extracting_accession_numbers(virus_names[2])
    if a.person == "fritz":
        viruses = ViralZone(fritz, a.id)
        virus_names = viruses.crawler()
        baltimore_1 = viruses.extracting_accession_numbers(virus_names[0])

    if not os.path.exists("oriol/data/katrin"):
        os.mkdir("oriol/data/katrin")
    if not os.path.exists("oriol/data/katrin/ncbi_xml"):
        os.mkdir("oriol/data/katrin/ncbi_xml")
    if not os.path.exists("oriol/data/katrin/genbank_xml"):
        os.mkdir("oriol/data/katrin/genbank_xml")

    if not len(os.listdir("oriol/data/katrin/ncbi_xml")) == np.sum(np.array([len(baltimore_4), len(baltimore_5), len(baltimore_6)])):
        xml_download_baltimore_4 = GetXml(baltimore_4, "katrin", "ncbi")
        xml_download_baltimore_4.download()
        xml_download_baltimore_5 = GetXml(baltimore_5, "katrin", "ncbi")
        xml_download_baltimore_5.download()
        xml_download_baltimore_6 = GetXml(baltimore_6, "katrin", "ncbi")
        xml_download_baltimore_6.download()

    xml_files = glob.glob(os.path.join("oriol/data/katrin/ncbi_xml/", "*.xml"))

    splitted = [os.path.split(x) for x in xml_files]
    "Only NCBI optimized"
    ids = [x[:9] for _, x in splitted]

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

    os.mknod("oriol/data/katrin/ncbi.json")
    with open('oriol/data/katrin/ncbi.json', 'w') as f:
        json.dump(json_dict, f, indent=3)


