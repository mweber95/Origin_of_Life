import argparse
import json
import numpy as np
import re
import requests
import os
from Bio import Entrez, SeqIO
from bs4 import BeautifulSoup as bs
from urllib import request as urlreq


parser = argparse.ArgumentParser()
parser.add_argument("--id", required=True, help="Choose ID type!", choices=["ncbi", "genbank"])

a = parser.parse_args()

Entrez.email = 'weber7@hs-mittweida.de'


def get_virus_names(viral_zone):
    baltimore_groups = []
    for url in viral_zone:
        html_code = requests.get(url)
        soup = bs(html_code.content, "html.parser")
        c = [viruses for viruses in soup.getText().split() if any(x in viruses for x in ["virus", "viridae"])]
        baltimore_groups.append(c)
    return baltimore_groups


def extracting_accession_numbers(virus_names):
    accession_dict = {}
    with open("sequences/processing/viral.txt") as file:
        for line in file:
            for virus in virus_names:
                if virus in line:
                    if a.id == "ncbi":
                        ncbi_id = line.strip()[:9]
                        accession_dict[ncbi_id] = virus
                    elif a.id == "genbank":
                        genbank_id = re.split(r'\t+', line.strip())[1]
                        accession_dict[genbank_id] = virus
    return accession_dict


def download_xml(ncbi_id):
    splitted_array = np.array_split(list(ncbi_id.keys()), int((len(ncbi_id)/400) + 1))
    print(splitted_array)
    #all_ids = ",".join(list(ncbi_id.keys())[:400])
    #url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={all_ids}&format=xml"
    #urlreq.urlretrieve(url, f"data.xml")


if __name__ == "__main__":

    """checking whether folder sequences exists or not --> if not create needed folders"""
    if not os.path.exists("sequences"):
        os.mkdir("sequences")
    if not os.path.exists("sequences/processing"):
        os.mkdir("sequences/processing")
    if not os.path.exists("sequences/ncbi"):
        os.mkdir("sequences/ncbi")
    if not os.path.exists("sequences/genbank"):
        os.mkdir("sequences/genbank")

    if not os.path.exists("sequences/ncbi/ssRNA_plus"):
        os.mkdir("sequences/ncbi/ssRNA_plus")
    if not os.path.exists("sequences/ncbi/ssRNA_minus"):
        os.mkdir("sequences/ncbi/ssRNA_minus")
    if not os.path.exists("sequences/ncbi/ssRNA_plus_rev"):
        os.mkdir("sequences/ncbi/ssRNA_plus_rev")

    if not os.path.exists("sequences/genbank/ssRNA_plus"):
        os.mkdir("sequences/genbank/ssRNA_plus")
    if not os.path.exists("sequences/genbank/ssRNA_minus"):
        os.mkdir("sequences/genbank/ssRNA_minus")
    if not os.path.exists("sequences/genbank/ssRNA_plus_rev"):
        os.mkdir("sequences/genbank/ssRNA_plus_rev")

    """checking whether viral.txt exists or not --> if not download the file from ncbi"""
    if not os.path.exists("sequences/processing/viral.txt"):
        print("Beginning file download...")
        url_ncbi = "https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=10239&cmd=download2"
        urlreq.urlretrieve(url_ncbi, "sequences/processing/viral.txt")
        print("Finished file download!")

    """first function --> crawling ViralZone (https://viralzone.expasy.org/) for ssRNA virus names"""
    viral_zone_url = ["https://viralzone.expasy.org/245",
                      "https://viralzone.expasy.org/240",
                      "https://viralzone.expasy.org/246"]

    baltimore = get_virus_names(viral_zone_url)

    """second function --> searching the viral.txt with extracted virus names for NCBI or Genbank ID's"""

    accession_plus = extracting_accession_numbers(baltimore[0])
    accession_minus = extracting_accession_numbers(baltimore[1])
    accession_plus_rev = extracting_accession_numbers(baltimore[2])

    """third function --> downloading all sequences with biopython from ncbi"""

    download_xml(accession_plus)
    print(len(accession_plus))

    if a.id == "ncbi":
        folder_name_plus = "sequences/ncbi/ssRNA_plus/"
        folder_name_minus = "sequences/ncbi/ssRNA_minus/"
        folder_name_plus_rev = "sequences/ncbi/ssRNA_plus_rev/"

    if a.id == "genbank":
        folder_name_plus = "sequences/genbank/ssRNA_plus/"
        folder_name_minus = "sequences/genbank/ssRNA_minus/"
        folder_name_plus_rev = "sequences/genbank/ssRNA_plus_rev/"

    #download_fasta(accession_plus, folder_name_plus)
    #download_fasta(accession_minus, folder_name_minus)
    #download_fasta(accession_plus_rev, folder_name_plus_rev)

    """fourth function --> crawling each site for additional information"""



    #with open('data.json', 'w') as f:
        #json.dump(json_dict, f, indent=3)

# def download_fasta(accession_ids, foldername):
#     for i, accession_id in enumerate(accession_ids):
#         handle = Entrez.efetch(db="nuccore", id=accession_id, rettype="fasta", retmode="text")
#         os.mknod(foldername + str(accession_id) + ".fasta")
#         with open(foldername + str(accession_id) + ".fasta", "w") as file:
#             file.write(handle.read())
#         print(f"{i+1}/{len(accession_ids)} {accession_id}.fasta successfully created")


# def add_information(dict, baltimore_group):
#     url = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id="
#     for c, key in enumerate(dict):
#         html_code = requests.get(f"{url}{key}")
#         soup = bs(html_code.content, "html.parser")
#         plain_str = soup.getText()
#         lines = re.split(r'\n+', plain_str)
#         json_dict[key] = {"baltimore_group": baltimore_group}
#         json_dict[key]["CDS"] = {}
#         cds = []
#         for counter, element in enumerate(lines):
#             if "LOCUS" in element:
#                 length_regex = re.findall(r"[0-9]{3,5}.{1}[bp]{2}", element)
#                 length_regex_int = re.findall(r"[0-9]{3,5}", length_regex[0])
#                 if not key in json_dict:
#                     json_dict[key] = {"length": int(length_regex_int[0])}
#                 else:
#                     json_dict[key]["length"] = int(length_regex_int[0])
#             if "mol_type" in element:
#                 moltype_regex = re.findall(r"[genomicvral]{5,7}.{1}[cDRNA]{3,4}", element)
#                 if not key in json_dict:
#                     json_dict[key] = {"mol_type": moltype_regex[0]}
#                 else:
#                     json_dict[key]['mol_type'] = moltype_regex[0]
#             if "ORGANISM" in element:
#                 lineage_update = []
#                 element = lines[counter+1]
#                 lineage = element.split("; ")
#                 lineage[0] = lineage[0].replace(" ", "")
#                 if not re.search(r"\.", lineage[-1]):
#                     element = lines[counter+2]
#                     element = element[12:]
#                     lineage += element.split("; ")
#                 for virus_name in lineage:
#                     virus_name = virus_name.replace(";", "")
#                     virus_name = virus_name.replace(".", "")
#                     lineage_update.append(virus_name)
#                 if lineage[-2] == "unclassified":
#                     placeholder = f"{lineage_update[-2]} {lineage_update[-1]}"
#                     lineage_update = lineage_update[:-2]
#                     lineage_update.append(placeholder)
#                 if not key in json_dict:
#                     json_dict[key] = {"lineage": lineage_update}
#                 else:
#                     json_dict[key]['lineage'] = lineage_update
#             if " CDS " in element:
#                 element = element.split()
#                 cds.append(element[:][1])
#         #print(cds, key)
#         final_list = []
#         print(cds, key)
#         for element in cds:
#             element = element.replace("<", "")
#             if "join" in element:
#                 string = "join"
#                 tupel_pre = element[5:-1]
#                 tupel_pre = tupel_pre.split(",")
#                 new_list = []
#                 for area in tupel_pre:
#                     area = area.split("..")
#                     new_list.append((int(area[0]), int(area[1])))
#                     print(new_list)
#             elif "complement" in element:
#                 string = "complement"
#                 tupel_pre = element[11:-1].split("..")
#                 final_list.append([string, (int(tupel_pre[0]), int(tupel_pre[1]))])
#             else:
#                 tupel_pre = element.split("..")
#                 string = "regular"
#                 final_list.append([string, (int(tupel_pre[0]), int(tupel_pre[1]))])
#         print(final_list)
#             #element.replace("(", "")
#             #element.replace(")", "")
#             #print(element)
#
#                 cds = re.findall(r"\s+CDS\s+(join)?(complement)?\(?(\d+)\.\.(\d+)(,?(\d+)\.\.(\d+))?\)?", element)
#                 #print(f"{cds} ... {key}")

