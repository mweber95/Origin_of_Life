import os
import argparse
import glob
import re
import requests
import textwrap
from Bio import Entrez
from bs4 import BeautifulSoup as bs
from urllib import request as urlreq

parser = argparse.ArgumentParser()
parser.add_argument("--operation", required=False, help="operation to be executed",
                    choices=["download", "extract", "combine"])
parser.add_argument("--seq_type", required=False, help="which sequence typ should be filtered",
                    choices=["genome", "cds", "sequence"])

a = parser.parse_args()

Entrez.email = 'weber7@hs-mittweida.de'


def crawler(url, list_for_virus_names):
    html_code = requests.get(url)
    soup = bs(html_code.content, "html.parser")
    for tag in soup.find_all("a", {"href": re.compile(r"/by_species/[0-9]")})[7:]:
        list_for_virus_names.append(tag.getText())
    return list_for_virus_names


def extracting_accession_numbers(virus_names, accession_dict):
    with open("viral.txt") as file:
    #with open("sequences/processing/viral.txt") as file:
        for line in file:
            for virus in virus_names:
                if virus in line:
                    prepro = line.strip() #
                    ncbi_id = re.split(r'\t+', prepro)[1] #
                    print(ncbi_id) #
                    #ncbi_id = line.strip()[:9]
                    accession_dict[ncbi_id] = virus
    return accession_dict


def download_fasta(accession_ids, foldername, orientation):
    for accession_id in accession_ids:
        handle = Entrez.efetch(db="nuccore", id=accession_id, rettype="fasta", retmode="text")
        os.mknod(foldername + str(accession_id) + ".fasta")
        with open(foldername + str(accession_id) + ".fasta", "w") as file:
            file.write(str(accession_id + " " + accession_ids[accession_id] + " " + orientation + "\n" + handle.read()))
        print(str(accession_id) + ".fasta successfully created")


def convert_dna_to_rna(list_of_files, dict):
    for file in list_of_files:
        with open(str(file), "r+") as f:
            all_lines = f.readlines()
            header_art = all_lines[0]
            header_ori = all_lines[1]
            sequence_splitted = all_lines[2:]
            sequence_splitted_updated = [x[:-1] for x in sequence_splitted]
            sequence = "".join(sequence_splitted_updated)
            rna = "".join(dict.get(base) for base in sequence)
            new_fasta_rna = header_art + header_ori + textwrap.fill(rna, 70)
            f.seek(0)
            f.write(new_fasta_rna)


def extract(list_of_files):
    if a.seq_type == "genome":
        for file in list_of_files:
            with open(str(file), "r+") as f:
                all_lines = f.readlines()
                header = all_lines[0]
                checking_type = all_lines[1]
                sequence = all_lines[2:]
                sequence_updated = "".join(sequence)
                if a.seq_type in checking_type:
                    os.mknod("sequences/genome/" + os.path.split(file)[1])
                    with open(str("sequences/genome/" + os.path.split(file)[1]), "w") as new_f:
                        new_f.write(">" + header + textwrap.fill(sequence_updated, 70))


def combine_files(list_of_files):
    os.mknod('sequences/genome/all.fasta')
    with open('sequences/genome/all.fasta', 'w') as outfile:
        for file in list_of_files:
            with open(file) as infile:
                for line in infile:
                    outfile.write(line)
                outfile.write("\n\n")


def main():

    if a.operation == "download":

        """checking whether folder sequences exists --> if not create needed folders"""
        if not os.path.exists("sequences"):
            os.mkdir("sequences")
        if not os.path.exists("sequences/processing"):
            os.mkdir("sequences/processing")
        if not os.path.exists("sequences/ssRNA_plus"):
            os.mkdir("sequences/ssRNA_plus")
        if not os.path.exists("sequences/ssRNA_minus"):
            os.mkdir("sequences/ssRNA_minus")

        """checking whether viral.txt exists or not --> if not download the file from ncbi"""
        if not os.path.exists("sequences/processing/viral.txt"):
            print("Beginning file download...")
            url_ncbi = "https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=10239&cmd=download2"
            urlreq.urlretrieve(url_ncbi, "sequences/processing/viral.txt")
            print("Finished file download!")

        """first function --> crawling viralzone (https://viralzone.expasy.org/) for ssRNA virus names"""
        list_for_virus_names_plus = []
        list_for_virus_names_minus = []
        url_list = ["https://viralzone.expasy.org/245", "https://viralzone.expasy.org/240"]
        virus_names_plus = crawler(url_list[0], list_for_virus_names_plus)
        virus_names_minus = crawler(url_list[1], list_for_virus_names_minus)

        """second function --> searching the viral.txt with extracted virus names"""
        accession_dict_plus = {}
        accession_dict_minus = {}
        if not os.path.exists("sequences/processing/accession_numbers.txt"):
            os.mknod("sequences/processing/accession_numbers.txt")
        unique_accession_plus = extracting_accession_numbers(virus_names_plus, accession_dict_plus)
        unique_accession_minus = extracting_accession_numbers(virus_names_minus, accession_dict_minus)

        """third function --> downloading all sequences with biopython from ncbi"""
        if not len(os.listdir("sequences/ssRNA_plus")) == len(unique_accession_plus):
            foldername_plus = "sequences/ssRNA_plus/"
            orientation_plus = "(+)"
            download_fasta(unique_accession_plus, foldername_plus, orientation_plus)
        if not len(os.listdir("sequences/ssRNA_minus")) == len(unique_accession_minus):
            foldername_minus = "sequences/ssRNA_minus/"
            orientation_minus = "(-)"
            download_fasta(unique_accession_minus, foldername_minus, orientation_minus)

        """fourth function --> convert dna to rna"""
        all_files_plus = glob.glob(os.path.join("sequences/ssRNA_plus", "*.fasta"))
        all_files_minus = glob.glob(os.path.join("sequences/ssRNA_minus", "*.fasta"))

        dict = {"A": "A", "T": "U", "G": "G", "C": "C",
                "W": "W", "S": "S", "M": "M", "K": "K", "R": "R", "Y": "Y",
                "B": "B", "V": "V", "D": "D", "H": "H",
                "N": "N", "Z": "Z"}

        convert_dna_to_rna(all_files_minus, dict)
        convert_dna_to_rna(all_files_plus, dict)

    if a.operation == "extract":

        if not os.path.exists("sequences/genome"):
            os.mkdir("sequences/genome")
        if not os.path.exists("sequences/cds"):
            os.mkdir("sequences/cds")
        if not os.path.exists("sequences/sequence"):
            os.mkdir("sequences/sequence")

        all_files_plus = glob.glob(os.path.join("sequences/ssRNA_plus", "*.fasta"))
        all_files_minus = glob.glob(os.path.join("sequences/ssRNA_minus", "*.fasta"))

        extract(all_files_plus)
        extract(all_files_minus)

    if a.operation == "combine":

        all_files = glob.glob(os.path.join("sequences/genome", "*.fasta"))

        combine_files(all_files)


#main()

list_for_virus_names_plus = []
list_for_virus_names_minus = []
url_list = ["https://viralzone.expasy.org/245", "https://viralzone.expasy.org/240"]
virus_names_plus = crawler(url_list[0], list_for_virus_names_plus)
virus_names_minus = crawler(url_list[1], list_for_virus_names_minus)

accession_dict_plus = {}
accession_dict_minus = {}

unique_accession_plus = extracting_accession_numbers(virus_names_plus, accession_dict_plus)
unique_accession_minus = extracting_accession_numbers(virus_names_minus, accession_dict_minus)

foldername_minus = "genbank/viral/ssRNA_minus/"
orientation_minus = "(-)"
download_fasta(unique_accession_minus, foldername_minus, orientation_minus)

foldername_plus = "genbank/viral/ssRNA_plus/"
orientation_plus = "(+)"
download_fasta(unique_accession_plus, foldername_plus, orientation_plus)

all_files_plus = glob.glob(os.path.join("genbank/viral/ssRNA_plus", "*.fasta"))
all_files_minus = glob.glob(os.path.join("genbank/viral/ssRNA_minus", "*.fasta"))


def test(all_files):
    counter = 0
    for file in all_files:
        with open(file, "r") as f:
            if "genome" in f.readlines()[1]:
                counter += 1
    print(counter)


test(all_files_plus)
test(all_files_minus)

with open("sequences.fasta", "r") as file:
    counter = 0
    for line in file:
        if "genome" in line:
            counter += 1

print(counter)

with open("sequences_minus.fasta", "r") as file:
    counter = 0
    for line in file:
        if "genome" in line:
            counter += 1

print(counter)


