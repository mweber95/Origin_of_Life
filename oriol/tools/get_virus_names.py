import re
import requests
from bs4 import BeautifulSoup as bs


class ViralZone:
    def __init__(self, id_type):
        self.id_type = id_type

    def crawler(self, url_list):
        baltimore_groups = []
        for url in url_list:
            html_code = requests.get(url_list[url])
            soup = bs(html_code.content, "html.parser")
            c = [viruses for viruses in soup.getText().split() if any(x in viruses for x in ["virus", "viridae"])]
            baltimore_groups.append(c)
        return baltimore_groups

    def extracting_accession_numbers(self, virus_names):
        accession_dict = {}
        with open("oriol/data/viral.txt") as file:
            for line in file:                           ### for line, virus in zip()
                for virus in virus_names:
                    if virus in line:
                        if self.id_type == "ncbi":
                            ncbi_id = line.strip()[:9]
                            accession_dict[ncbi_id] = virus
                        elif self.id_type == "genbank":
                            genbank_id = re.split(r'\t+', line.strip())[1]
                            accession_dict[genbank_id] = virus
        return accession_dict
