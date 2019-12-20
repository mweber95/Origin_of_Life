import os
import re
import requests
import numpy as np
from bs4 import BeautifulSoup as bs
from sklearn.datasets import load_digits
from sklearn.manifold import MDS, TSNE
import matplotlib.pyplot as plt
import umap



def crawler(url, list_for_virus_names):
    html_code = requests.get(url)
    soup = bs(html_code.content, "html.parser")
    for tag in soup.find_all("a", {"href": re.compile(r"/by_species/[0-9]")})[7:]:
        list_for_virus_names.append(tag.getText())
    return list_for_virus_names


def extracting_accession_numbers(virus_names, accession_dict, plus_minus):
    with open("sequences/processing/viral.txt") as file:
        for line in file:
            for virus in virus_names:
                if virus in line:
                    ncbi_id = line.strip()[:9]
                    accession_dict[ncbi_id] = (virus, plus_minus)
    return accession_dict


def distmat_alt(plus, minus):
    with open("sequences/genome/distmat_genome.mat", "r+") as file:
        all_lines = file.readlines()
        with open("sequences/genome/distmat_altered.mat", "w") as outfile:
            for line in all_lines[1:]:
                if line[:9] in plus:
                    value, sign = plus.get(line[:9])
                elif line[:9] in minus:
                    value, sign = minus.get(line[:9])
                else:
                    continue
                outfile.write(f"{line[:9]} {value} ({sign}) {line[10:]}")


def distmat_2():
    with open("sequences/genome/distmat_genome.mat", "r+") as file:
        all_lines = file.readlines()
        with open("sequences/genome/distmat_altered_2.mat", "w") as outfile:
            for line in all_lines[1:]:
                outfile.write(line[10:])


def reading_distmat():
    array = np.genfromtxt("sequences/genome/distmat_altered_2.mat", skip_header=0)
    #embedding = umap.UMAP(n_neighbors=10, min_dist=0.1, metric='euclidean')
    embedding = TSNE(n_components=2, perplexity=40.0, n_iter=1000)
    coord = embedding.fit_transform(array)
    return coord


def main():

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

    plus = "+"
    minus = "-"
    unique_accession_plus = extracting_accession_numbers(virus_names_plus, accession_dict_plus, plus)
    unique_accession_minus = extracting_accession_numbers(virus_names_minus, accession_dict_minus, minus)

    print(len(unique_accession_plus))
    print(len(unique_accession_minus))

    #distmat_alt(unique_accession_plus, unique_accession_minus)
    #distmat_2()
    coord = reading_distmat()
    print(coord)
    print(coord.shape)
    np.savetxt("ohne_label_TSNE_perp40.csv", coord, delimiter=",")

    #plt.scatter(coord[0][:], coord[1][:])
    #plt.show()

main()