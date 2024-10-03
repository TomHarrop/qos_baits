#!/usr/bin/env python3

import os
import csv
from Bio import SeqIO


def read_merge_info(directory):
    merge_info = {}
    for file_name in os.listdir(directory):
        if file_name.endswith(".csv") and not file_name.startswith("."):
            with open(os.path.join(directory, file_name), mode="r") as file:
                reader = csv.DictReader(file)
                for row in reader:
                    if row["ref_hit"] not in merge_info:
                        merge_info[row["ref_hit"]] = []
                    merge_info[row["ref_hit"]].append(row["query_hit"])
    return merge_info



ref_targets = "output/010_captus/qos.mega353/min1000000/03_extractions/qos.1000000__captus-ext/01_coding_NUC/NUC_coding_NT.fna"
query_targets = "output/010_captus/qos.peakall/min1000000/03_extractions/qos.1000000__captus-ext/01_coding_NUC/NUC_coding_NT.fna"
loci_to_merge_dir = "output/020_overlaps/qos_min1000000.mega353.peakall/loci_to_merge"

combined_output = "combined.fasta"

merge_info = read_merge_info(loci_to_merge_dir)

file_handles = {}

# open a file handle for each locus that needs to be merged

# open a 

for 

