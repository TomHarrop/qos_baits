#!/usr/bin/env python3

from Bio import SeqIO
from functools import cache
from pathlib import Path
from snakemake import logger
import csv
import logging
import os


def generate_new_spec(record_id, species_name):
    new_id, locus = get_new_id(record_id, species_name)
    new_spec = new_id.replace("-", "_")
    return new_spec


@cache
def get_new_id(record_id, species_name):
    try:
        old_spec, locus, paralog_number = record_id.split("__")
        new_id = f"{species_name}_paralog{paralog_number}-{locus}"
    except ValueError:
        old_spec, locus = record_id.split("__")
        new_id = f"{species_name}-{locus}"
    return new_id, locus


def read_merge_info(directory):
    merge_info = {}
    for file_name in os.listdir(directory):
        if file_name.endswith(".csv") and not file_name.startswith("."):
            with open(os.path.join(directory, file_name), mode="r") as file:
                reader = csv.DictReader(file)
                for row in reader:
                    merge_info[row["query_hit"]] = row["ref_hit"]
    return merge_info


def read_orthogroup_info(orthogroup_file):
    seq_id_to_orthogroup = {}
    not_an_orthogroup = []
    with open(orthogroup_file, mode="r") as file:
        for line in file:
            og, seq_ids = line.strip().split(":")
            split_ids = seq_ids.split()
            # it's not a "group" if there's only one sequence in it
            if len(split_ids) > 1:
                for seq_id in split_ids:
                    seq_id_to_orthogroup[seq_id] = og
            else:
                not_an_orthogroup.append(split_ids[0])
    return seq_id_to_orthogroup, not_an_orthogroup


# dev
# inputs
# species_name = "qos"
# ref_targets = "output/010_captus/qos.mega353/min1000000/03_extractions/qos.1000000__captus-ext/01_coding_NUC/NUC_coding_NT.fna"
# query_targets = "output/010_captus/qos.peakall/min1000000/03_extractions/qos.1000000__captus-ext/01_coding_NUC/NUC_coding_NT.fna"
# loci_to_merge_dir = "output/020_overlaps/qos_min1000000.mega353.peakall/loci_to_merge"
# orthofinder_dir = "output/005_grouped-targets/peakall/orthofinder"

# outputs
# combined_output = "test/merged_targets.fasta"
# renamed_sequences = "test/renamed_sequences.csv"


def main():
    # This is a dict of sequence ID to locus to merge to. Any sequence that is
    # recorded as a "locus to merge" will be renamed according to the locus it
    # is to be merged with
    logger.info(f"Reading merge info from {loci_to_merge_dir}")
    merge_info = read_merge_info(loci_to_merge_dir)
    logger.info(f"Loci matching {sorted(set(merge_info.values()))} will be merged.")

    # Any sequence that is in an orthogroup will be renamed according to that
    # orthogroup. The orthogroup text file also contains singletons.
    logger.info(f"Reading orthogroup info from {orthofinder_dir}")
    orthogroup_file = Path(orthofinder_dir, "Orthogroups", "Orthogroups.txt")
    orthogroup_info, og_singletons = read_orthogroup_info(orthogroup_file)

    logger.info(f"Ignoring {len(og_singletons)} OrthoGroups with only one member.")

    # we will accumulate the records per-locus so they can be output in order
    records_by_locus = {}

    # start with the reference loci
    logger.info(f"Reading reference targets")
    ref_records = SeqIO.parse(ref_targets, "fasta")

    i = 0
    for record in ref_records:
        i += 1
        new_id, locus = get_new_id(record.id, species_name)
        record.id = new_id
        if locus not in records_by_locus:
            records_by_locus[locus] = []
        records_by_locus[locus].append(record)

    logger.info(
        f"{ref_targets} contains {i} reference loci in {len(records_by_locus)} loci."
    )

    # sort the query loci
    records_that_will_be_merged_by_locus = {}
    records_that_will_be_merged_by_orthogroup = {}
    i, j, k = 0, 0, 0

    query_records = SeqIO.parse(query_targets, "fasta")
    for record in query_records:
        new_id, locus = get_new_id(record.id, species_name)
        if locus in merge_info:
            i += 1
            # these hits had overlaps with the mega353 hits
            logger.info(f"Locus {locus} overlaps with a mega353 locus.")
            if locus not in records_that_will_be_merged_by_locus:
                records_that_will_be_merged_by_locus[locus] = []
            records_that_will_be_merged_by_locus[locus].append(record)
        elif locus in orthogroup_info:
            j += 1
            logger.info(f"Locus {locus} was placed in an orthogroup.")
            # these hits appear to be orthologs, based on orthofinder results from
            # the peakall target file
            if locus not in records_that_will_be_merged_by_orthogroup:
                records_that_will_be_merged_by_orthogroup[locus] = []
            records_that_will_be_merged_by_orthogroup[locus].append(record)
        else:
            # these had hits in the query, didn't overlap mega353 hits, and
            # aren't orthlogs, so we will keep the current locus name.
            k += 1
            logger.info(f"Locus {locus} won't be merged.")
            record.id = new_id
            if locus not in records_by_locus:
                records_by_locus[locus] = []
            records_by_locus[locus].append(record)

    # what just happened
    logger.info(f"Found {i} loci to merge by overlap.")
    logger.info(f"Found {j} loci to merge by orthogroup.")
    logger.info(f"Found {k} loci to keep as-is.")

    # keep track of what we rename
    renamed_items = []

    # merge the records that had overlaps with the mega353 targets
    for locus, records in records_that_will_be_merged_by_locus.items():
        new_locus = merge_info[locus]
        for record in records:
            # we are going to retain the original locus but move it to the
            # "species" field, because multiple targets from peakall can overlap
            # with the same mega353 target and they need to be identified
            # separately.
            new_spec = generate_new_spec(record.id, species_name)
            output_id = f"{new_spec}-{new_locus}"
            renamed_items.append((record.id, output_id))
            record.id = output_id
            if new_locus not in records_by_locus:
                records_by_locus[new_locus] = []
            records_by_locus[new_locus].append(record)

    # merge the records that were in orthogroups
    for locus, records in records_that_will_be_merged_by_orthogroup.items():
        new_locus = orthogroup_info[locus]
        for record in records:
            # we are going to retain the original locus but move it to the
            # "species" field to keep track of which target from peakall generated
            # the hit.
            new_spec = generate_new_spec(record.id, species_name)
            output_id = f"{new_spec}-{new_locus}"
            renamed_items.append((record.id, output_id))
            record.id = output_id
            if new_locus not in records_by_locus:
                records_by_locus[new_locus] = []
            records_by_locus[new_locus].append(record)

    # deduplicate the output
    for locus, records in records_by_locus.items():
        records_by_locus[locus] = list(
            {record.id: record for record in records}.values()
        )

    with open(combined_output, "w") as f:
        for locus, records in records_by_locus.items():
            SeqIO.write(records, f, "fasta")

    with open(renamed_sequences, "w") as f:
        writer = csv.writer(f)
        writer.writerow(["old_id", "new_id"])
        for old_id, new_id in renamed_items:
            writer.writerow([old_id, new_id])


if __name__ == "__main__":

    # inputs
    ref_targets = snakemake.input["ref_targets"]
    query_targets = snakemake.input["query_targets"]
    loci_to_merge_dir = snakemake.input["loci_to_merge"]
    orthofinder_dir = snakemake.input["orthofinder"]

    # outputs
    combined_output = snakemake.output["merged_targets"]
    renamed_sequences = snakemake.output["renamed_sequences"]

    # params
    species_name = snakemake.wildcards["ref_dataset"]

    # log
    logfile = snakemake.log[0]

    file_handler = logging.FileHandler(logfile)
    logger.logfile_handler = file_handler
    logger.logger.addHandler(logger.logfile_handler)

    try:
        main()
    except Exception as e:
        logger.error(e)
        raise e
