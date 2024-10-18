#!/usr/bin/env python3

from Bio import SeqIO
from pathlib import Path
from snakemake import logger
import csv
import logging


def read_merge_info(directory):
    merge_info = {}
    for file in Path(directory).glob("*.csv"):
        if file.name.endswith(".csv") and not file.name.startswith("."):
            with open(file, mode="r") as csv_file:
                reader = csv.DictReader(csv_file)
                for row in reader:
                    query_hit = row["query_hit"]
                    if "__" in query_hit:
                        query_hit = query_hit.split("__")[0]
                    merge_info[query_hit] = row["ref_hit"]
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


def rename_record(record, locus, records_by_locus, renamed_items):
    output_id = f"{record.id}-{locus}"
    logger.info(f"Record {record.id} will be added to locus {locus} as {output_id}.")
    renamed_items.append((record.id, output_id))
    if locus not in records_by_locus:
        records_by_locus[locus] = []
    record.id = output_id
    records_by_locus[locus].append(record)


# dev
# # inputs
# query_targets = "output/000_reference/query/peakall12.deduplicated_renamed.fasta"
# loci_to_merge_dirs = [
#     "output/020_overlaps/pzijinensis_min1000000.mega353.peakall12/loci_to_merge",
#     "output/020_overlaps/qos_min1000000.mega353.peakall12/loci_to_merge",
# ]
# orthofinder_dir = "output/005_grouped-targets/peakall12/orthofinder"

# # outputs
# renamed_targetfile = "test/renamed_original_targets.fasta"
# renamed_sequences = "test/renamed_original_targets.csv"


def main():
    all_merge_info = {}
    for loci_to_merge_dir in loci_to_merge_dirs:
        logger.info(f"Reading merge info from {loci_to_merge_dir}")
        merge_info = read_merge_info(loci_to_merge_dir)
        duplicate_keys = set(all_merge_info.keys()) & set(merge_info.keys())
        if duplicate_keys:
            mismatched_keys = {
                k for k in duplicate_keys if merge_info[k] != all_merge_info[k]
            }
            if mismatched_keys:
                raise ValueError(f"Sequence {mismatched_keys} matched different loci.")
        all_merge_info.update(merge_info)

    logger.info(
        f"Sequences overlapping loci {sorted(set(all_merge_info.values()))} were found. The sequences will be renamed."
    )

    logger.info(f"Reading orthogroup info from {orthofinder_dir}")
    orthogroup_file = Path(orthofinder_dir, "Orthogroups", "Orthogroups.txt")
    orthogroup_info, og_singletons = read_orthogroup_info(orthogroup_file)

    logger.info(f"Ignoring {len(og_singletons)} OrthoGroups with only one member.")

    # we will accumulate the records per-locus so they can be output in order
    records_by_locus = {}
    records_by_locus["__unassigned"] = []
    renamed_items = []
    h, i, j, k = 0, 0, 0, 0

    logger.info(f"Reading query targets from {query_targets}")
    query_records = SeqIO.parse(query_targets, "fasta")

    for record in query_records:
        # this already has locus information! ignore.
        if "-" in record.id:
            h += 1
            logger.info(f"Record {record.id} already has locus information.")
            locus = record.id.split("-")[-1]
            if locus not in records_by_locus:
                records_by_locus[locus] = []
            records_by_locus[locus].append(record)
        # these hits had overlaps with the mega353 hits
        elif record.id in all_merge_info:
            i += 1
            locus = all_merge_info[record.id]
            rename_record(record, locus, records_by_locus, renamed_items)
        # these hits appear to be orthologs, based on orthofinder results from the
        # peakall target file
        elif record.id in orthogroup_info:
            j += 1
            locus = orthogroup_info[record.id]
            rename_record(record, locus, records_by_locus, renamed_items)
        # neither overlaps or orthogroups
        else:
            k += 1
            logger.info(f"Record {record.id} won't be renamed.")
            records_by_locus["__unassigned"].append(record)

    # put records without locus info at the end
    records_by_locus["__unassigned"] = records_by_locus.pop("__unassigned")

    # what just happened
    logger.info(f"Found {h} records with existing locus information.")
    logger.info(f"Found {i} records to merge by overlap.")
    logger.info(f"Found {j} records to merge by orthogroup.")
    logger.info(f"Found {k} records to keep as-is.")

    # deduplicate the output
    for locus, records in records_by_locus.items():
        records_by_locus[locus] = list(
            {record.id: record for record in records}.values()
        )

    with open(renamed_targetfile, "w") as f:
        for locus, records in records_by_locus.items():
            SeqIO.write(records, f, "fasta")

    with open(renamed_sequences, "w") as f:
        writer = csv.writer(f)
        writer.writerow(["old_id", "new_id"])
        for old_id, new_id in renamed_items:
            writer.writerow([old_id, new_id])


if __name__ == "__main__":

    # inputs
    query_targets = snakemake.input["targetfile"]
    loci_to_merge_dirs = snakemake.input["loci_to_merge"]
    orthofinder_dir = snakemake.input["orthofinder"]

    # outputs
    renamed_targetfile = snakemake.output["renamed_targetfile"]
    renamed_sequences = snakemake.output["renamed_sequences"]

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
