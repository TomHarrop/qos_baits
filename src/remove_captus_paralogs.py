#!/usr/bin/env python3

from Bio import SeqIO
from pathlib import Path
from snakemake import logger
import logging


def append_to_locus(record, locus, locus_dict):
    if locus not in locus_dict:
        locus_dict[locus] = []
    locus_dict[locus].append(record)


def split_record_id(record_id):
    return record_id.split("-")


# dev
# targets_with_paralogs = Path(
#     "output/030_merged-target-sequences/pzijinensis_min1000000.mega353.peakall/merged_targets.fasta"
# )
# logfile = "test.log"


def main():
    # we will accumulate the records per-locus so they can be output in order
    records_by_locus = {}
    discarded_paralogs = {}

    logger.info(f"Reading input from {targets_with_paralogs}")
    input_target_handle = SeqIO.parse(targets_with_paralogs, "fasta")

    i, j, k = 0, 0, 0

    for record in input_target_handle:
        seq_name, locus = split_record_id(record.id)
        if "paralog" not in seq_name:
            i += 1
            append_to_locus(record, locus, records_by_locus)
        # we're only trying to resolve sequences that Captus *knows* are
        # paralogs
        elif "paralog" in seq_name:
            if "[hit=00]" in record.description:
                j += 1
                append_to_locus(record, locus, records_by_locus)
            else:
                k += 1
                append_to_locus(record, locus, discarded_paralogs)
        else:
            logger.error("What is this place?")
            raise ValueError(f"Can't do anything with {record.id}")

    logger.info(f"Outputting records for {len(records_by_locus)} loci.")
    logger.info(
        f"{i} records were not identified as paralogs by Captus. These are kept as-is."
    )
    logger.info(
        f'{j} records were identified as paralogs by Captus and are kept because they had "[hit=00]" in the description.'
    )
    logger.info(
        f"{k} records with lower hit scores were discarded for {len(discarded_paralogs)} loci."
    )

    with open(output_target, "w") as f:
        for locus, record_list in records_by_locus.items():
            
            for record in record_list:
                SeqIO.write(record, f, "fasta")


if __name__ == "__main__":

    targets_with_paralogs = Path(snakemake.input["targets"])
    output_target = Path(snakemake.output["targets"])
    logfile = Path(snakemake.log[0])

    file_handler = logging.FileHandler(logfile)
    logger.logfile_handler = file_handler
    logger.logger.addHandler(logger.logfile_handler)

    try:
        main()
    except Exception as e:
        logger.error(e)
        raise e
