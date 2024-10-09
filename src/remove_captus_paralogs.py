#!/usr/bin/env python3

from Bio import SeqIO
from pathlib import Path
from snakemake import logger
import logging
import re


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
    loci_with_paralogs = {}

    input_target_handle = SeqIO.parse(targets_with_paralogs, "fasta")

    for record in input_target_handle:
        seq_name, locus = split_record_id(record.id)
        if "paralog" in seq_name:
            paralog_number = re.search(r"paralog(\d+)", seq_name).group(1)
            if locus not in loci_with_paralogs:
                loci_with_paralogs[locus] = {}
            loci_with_paralogs[locus][paralog_number] = record
        else:
            if locus in records_by_locus:
                logger.error(f"Duplicate locus: {locus}")
                logger.error(record)
                raise KeyError(f"Duplicate locus: {locus}")
            records_by_locus[locus] = record

    for locus, paralog_dict in loci_with_paralogs.items():
        paralog_numbers = list(paralog_dict.keys())
        lowest_paralog = min(paralog_numbers)
        logger.info(f"Choosing paralog {lowest_paralog} for locus {locus}")
        chosen_record = paralog_dict[lowest_paralog]
        logger.info(chosen_record.id)
        seq_name, locus = split_record_id(chosen_record.id)
        if locus in records_by_locus:
            logger.error(f"Duplicate locus: {locus}")
            logger.error(record)
            raise KeyError(f"Duplicate locus: {locus}")
        records_by_locus[locus] = chosen_record

    with open(output_target, "w") as f:
        for locus, records in records_by_locus.items():
            SeqIO.write(records, f, "fasta")


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
