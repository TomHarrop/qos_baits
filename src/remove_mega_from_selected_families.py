#!/usr/bin/env python3

from Bio import SeqIO
from pathlib import Path
from snakemake import logger
import logging


def main():
    with open(select_file, mode="r") as file:
        select_families = [
            line
            for line in file.read().splitlines()
            if line and not line.startswith("#") and not line.startswith("[")
        ]

    logger.info(f"Keeping records that match any of {select_families}")

    records = SeqIO.parse(input_targets, format="fasta")

    with open(output_targets, mode="w") as file:
        for record in records:
            if record.id.split("-")[0] in select_families:
                logger.info(f"Selecting {record.id}")
                SeqIO.write(record, file, format="fasta")
            else:
                logger.debug(f"Discarding {record.id}")


if __name__ == "__main__":
    select_file = snakemake.input["select_file"]
    input_targets = snakemake.input["targets"]

    output_targets = snakemake.output["targets"]

    logfile = snakemake.log[0]

    file_handler = logging.FileHandler(logfile)
    logger.logfile_handler = file_handler
    logger.logger.addHandler(logger.logfile_handler)

    try:
        main()
    except Exception as e:
        logger.error(e)
        raise e
