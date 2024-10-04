#!/usr/bin/env python3

from Bio import SeqIO
from pathlib import Path


def get_record_prefix(record):
    return record.id.split("_")[0]


# dev
# ungrouped_target_file = Path(
#     "output", "000_reference", "query", "peakall.deduplicated_renamed.fasta"
# )
# output_directory = Path("test", "grouped_targets")


def main(ungrouped_target_file, output_directory):
    input_handle = SeqIO.parse(ungrouped_target_file, "fasta")

    if not Path(output_directory).exists():
        Path(output_directory).mkdir(parents=True)

    output_handles = {}

    for record in input_handle:
        prefix = get_record_prefix(record)
        if prefix not in output_handles:
            output_handles[prefix] = open(
                Path(output_directory, f"{prefix}.fasta"), "w"
            )
        SeqIO.write(record, output_handles[prefix], "fasta")

    for handle in output_handles.values():
        handle.close()


if __name__ == "__main__":
    ungrouped_target_file = Path(snakemake.input["targets"])
    output_directory = Path(snakemake.output["outdir"])

    main(ungrouped_target_file, output_directory)
