#!/usr/bin/env python3

from pathlib import Path
import tempfile


# set up a temporary directory for this run
try:
    run_tmpdir = config["run_tmpdir"]
    print(f"Caught run_tmpdir {run_tmpdir}")
except KeyError as e:
    print(f"{e} not set in config")
    run_tmpdir = tempfile.mkdtemp()
    print(f"Setting run_tmpdir to {run_tmpdir}")
    print("This probably won't work on a cluster!")

# containers

# modules
module_tag = "0.0.44"
rm_snakefile = github(
    "tomharrop/smk-modules",
    path="modules/repeatmasker/Snakefile",
    tag=module_tag,
)

# globals
outdir = Path("output")
qos_genome = Path("data", "reference", "QOS_assembly_hifi.fasta")


module repeatmasker:
    snakefile:
        rm_snakefile
    config:
        {
            "outdir": Path(outdir, "010_repeatmasker"),
            "query_genome": qos_genome,
            "rm_output": Path(
            outdir, "010_repeatmasker", "genome_masked.fa.gz"
            ),
            "run_tmpdir": run_tmpdir,
        }


use rule * from repeatmasker as rm_*
