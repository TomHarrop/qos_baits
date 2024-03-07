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
bbmap = "docker://quay.io/biocontainers/bbmap:39.01--h92535d8_1"

# modules
module_tag = "0.0.44"
rm_snakefile = github(
    "tomharrop/smk-modules",
    path="modules/repeatmasker/Snakefile",
    tag=module_tag,
)

# globals
outdir = Path("output")
logdir = Path(outdir, "logs")
qos_genome = Path("data", "reference", "QOS_assembly_hifi.fasta")
qos_genome_5k = Path(outdir, "000_reference", "assembly.5000.fasta")


rule target:
    input:
        expand(
            Path(
                outdir, "010_repeatmasker.{assembly}", "genome_masked.fa.gz"
            ),
            assembly=["all", "5000"],
        ),


# keep this module separate
module rm_all:
    snakefile:
        rm_snakefile
    config:
        {
            "outdir": Path(outdir, "010_repeatmasker.all"),
            "query_genome": qos_genome,
            "rm_output": Path(
            outdir, "010_repeatmasker.all", "genome_masked.fa.gz"
            ),
            "run_tmpdir": run_tmpdir,
        }


use rule * from rm_all as rm_all_*


module rm_subset:
    snakefile:
        rm_snakefile
    config:
        {
            "outdir": Path(outdir, "010_repeatmasker.{assembly}"),
            "query_genome": Path(
            outdir, "000_reference", "assembly.{assembly}.fasta"
            ),
            "rm_output": Path(
            outdir, "010_repeatmasker.{assembly}", "genome_masked.fa.gz"
            ),
            "run_tmpdir": run_tmpdir,
        }


use rule * from rm_subset as rm_subset_*


# this genome is highly fragmented
rule reformat:
    input:
        qos_genome,
    output:
        Path(outdir, "000_reference", "assembly.{minlength}.fasta"),
    log:
        Path(logdir, "reformat.{minlength}.log"),
    container:
        bbmap
    shell:
        "reformat.sh "
        "in={input} "
        "minlength={wildcards.minlength} "
        "out={output} "
        "2>{log}"
