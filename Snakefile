#!/usr/bin/env python3

from pathlib import Path


# containers
bbmap = "docker://quay.io/biocontainers/bbmap:39.01--h92535d8_1"
biopython = "docker://quay.io/biocontainers/biopython:1.78"
captus = "docker://quay.io/biocontainers/captus:1.0.1--pyhdfd78af_2"


# globals
outdir = Path("output")
logdir = Path(outdir, "logs")
qos_genome = Path("data", "reference", "QOS_assembly_hifi.fasta")
ccs_reads = Path("data", "raw_reads").glob("*.bam")
qos_genome_5k = Path(outdir, "000_reference", "assembly.5000.fasta")
flye_directory = Path(outdir, "020_flye")


#####################
# ASSEMBLY PIPELINE #
#####################

# currently flye is choking on this data


##########
# CAPTUS #
##########


# Try to extract the mega353 targets bundled with Captus.  Later we can do any
# targets by adding the target file and using
# `"--nuc_refs {input.target_file}"`
rule captus_extract:
    input:
        external_fasta=Path(
            outdir, "000_reference", "assembly.{minlength}.fasta"
        ),
    output:
        outdir=directory(
            Path(outdir, "040_captus", "min{minlength}", "03_extractions")
        ),
        refs_json=Path(
            "040_captus",
            "min{minlength}",
            "captus-assembly_extract.refs.json",
        ),
    log:
        Path(logdir, "extract.{minlength}.log"),
    benchmark:
        Path(logdir, "benchmark.extract.{minlength}.log")
    threads: lambda wildcards, attempt: 32
    resources:
        time=lambda wildcards, attempt: "5-00",
        mem_mb=lambda wildcards, attempt: 32e3,
    shadow:
        "minimal"
    container:
        captus
    shell:
        "captus_assembly extract "
        "--captus_assemblies_dir 02_assemblies "
        "--fastas {input.external_fasta} "
        "--out {output.outdir}/. "
        "--nuc_refs Mega353 "
        "--mit_refs SeedPlantsMIT "
        "--ptd_refs SeedPlantsPTD "
        '--ram "$(( {resources.mem_mb}/1000 ))" '
        "--threads {threads} "
        "--concurrent 3 "
        "&> {log} "
        "; mv {output.outdir}/captus-assembly_extract.refs.json "
        "{output.refs_json}"


rule captus_targets:
    default_target: True
    input:
        expand(
            [str(x) for x in rules.captus_extract.output],
            minlength=["1000000", "100000"],
        ),


############
# REF DATA #
############


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
