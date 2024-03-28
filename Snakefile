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
minimap2 = "docker://quay.io/biocontainers/minimap2:2.27--he4a0461_1"
samtools = "docker://quay.io/biocontainers/samtools:1.19--h50ea8bc_0"


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
ccs_reads = Path("data", "raw_reads").glob("*.bam")
qos_genome_5k = Path(outdir, "000_reference", "assembly.5000.fasta")


rule target:
    input:
        expand(
            Path(
                outdir, "010_repeatmasker.{assembly}", "genome_masked.fa.gz"
            ),
            assembly=["all", "100000"],
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


# map the CCS reads back and try to run purge haplotigs
rule map_target:
    input:
        expand(
            Path(outdir, "010_purge-haplotigs", "{minlength}", "aligned.bam"),
            minlength=["100000"],
        ),


rule sort_ccs_bamfile:
    input:
        Path(
            run_tmpdir,
            "010_purge-haplotigs",
            "{minlength}",
            "aligned.sam",
        ),
    output:
        Path(outdir, "010_purge-haplotigs", "{minlength}", "aligned.bam"),
    params:
        wd=Path(run_tmpdir, "010_purge-haplotigs"),
        mem_mb_per_thread=8e3,
    log:
        Path(logdir, "sort_ont_bamfile.{minlength}.log"),
    threads: 4
    resources:
        time=lambda wildcards, attempt: 480 * (2**attempt),
        mem_mb=lambda wildcards, threads: int(8e3) * threads,
    container:
        samtools
    shell:
        "samtools sort "
        "-m {params.mem_mb_per_thread}M "
        "-o {output} "
        "-T {params.wd} "
        "< {input} "
        "2> {log}"


rule map_ccs_reads:
    input:
        ref=Path(outdir, "000_reference", "assembly.{minlength}.fasta"),
        reads=Path(run_tmpdir, "ccs_reads.fastq"),
    output:
        pipe(
            Path(
                run_tmpdir,
                "010_purge-haplotigs",
                "{minlength}",
                "aligned.sam",
            )
        ),
    log:
        Path(logdir, "map_ont_reads.{minlength}.log"),
    threads: 24
    resources:
        time=lambda wildcards, attempt: 480 * (2**attempt),
        mem_mb=int(32e3),
    container:
        minimap2
    shell:
        "minimap2 "
        "-t {threads} "
        "-ax map-hifi "
        "{input.ref} "
        "{input.reads} "
        "--secondary=no "
        "> {output} "
        "2> {log}"


rule bam_to_fastq:
    input:
        ccs_reads,
    output:
        pipe(Path(run_tmpdir, "ccs_reads.fastq")),
    log:
        Path(logdir, "bam_to_fastq.log"),
    container:
        samtools
    shell:
        "samtools fastq {input} >> {output} 2> {log}"


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
