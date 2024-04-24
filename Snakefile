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
flye = "docker://quay.io/biocontainers/flye:2.9.3--py310h2b6aa90_1"
minimap2 = "docker://quay.io/biocontainers/minimap2:2.27--he4a0461_1"
samtools = "docker://quay.io/biocontainers/samtools:1.19--h50ea8bc_0"

# modules
module_tag = "0.0.50"
rm_snakefile = github(
    "tomharrop/smk-modules",
    path="modules/repeatmasker/Snakefile",
    tag=module_tag,
)
purge_snakefile = github(
    "tomharrop/smk-modules",
    path="modules/purge_haplotigs/Snakefile",
    tag=module_tag,
)

# globals
outdir = Path("output")
logdir = Path(outdir, "logs")
qos_genome = Path("data", "reference", "QOS_assembly_hifi.fasta")
ccs_reads = Path("data", "raw_reads").glob("*.bam")
qos_genome_5k = Path(outdir, "000_reference", "assembly.5000.fasta")
flye_directory = Path(outdir, "020_flye")


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
            Path(
                outdir, "010_purge-haplotigs", "{minlength}", "histogram.png"
            ),
            minlength=["100000"],
        ),


module purge_haplotigs:
    snakefile:
        purge_snakefile
    config:
        {
            "bamfile": Path(
            outdir, "010_purge-haplotigs", "{minlength}", "aligned.bam"
            ),
            "contigs": Path(
            outdir, "000_reference", "assembly.{minlength}.fasta"
            ),
            "outdir": Path(outdir, "010_purge-haplotigs", "{minlength}"),
            "run_tmpdir": Path(
            run_tmpdir, "010_purge-haplotigs", "{minlength}"
            ),
        }


use rule * from purge_haplotigs as purge_haplotigs_*


use rule purge from purge_haplotigs as purge_haplotigs_purge with:
    threads: lambda wildcards, attempt: 48 * attempt
    resources:
        time=lambda wildcards, attempt: f"{attempt - 1}-23",
        mem_mb=lambda wildcards, attempt: int(256e3),


use rule hist from purge_haplotigs as purge_haplotigs_hist with:
    threads: lambda wildcards, attempt: 12 * attempt
    resources:
        time=lambda wildcards, attempt: 90 * attempt,


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
        mem_mb_per_thread=int(8e3),
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


rule flye:
    input:
        ccs=Path(run_tmpdir, "ccs_reads.fastq"),
    output:
        Path(flye_directory, "assembly.fasta"),
    params:
        outdir=lambda wildcards, input: Path(input.ccs[0]).parent,
    log:
        Path(logdir, "flye.log"),
    threads: min(130, workflow.cores) - 2
    resources:
        time=167 * 60,  # 6 days 23 hours
        mem_mb=int(384e3),
    container:
        flye
    shell:
        "flye "
        "--pacbio-hifi "
        "{input.ccs} "
        "--out-dir {params.outdir} "
        "--threads {threads} "
        "&>> {log}"



rule bam_to_fastq:
    input:
        Path(run_tmpdir, "ccs_reads.bam"),
    output:
        temp(Path(run_tmpdir, "ccs_reads.fastq")),
    log:
        Path(logdir, "bam_to_fastq.log"),
    threads: 1
    container:
        samtools
    shell:
        "samtools fastq "
        "-o {output} "
        "&>{log} "
        "<{input} "

rule samtools_cat:
    input:
        ccs_reads,
    output:
        pipe(Path(run_tmpdir, "ccs_reads.bam")),
    log:
        Path(logdir, "samtools_cat.log"),
    threads: 1
    container:
        samtools
    shell:
        "samtools cat "
        "{input} "
        "2>{log} "
        ">> {output} "


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
