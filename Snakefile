#!/usr/bin/env python3

from pathlib import Path
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

# containers
bbmap = "docker://quay.io/biocontainers/bbmap:39.01--h92535d8_1"
biopython = "docker://quay.io/biocontainers/biopython:1.78"
captus = "docker://quay.io/biocontainers/captus:1.0.1--pyhdfd78af_2"


# globals
outdir = Path("output")
logdir = Path(outdir, "logs")
qos_genome = Path("data", "reference", "QOS_assembly_hifi.fasta")

mega353_file_url = (
    "https://github.com/chrisjackson-pellicle/"
    "NewTargets/raw/v1.0.0/mega353.fasta.zip"
)

target_files = [
    Path(outdir, "000_reference", "mega353.fasta.gz"),
    # This is Lalita's target file, i *think* this is what lab members refer to
    # as "Esserman"
    Path("data", "reference", "Orchidaceae963TargetFile_LSv1.fasta.gz"),
]

# snakemake's downloader
HTTP = HTTPRemoteProvider()


#####################
# ASSEMBLY PIPELINE #
#####################

# currently flye is choking on this data


##########
# CAPTUS #
##########


def format_extract_target(wildcards, input):
    if wildcards.targetset == "mega353":
        return "--nuc_refs Mega353"
    elif wildcards.targetset == "mega353_plus_orchidaceae":
        target_path = input.targets
        return f"--nuc-refs {target_path}"
    else:
        raise valueError(f"wtf {wildcards}")


def get_extract_input(wildcards):
    input_dict = {
        "external_fasta": Path(
            outdir, "000_reference", "assembly.{minlength}.fasta"
        ),
    }
    if wildcards.targetset == "mega353_plus_orchidaceae":
        input_dict["targets"] = Path(
            outdir,
            "000_reference",
            "combined_targetfiles.deduplicated_renamed.fasta",
        )
    return input_dict


# Try to extract the mega353 targets bundled with Captus.  Later we can do any
# targets by adding the target file and using
# `"--nuc_refs {input.target_file}"`
rule captus_extract:
    input:
        unpack(get_extract_input),
    output:
        outdir=directory(
            Path(
                outdir,
                "040_captus",
                "min{minlength}.{targetset}",
                "03_extractions",
            )
        ),
        refs_json=Path(
            "040_captus",
            "min{minlength}.{targetset}",
            "captus-assembly_extract.refs.json",
        ),
    params:
        targets=format_extract_target,
    log:
        Path(logdir, "extract.{minlength}.{targetset}.log"),
    benchmark:
        Path(logdir, "benchmark.extract.{minlength}.{targetset}.log")
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
        "{params.targets} "
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
            targetset=["mega353", "mega353_plus_orchidaceae"],
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


rule rename_target_sequences:
    input:
        Path(
            outdir, "000_reference", "combined_targetfiles.deduplicated.fasta"
        ),
    output:
        Path(
            outdir,
            "000_reference",
            "combined_targetfiles.deduplicated_renamed.fasta",
        ),
    log:
        Path(logdir, "rename_target_sequences.log"),
    resources:
        time=1,
    container:
        bbmap
    shell:
        "reformat.sh "
        "in=stdin.fasta "
        "out={output} "
        "trimreaddescription=t "
        "fixheaders=t "
        "fixjunk=t "
        "< {input} "
        "2> {log}"


rule combine_targetfile:
    input:
        target_files,
    output:
        pipe=pipe(
            Path(
                outdir,
                "000_reference",
                "combined_targetfiles.deduplicated.fasta",
            )
        ),
        duplicates=Path(
            outdir,
            "000_reference",
            "combined_targetfiles.discarded_duplicates.fasta",
        ),
    log:
        Path(logdir, "combine_targetfile.log"),
    threads: 2
    resources:
        time=1,
    container:
        bbmap
    shell:
        "cat {input} | "
        "dedupe.sh "
        "in=stdin.fasta.gz "
        "out=stdout.fasta "
        "outd={output.duplicates} "
        "uniquenames=t "
        "sort=name "
        "ascending=t "
        "exact=t "
        "touppercase=t "
        "maxsubs=0 "
        "maxedits=0 "
        "> {output.pipe} "
        "2> {log} "


rule download_targetfile:
    input:
        HTTP.remote(mega353_file_url, keep_local=True),
    output:
        Path(outdir, "000_reference", "mega353.fasta.gz"),
    threads: 1
    resources:
        partition_flag="--partition=io",
    shell:
        "unzip -p {input} | gzip -1 >> {output}"
