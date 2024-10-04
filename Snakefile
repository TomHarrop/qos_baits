#!/usr/bin/env python3

from pathlib import Path
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider


def clobber_type(file):
    if isinstance(file, list):
        return file[0]
    else:
        return file


def generate_conversion_code(wildcards, input, output):
    my_string = clobber_type(input)
    my_suffix = Path(my_string).suffix
    if my_suffix == ".zip":
        return f"unzip -p {input} | gzip -1 >> {output}"
    elif my_suffix == ".gz":
        return f"cp {input} {output}"
    elif my_suffix in [".fasta", ".fa", ".fna"]:
        return f"gzip -1 < {input} > {output}"
    else:
        return f"cp {input} > {output}"


def get_remote_file(wildcards):
    if wildcards.datatype == "query":
        my_file = query_file_locations[wildcards.dataset]
        return my_file
    elif wildcards.datatype == "reference":
        my_file = reference_genome_locations[wildcards.dataset]
        return my_file
    else:
        raise ValueError(f"Unknown datatype {wildcards.datatype}")


# containers
bbmap = "docker://quay.io/biocontainers/bbmap:39.01--h92535d8_1"
biopython = "docker://quay.io/biocontainers/biopython:1.78"
captus = "docker://quay.io/biocontainers/captus:1.0.1--pyhdfd78af_2"
orthofinder = "docker://davidemms/orthofinder:2.5.5.2"
r = "docker://ghcr.io/tomharrop/r-containers:r2u_24.04_cv1"
samtools = "docker://quay.io/biocontainers/samtools:1.21--h50ea8bc_0"

# globals
outdir = Path("output")
logdir = Path(outdir, "logs")

# snakemake's downloader
HTTP = HTTPRemoteProvider()

# input files
reference_genome_locations = {
    "qos": Path("data", "reference", "QOS_assembly_hifi.fasta"),
    "pzijinensis": HTTP.remote(
        (
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/"
            "GCA/039/513/925/GCA_039513925.1_PZIJ_v1.0/"
            "GCA_039513925.1_PZIJ_v1.0_genomic.fna.gz"
        ),
        keep_local=True,
    ),
    "pguangdongensis": HTTP.remote(
        (
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/"
            "GCA/039/583/875/GCA_039583875.1_PGUA_v1.0/"
            "GCA_039583875.1_PGUA_v1.0_genomic.fna.gz"
        ),
        keep_local=True,
    ),
}
all_reference_genomes = reference_genome_locations.keys()

query_file_locations = {
    "mega353": HTTP.remote(
        (
            "https://github.com/chrisjackson-pellicle/"
            "NewTargets/raw/v1.0.0/mega353.fasta.zip"
        ),
        keep_local=True,
    ),
    # This is Lalita's target file, i *think* this is what lab members refer to
    # as "Esserman"
    "orchidaceae963": Path(
        "data", "reference", "Orchidaceae963TargetFile_LSv1.fasta.gz"
    ),
    # Downloaded from
    # https://datadryad.org/stash/dataset/doi:10.5061/dryad.z08kprrbj
    "peakall": Path("data", "reference", "S1-2FinalSeq.fa"),
    # Downloaded from
    # https://datadryad.org/stash/dataset/doi:10.5061/dryad.sj3tx96bn
    "orchidinae205": Path(
        "data", "reference", "orchidinae-205_reference.fasta"
    ),
}

all_query_datasets = query_file_locations.keys()


wildcard_constraints:
    ref_targets="|".join(all_query_datasets),
    query_targets="|".join(all_query_datasets),
    ref_dataset="|".join(all_reference_genomes),
    query_dataset="|".join(all_query_datasets),


#####################
# ASSEMBLY PIPELINE #
#####################

# currently flye is choking on this data


#########################
# FIND OVERLAPPING LOCI #
#########################


rule find_overlaps_target:
    input:
        # includes check for overlaps among reference loci
        expand(
            Path(
                outdir,
                "020_overlaps",
                "{ref_dataset}_min{minlength}.{ref_targets}.{query_targets}",
                "loci_to_merge",
            ),
            ref_dataset=all_reference_genomes,
            minlength=["1000000"],
            ref_targets=["mega353"],
            query_targets=[x for x in all_query_datasets if x != "mega353"],
        ),
        expand(
            Path(
                outdir,
                "030_merged-target-sequences",
                "{ref_dataset}_min{minlength}.{ref_targets}.{query_targets}",
                "merged_targets.fasta",
            ),
            ref_dataset=["qos", "pzijinensis"],
            minlength=["1000000"],
            ref_targets=["mega353"],
            query_targets=["peakall"],
        ),


# only implemented for peakall vs mega353
rule merge_extracted_sequences:
    input:
        ref_targets=Path(
            outdir,
            "010_captus",
            "{ref_dataset}.{ref_targets}",
            "min{minlength}",
            "03_extractions",
            "{ref_dataset}.{minlength}__captus-ext",
            "01_coding_NUC",
            "NUC_coding_NT.fna",
        ),
        query_targets=Path(
            outdir,
            "010_captus",
            "{ref_dataset}.{ref_targets}",
            "min{minlength}",
            "03_extractions",
            "{ref_dataset}.{minlength}__captus-ext",
            "01_coding_NUC",
            "NUC_coding_NT.fna",
        ),
        loci_to_merge=Path(
            outdir,
            "020_overlaps",
            "{ref_dataset}_min{minlength}.{ref_targets}.{query_targets}",
            "loci_to_merge",
        ),
        orthofinder=Path(
            outdir,
            "005_grouped-targets",
            "{query_targets}",
            "orthofinder",
        ),
    output:
        merged_targets=Path(
            outdir,
            "030_merged-target-sequences",
            "{ref_dataset}_min{minlength}.{ref_targets}.{query_targets}",
            "merged_targets.fasta",
        ),
        renamed_sequences=Path(
            outdir,
            "030_merged-target-sequences",
            "{ref_dataset}_min{minlength}.{ref_targets}.{query_targets}",
            "renamed_sequences.csv",
        ),
    log:
        Path(
            logdir,
            "merge_extracted_sequences",
            "{ref_dataset}_min{minlength}.{ref_targets}.{query_targets}.log",
        ),
    container:
        biopython
    script:
        "src/merge_extracted_sequences.py"


rule generate_lists_of_loci_to_merge:
    input:
        overlapping_loci=Path(
            outdir,
            "020_overlaps",
            "{ref_dataset}_min{minlength}.{ref_targets}.{query_targets}",
            "overlapping_loci.csv",
        ),
        ref_self_overlaps=Path(
            outdir,
            "020_overlaps",
            "{ref_dataset}_min{minlength}.{ref_targets}.{ref_targets}",
            "overlapping_loci.csv",
        ),
    output:
        outdir=directory(
            Path(
                outdir,
                "020_overlaps",
                "{ref_dataset}_min{minlength}.{ref_targets}.{query_targets}",
                "loci_to_merge",
            )
        ),
    log:
        Path(
            logdir,
            "generate_lists_of_loci_to_merge",
            "{ref_dataset}_min{minlength}.{ref_targets}.{query_targets}.log",
        ),
    container:
        r
    script:
        "src/generate_lists_of_loci_to_merge.R"


rule find_overlapping_loci:
    input:
        query_gff=Path(
            outdir,
            "010_captus",
            "{ref_dataset}.{query_targets}",
            "min{minlength}",
            "03_extractions",
            "{ref_dataset}.{minlength}__captus-ext",
            "06_assembly_annotated",
            "{ref_dataset}.{minlength}_hit_contigs.gff",
        ),
        ref_gff=Path(
            outdir,
            "010_captus",
            "{ref_dataset}.{ref_targets}",
            "min{minlength}",
            "03_extractions",
            "{ref_dataset}.{minlength}__captus-ext",
            "06_assembly_annotated",
            "{ref_dataset}.{minlength}_hit_contigs.gff",
        ),
        fai=Path(
            outdir, "000_reference", "reference", "{ref_dataset}.fasta.fai"
        ),
    output:
        overlapping_loci=Path(
            outdir,
            "020_overlaps",
            "{ref_dataset}_min{minlength}.{ref_targets}.{query_targets}",
            "overlapping_loci.csv",
        ),
        proximate_loci=Path(
            outdir,
            "020_overlaps",
            "{ref_dataset}_min{minlength}.{ref_targets}.{query_targets}",
            "proximate_loci.csv",
        ),
    params:
        maxgap=int(10000),  # maxgap between "proximate" loci
    log:
        Path(
            logdir,
            "find_overlapping_loci",
            "{ref_dataset}_min{minlength}.{ref_targets}.{query_targets}.log",
        ),
    container:
        r
    script:
        "src/find_overlapping_loci.R"


# group loci in the peakall file
rule orthofinder:
    input:
        Path(
            outdir,
            "000_reference",
            "query",
            "{query_dataset}.grouped_targets",
        ),
    output:
        outdir=directory(
            Path(
                outdir,
                "005_grouped-targets",
                "{query_dataset}",
                "orthofinder",
            )
        ),
    log:
        Path(logdir, "orthofinder.{query_dataset}.log"),
    resources:
        mem_mb=int(32e3),
        time=30,
    threads: 12
    container:
        orthofinder
    shadow:
        "minimal"
    shell:
        "orthofinder "
        "-t {threads} "
        "-a {threads} "
        "-d "
        "-M msa "
        "-S diamond_ultra_sens "
        "-oa "
        "-o OrthoFinder "
        "-f {input} "
        "&> {log} "
        "&& "
        "find OrthoFinder -mindepth 1 -maxdepth 1 -type d "
        "-name 'Results_*' "
        "-exec mv {{}} {output.outdir} \;"


##########
# CAPTUS #
##########


rule captus_extract:
    wildcard_constraints:
        ref_dataset="|".join(all_reference_genomes),
        query_dataset="|".join(all_query_datasets),
    input:
        external_fasta=Path(
            outdir,
            "000_reference",
            "reference",
            "{ref_dataset}.{minlength}.fasta",
        ),
        targets=Path(
            outdir,
            "000_reference",
            "query",
            "{query_dataset}.deduplicated_renamed.fasta",
        ),
    output:
        outdir=directory(
            Path(
                outdir,
                "010_captus",
                "{ref_dataset}.{query_dataset}",
                "min{minlength}",
                "03_extractions",
            )
        ),
        refs_json=Path(
            outdir,
            "010_captus",
            "{ref_dataset}.{query_dataset}",
            "min{minlength}",
            "captus-assembly_extract.refs.json",
        ),
        annotated_assembly=Path(
            outdir,
            "010_captus",
            "{ref_dataset}.{query_dataset}",
            "min{minlength}",
            "03_extractions",
            "{ref_dataset}.{minlength}__captus-ext",
            "06_assembly_annotated",
            "{ref_dataset}.{minlength}_hit_contigs.gff",
        ),
    log:
        Path(
            logdir, "extract", "{ref_dataset}.{query_dataset}.{minlength}.log"
        ),
    benchmark:
        Path(
            logdir,
            "benchmark.extract.{ref_dataset}.{query_dataset}.{minlength}.log",
        )
    threads: lambda wildcards, attempt: 3  # bc BLAT is single-threaded
    resources:
        time=lambda wildcards, attempt: "7-00",
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
        "--nuc_refs {input.targets} "
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
            ref_dataset=all_reference_genomes,
            query_dataset=all_query_datasets,
            minlength=["1000000", "100000"],
        ),


############
# REF DATA #
############


# do the same thing for all reference genomes
rule index_ref_dataset:
    input:
        Path(outdir, "000_reference", "reference", "{ref_dataset}.fasta.gz"),
    output:
        Path(outdir, "000_reference", "reference", "{ref_dataset}.fasta.fai"),
    log:
        Path(logdir, "index_ref_dataset.{ref_dataset}.log"),
    shadow:
        "minimal"
    container:
        samtools
    shell:
        "samtools faidx "
        "--fai-idx {output} "
        "<( gunzip -c {input} ) "
        "2> {log}"


rule reformat:
    wildcard_constraints:
        ref_dataset="|".join(all_reference_genomes),
    input:
        Path(outdir, "000_reference", "reference", "{ref_dataset}.fasta.gz"),
    output:
        Path(
            outdir,
            "000_reference",
            "reference",
            "{ref_dataset}.{minlength}.fasta",
        ),
    log:
        Path(logdir, "reformat.{ref_dataset}.{minlength}.log"),
    container:
        bbmap
    shell:
        "reformat.sh "
        "in={input} "
        "minlength={wildcards.minlength} "
        "out={output} "
        "2>{log}"


# this only needs to be done for the Peakall target file
rule group_targets_by_prefix:
    input:
        targets=Path(
            outdir,
            "000_reference",
            "query",
            "{query_dataset}.deduplicated_renamed.fasta",
        ),
    output:
        outdir=directory(
            Path(
                outdir,
                "000_reference",
                "query",
                "{query_dataset}.grouped_targets",
            )
        ),
    container:
        biopython
    script:
        "src/group_targets_by_prefix.py"


rule rename_target_sequences:
    wildcard_constraints:
        query_dataset="|".join(all_query_datasets),
    input:
        Path(
            outdir,
            "000_reference",
            "query",
            "{query_dataset}.deduplicated.fasta",
        ),
    output:
        Path(
            outdir,
            "000_reference",
            "query",
            "{query_dataset}.deduplicated_renamed.fasta",
        ),
    log:
        Path(logdir, "rename_target_sequences.{query_dataset}.log"),
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


rule dedupe_targetfile:
    input:
        Path(outdir, "000_reference", "query", "{query_dataset}.fasta.gz"),
    output:
        pipe=pipe(
            Path(
                outdir,
                "000_reference",
                "query",
                "{query_dataset}.deduplicated.fasta",
            )
        ),
        duplicates=Path(
            outdir,
            "000_reference",
            "query",
            "{query_dataset}.discarded_duplicates.fasta",
        ),
    log:
        Path(logdir, "dedupe_targetfile.{query_dataset}.log"),
    threads: 2
    resources:
        time=1,
    container:
        bbmap
    shell:
        "dedupe.sh "
        "in={input} "
        "out=stdout.fasta "
        "outd={output.duplicates} "
        "ascending=t "
        "exact=t "
        "fixjunk=t "
        "maxedits=0 "
        "maxsubs=0 "
        "sort=name "
        "touppercase=t "
        "uniquenames=t "
        "> {output.pipe} "
        "2> {log} "


rule collect_remote_file:
    input:
        get_remote_file,
    output:
        Path(outdir, "000_reference", "{datatype}", "{dataset}.fasta.gz"),
    params:
        conversion_code=generate_conversion_code,
    threads: 1
    resources:
        partition_flag="--partition=io",
    shell:
        "{params.conversion_code}"


###########
# TARGETS #
###########


rule target:
    default_target: True
    input:
        rules.find_overlaps_target.input,
