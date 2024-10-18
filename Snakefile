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


def get_remote_script(wildcards):
    remote_script = remote_scripts[wildcards.script_name]
    if isinstance(remote_script, snakemake.sourcecache.GithubFile):
        return HTTP.remote(remote_script.get_path_or_uri(), keep_local=True)
    else:
        raise ValueError(f"Implement handler for {type(remote_script)}")


# containers
bbmap = "docker://quay.io/biocontainers/bbmap:39.01--h92535d8_1"
biopython = "docker://quay.io/biocontainers/biopython:1.78"
captus = "docker://quay.io/biocontainers/captus:1.0.1--pyhdfd78af_2"
orthofinder = "docker://davidemms/orthofinder:2.5.5.2"
r = "docker://ghcr.io/tomharrop/r-containers:r2u_24.04_cv1"
samtools = "docker://quay.io/biocontainers/samtools:1.21--h50ea8bc_0"
qcat = "docker://quay.io/biocontainers/qcat:1.1.0--py_0"

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
    "peakall12": Path("data", "reference", "S1-2FinalSeq.fa"),
    "peakall35": Path("data", "reference", "S3-5FinalSeq.fa"),
    # Downloaded from
    # https://datadryad.org/stash/dataset/doi:10.5061/dryad.sj3tx96bn
    "orchidinae205": Path(
        "data", "reference", "orchidinae-205_reference.fasta"
    ),
}

all_query_datasets = query_file_locations.keys()


remote_scripts = {
    "filter_mega353.py": github(
        repo="chrisjackson-pellicle/NewTargets",
        path="filter_mega353.py",
        tag="v1.0.0",
    )
}


wildcard_constraints:
    ref_targets="|".join(all_query_datasets),
    query_targets="|".join(all_query_datasets),
    ref_dataset="|".join(all_reference_genomes),
    query_dataset="|".join(all_query_datasets),


#########################################################
# RENAME THE ORIGINAL TARGETFILES TO MATCH MERGING INFO #
#########################################################


rule rename_original_targetfile:
    input:
        targetfile=Path(
            outdir,
            "000_reference",
            "query",
            "{query_dataset}.deduplicated_renamed.fasta",
        ),
        orthofinder=Path(
            outdir,
            "005_grouped-targets",
            "{query_dataset}",
            "orthofinder",
        ),
        loci_to_merge=expand(
            Path(
                outdir,
                "020_overlaps",
                "{ref_dataset}_min{minlength}.{ref_targets}.{{query_dataset}}",
                "loci_to_merge",
            ),
            ref_dataset=["qos", "pzijinensis"],
            minlength=["1000000"],
            ref_targets=["mega353"],
        ),
    output:
        renamed_targetfile=Path(
            outdir,
            "050_updated-targetfiles",
            "{query_dataset}",
            "renamed_targetfile.fasta",
        ),
        renamed_sequences=Path(
            outdir,
            "050_updated-targetfiles",
            "{query_dataset}",
            "renamed_targets.csv",
        ),
    log:
        Path(logdir, "rename_original_targetfile", "{query_dataset}.log"),
    container:
        biopython
    script:
        "src/rename_original_targetfile.py"




#####################################
# RE-RUN CAPTUS WITH MERGED TARGETS #
#####################################


rule captus_extract_round2:
    input:
        external_fasta=Path(
            outdir,
            "000_reference",
            "reference",
            "{ref_dataset}.{minlength}.fasta",
        ),
        targetfile_list=expand(
            Path(
                outdir,
                "030_merged-target-sequences",
                "{ref_dataset}_min{minlength}.{ref_targets}.{query_targets}",
                "merged_targets.fasta",
            ),
            ref_dataset=["qos", "pzijinensis"],
            minlength=["1000000"],
            ref_targets=["mega353"],
            query_targets=["peakall12", "peakall35"],
        ),
    output:
        outdir=directory(
            Path(
                outdir,
                "040_captus_round2",
                "{ref_dataset}",
                "min{minlength}",
                "03_extractions",
            )
        ),
        refs_json=Path(
            outdir,
            "040_captus_round2",
            "{ref_dataset}",
            "min{minlength}",
            "captus-assembly_extract.refs.json",
        ),
    log:
        Path(logdir, "extract_round2", "{ref_dataset}.{minlength}.log"),
    benchmark:
        Path(
            logdir,
            "benchmark.extract.{ref_dataset}.{minlength}.log",
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
        "cat {input.targetfile_list} > my_targets.fa && "
        "captus_assembly extract "
        "--captus_assemblies_dir 02_assemblies "
        "--fastas {input.external_fasta} "
        "--out {output.outdir}/. "
        "--nuc_refs my_targets.fa "
        "--mit_refs SeedPlantsMIT "
        "--ptd_refs SeedPlantsPTD "
        '--ram "$(( {resources.mem_mb}/1000 ))" '
        "--threads {threads} "
        "--concurrent 3 "
        "&> {log} "
        "; mv {output.outdir}/captus-assembly_extract.refs.json "
        "{output.refs_json}"


#########################
# FIND OVERLAPPING LOCI #
#########################


rule find_overlaps_target:
    input:
        # includes check for overlaps among reference loci
        expand(
            Path(
                outdir,
                "030_merged-target-sequences",
                "{ref_dataset}_min{minlength}.{ref_targets}.{query_targets}",
                "merged_targets.no_captus_paralogs.fasta",
            ),
            ref_dataset=all_reference_genomes,
            minlength=["1000000"],
            ref_targets=["mega353"],
            query_targets=[x for x in all_query_datasets if x != "mega353"],
        ),


# only implemented for peakall vs mega353
rule remove_captus_paralogs:
    input:
        targets=Path(
            outdir,
            "030_merged-target-sequences",
            "{ref_dataset}_min{minlength}.{ref_targets}.{query_targets}",
            "merged_targets.fasta",
        ),
    output:
        targets=Path(
            outdir,
            "030_merged-target-sequences",
            "{ref_dataset}_min{minlength}.{ref_targets}.{query_targets}",
            "merged_targets.no_captus_paralogs.fasta",
        ),
    log:
        Path(
            logdir,
            "remove_captus_paralogs",
            "{ref_dataset}_min{minlength}.{ref_targets}.{query_targets}.log",
        ),
    container:
        biopython
    script:
        "src/remove_captus_paralogs.py"


rule merge_extracted_sequences:
    input:
        query_captus_directory=Path(
            outdir,
            "010_captus",
            "{ref_dataset}.{query_targets}",
            "min{minlength}",
            "03_extractions",
        ),
        ref_captus_directory=Path(
            outdir,
            "010_captus",
            "{ref_dataset}.{ref_targets}",
            "min{minlength}",
            "03_extractions",
        ),
        loci_to_merge=expand(
            Path(
                outdir,
                "020_overlaps",
                "{ref_dataset}_min{minlength}.{ref_targets}.{{query_targets}}",
                "loci_to_merge",
            ),
            # hardcoding for this particular pipeline
            # TODO: generalise
            ref_dataset=["qos", "pzijinensis"],
            minlength=["1000000"],
            ref_targets=["mega353"],
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
    params:
        ref_targets=lambda wildcards, input: Path(
            input.ref_captus_directory,
            f"{wildcards.ref_dataset}.{wildcards.minlength}__captus-ext",
            "01_coding_NUC",
            "NUC_coding_NT.fna",
        ),
        query_targets=lambda wildcards, input: Path(
            input.query_captus_directory,
            f"{wildcards.ref_dataset}.{wildcards.minlength}__captus-ext",
            "01_coding_NUC",
            "NUC_coding_NT.fna",
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
        query_self_overlaps=Path(
            outdir,
            "020_overlaps",
            "{ref_dataset}_min{minlength}.{query_targets}.{query_targets}",
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
        query_captus_directory=Path(
            outdir,
            "010_captus",
            "{ref_dataset}.{query_targets}",
            "min{minlength}",
            "03_extractions",
        ),
        ref_captus_directory=Path(
            outdir,
            "010_captus",
            "{ref_dataset}.{ref_targets}",
            "min{minlength}",
            "03_extractions",
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
        query_gff=lambda wildcards, input: Path(
            input.query_captus_directory,
            f"{wildcards.ref_dataset}.{wildcards.minlength}__captus-ext",
            "06_assembly_annotated",
            f"{wildcards.ref_dataset}.{wildcards.minlength}_hit_contigs.gff",
        ),
        ref_gff=lambda wildcards, input: Path(
            input.ref_captus_directory,
            f"{wildcards.ref_dataset}.{wildcards.minlength}__captus-ext",
            "06_assembly_annotated",
            f"{wildcards.ref_dataset}.{wildcards.minlength}_hit_contigs.gff",
        ),
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


rule remove_mega_from_selected_families:
    input:
        targets=Path(
            outdir,
            "000_reference",
            "query",
            "{query_dataset}.selected_families_withmega.fasta",
        ),
        select_file=Path("data", "select_file.txt"),
    output:
        targets=Path(
            outdir,
            "000_reference",
            "query",
            "{query_dataset}.selected_families.fasta",
        ),
    log:
        Path(logdir, "remove_mega_from_selected_families.{query_dataset}.log"),
    container:
        biopython
    script:
        "src/remove_mega_from_selected_families.py"


rule select_families:
    input:
        targets=Path(
            outdir,
            "000_reference",
            "query",
            "{query_dataset}.deduplicated_renamed.fasta",
        ),
        select_file=Path("data", "select_file.txt"),
        script=Path(
            outdir, "000_reference", "external_scripts", "filter_mega353.py"
        ),
    output:
        targets=Path(
            outdir,
            "000_reference",
            "query",
            "{query_dataset}.selected_families_withmega.fasta",
        ),
        report=Path(logdir, "select_families.{query_dataset}.report.txt"),
    log:
        Path(logdir, "select_families.{query_dataset}.log"),
    shadow:
        "minimal"
    container:
        qcat  # has biopython and pandas
    shell:
        "python3 {input.script} "
        "-filtered_target_file {output.targets} "
        "-report_filename {output.report} "
        "{input.targets} "
        "{input.select_file} "
        "&> {log}"


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


##########
# UTIILS #
##########


rule download_remote_script:
    input:
        get_remote_script,
    output:
        Path(outdir, "000_reference", "external_scripts", "{script_name}"),
    shell:
        "cp {input} {output}"


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
        Path(
            outdir,
            "000_reference",
            "query",
            "mega353.selected_families.fasta",
        ),
        expand(
            Path(
                outdir,
                "040_captus_round2",
                "{ref_dataset}",
                "min{minlength}",
                "03_extractions",
            ),
            ref_dataset=["qos", "pzijinensis"],
            minlength=["1000000"],
        ),
        # updated original targetfiles
        expand(
            rules.rename_original_targetfile.output,
            query_dataset=["peakall12", "peakall35"],
        ),
