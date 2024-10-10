#!/usr/bin/env Rscript

if (exists("snakemake")) {
  log <- file(snakemake@log[[1]], open = "wt")
  sink(log, type = "message")
  sink(log, append = TRUE, type = "output")

  # inputs
  overlapping_loci_file <- snakemake@input[["overlapping_loci"]]
  ref_self_overlaps_file <- snakemake@input[["ref_self_overlaps"]]
  query_self_overlaps_file <- snakemake@input[["query_self_overlaps"]]

  # outputs
  outdir <- snakemake@output[["outdir"]]
  # params
} else {
  overlapping_loci_file <-
    "output/020_overlaps/pzijinensis_min1000000.mega353.peakall/overlapping_loci.csv"
  ref_self_overlaps_file <-
    "output/020_overlaps/pzijinensis_min1000000.mega353.mega353/overlapping_loci.csv"
  query_self_overlaps_file <-
    "output/020_overlaps/pzijinensis_min1000000.peakall.peakall/overlapping_loci.csv"
  outdir <- "test/loci_to_merge"
}

check_if_paralog <- function(hit_group_name) {
  if (!grepl("__", hit_group_name)) {
    return(FALSE)
  }
  hit_split <- unlist(strsplit(hit_group_name, "__"))
  if (as.integer(hit_split[[length(hit_split)]]) != 0) {
    return(TRUE)
  }
  return(FALSE)
}

remove_paralogs_from_dt <- function(overlap_dt) {
  my_dt <- copy(overlap_dt)
  mycols <- names(overlap_dt)
  my_dt[
    , ref_hit_is_paralog := check_if_paralog(ref_hit_group_name),
    by = ref_hit_group_name
  ]
  my_dt[
    , query_hit_is_paralog := check_if_paralog(query_hit_group_name),
    by = query_hit_group_name
  ]
  return(
    my_dt[
      !(ref_hit_is_paralog)
    ][
      !(query_hit_is_paralog)
    ][, ..mycols]
  )
}

library(data.table)

# Check for self-overlaps in reference loci first. This is NOT HANDLED in this
# pipeline so just quit if we find any.
ref_self_overlaps <- fread(ref_self_overlaps_file)
ref_no_paralogs <- remove_paralogs_from_dt(ref_self_overlaps)

if (ref_no_paralogs[, any(ref_hit_group_name != query_hit_group_name)]) {
  stop("Giving up. There seem to be overlapping reference loci.")
}

# are any of the separate query loci actually the same
query_self_overlaps <- fread(query_self_overlaps_file)
query_no_paralogs <- remove_paralogs_from_dt(query_self_overlaps)
query_no_paralogs[ref_hit_group_name != query_hit_group_name]

# Now list loci that overlap each reference loci
overlapping_loci <- fread(overlapping_loci_file)
paralogs_exluded <- remove_paralogs_from_dt(overlapping_loci)

if (paralogs_exluded[, any(ref_hit_group_name == query_hit_group_name)]) {
  stop("Giving up. Some of the reference and query loci share names.")
}

query_loci_per_ref_locus <- paralogs_exluded[
  , .(
    query_hit = unique(query_hit_group_name),
    seqname = unique(query_hit_seqnames)
  ),
  by = ref_hit_group_name
]

# drop paralog numbers from the ref_hit, because these aren't in the target
# file.
query_loci_per_ref_locus[
  , ref_hit_group_name := gsub("__.*", "", ref_hit_group_name),
  by = ref_hit_group_name
]

setnames(query_loci_per_ref_locus, "ref_hit_group_name", "ref_hit")

query_locus_list <- split(
  query_loci_per_ref_locus, query_loci_per_ref_locus$ref_hit
)

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

lapply(names(query_locus_list), function(x) {
  fwrite(query_locus_list[[x]], paste0(outdir, "/", x, ".csv"))
})

sessionInfo()
